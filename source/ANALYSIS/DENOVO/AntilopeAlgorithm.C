/***************************************************************************
 *   Copyright (C) 2005 by Gunnar W. Klau, Markus Bauer, Patrick May       *
 *   gunnar@mi.fu-berlin.de, mbauer@inf.fu-berlin.de, patrick.may@zib.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

///This is an implementation of the k-shortest path algorithm form Jin Yen (1971)
//This is acually not one of the fastest k-shortest path algortihms but we can use it for our lagrange relaxation of
//the antisymmetric path problem

#include <OpenMS/ANALYSIS/DENOVO/AntilopeAlgorithm.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeLagrangeProblem.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeBranchBound.h>
#include <cmath>

#define FASTSUBOPT

#define FASTDEV
#undef DEBUG
//#define DEBUG

//#define MATRIX
namespace OpenMS
{  
  void YenAlgorithm::singleIteration_(Size k,
                                              const PathSolutionList& found_paths,
                                              PathSolution &path)
  {
    DeNovoLagrangeProblemBoost de_novo_lagrange(g_);
    if (k == 0)
    {
      //compute forward backward scores and backward traceback for efficient suboptimal solutions
      std::cout << "numVertices: " << seqan::numVertices(g_->graph) << std::endl;
#ifdef FASTSUBOPT
      seqan::String<DoubleReal> weights, dist_forward;
      seqan::String<VertexDescriptor> pred_forward;
      g_->getEdgeWeights(weights);
      seqan::String<VertexDescriptor> order;
      seqan::resizeVertexMap(g_->graph, order);
      for (Size i = 0; i < seqan::numVertices(g_->graph); ++i)
        order[i] = i;
      
      VertexDescriptor s = 0u;
      VertexDescriptor t = seqan::back(order);
      
      SpectrumGraphSeqan::dagShortestPathST(g_->graph, s, t, weights, pred_forward, dist_forward, order);
      
      seqan::transpose(g_->graph, graph_rev);
      
      seqan::reverse(order);
      SpectrumGraphSeqan::dagShortestPathST(graph_rev, t, s, weights, pred_backward, dist_backward, order);
#endif
      
      PathSolution cand_path;
      cand_path.score = 0;
      DoubleReal global_lower_bound = DeNovoLagrangeProblemBoost::MINUS_INF;

      BranchBoundDeNovo branch_bound(&global_lower_bound, de_novo_lagrange, &cand_path);
      branch_bound.solve();

      cand_path.deviation_node = 0;
      cand_path.parent_id = -1u;
      cand_heap_.insert(cand_path);
    }
    else
    {
      seqan::String<DoubleReal> weights;
      g_->getEdgeWeights(weights);
      
#ifdef FASTDEV
      //forbid edges and identify deviation node id
      Size deviation_node = found_paths[k-1].deviation_node;
      Size parent_id = k-1;
      while (found_paths[parent_id].parent_id != -1u && deviation_node <= found_paths[parent_id].deviation_node)
      {
        parent_id = found_paths[parent_id].parent_id;
        de_novo_lagrange.forbidEdge(found_paths[parent_id].path[deviation_node], found_paths[parent_id].path[deviation_node+1]);
      }
#else
      Size deviation_node = forbidEdges_(k, found_paths, de_novo_lagrange);
#endif
      
      const PathSolution::Path& parent_path = found_paths[k-1].path;
      std::vector<std::pair<HeuristicReturn, PathSolution> > candidates(parent_path.size());
      DoubleReal pre_score = 0;
      for(Size i = 0; i+1 < found_paths[k-1].path.size(); ++i)
      {
        candidates[i].second.score = pre_score;
        pre_score -= g_->getEdgeWeight(parent_path[i], parent_path[i+1]);
        candidates[i].first = INFEASIBLE;
      }

      PathSolution cand_path, prefix;
      
      cand_path.score = 0;
      prefix.score = 0;
      prefix.path.push_back(0);

      for(Size i = 1; i <= deviation_node;  ++i)
      {
        de_novo_lagrange.forbidConflictingNodes(parent_path[i]);
        prefix.path.push_back(parent_path[i]);
        prefix.score += g_->getEdgeWeight(parent_path[i-1], parent_path[i]);
      }
      
      for (Size i = deviation_node; i+1 < parent_path.size(); ++i)
      {
        cand_path.score = DeNovoLagrangeProblemBoost::MINUS_INF;
        cand_path.path.clear();
#ifdef DEBUG
        cout << "deviation node: " << i << endl;
#endif
        
//#ifdef MATRIX
//        if (parent_path[i]%2)
//        {
//          de_novo_lagrange.setPrefix(prefix);
//          de_novo_lagrange.forbidConflictingNodes(prefix.path.back());
//
//          prefix.path.push_back(parent_path[i+1]);
//          prefix.score += g_->getEdgeWeight(parent_path[i], parent_path[i+1]);
//          std::cout << "SKIPP " << parent_path[i] << std::endl;
//          continue;
//        }
//#endif
      
        de_novo_lagrange.setPrefix(prefix);
        de_novo_lagrange.forbidConflictingNodes(prefix.path.back());
        de_novo_lagrange.forbidEdge(parent_path[i], parent_path[i+1]);
        
        de_novo_lagrange.resetBestFeasible();
        de_novo_lagrange.resetBestInfeasible();
        
        DoubleReal lower_bound = cand_heap_.rbegin()->score;
        de_novo_lagrange.setLowerBound(lower_bound);
        
        cand_path.path = prefix.path;
        cand_path.score = -prefix.score;

        
#ifdef FASTSUBOPT
        DoubleReal timeFast = clock();
//        //get maximal weighted candidate successor using backward information
//        SpectrumGraphSeqan::OutEdgeIterator out_it(g_->graph, parent_path[i]);
//        VertexDescriptor succ = -1u;
//        DoubleReal back_score = SpectrumGraphSeqan::INFINITYdist / 2;
//        while (!seqan::atEnd(out_it))
//        {
//          VertexDescriptor target = seqan::targetVertex(g_->graph, *out_it);
//          if (!de_novo_lagrange.isEdgeForbidden(*out_it) && !de_novo_lagrange.isVertexForbidden(target))
//          {
//            std::cout << seqan::property(dist_backward, target) << std::endl;
//            if (seqan::property(dist_backward, target) + g_->getEdgeWeight(*out_it) < back_score)
//            {
//              back_score = seqan::property(dist_backward, target) + g_->getEdgeWeight(*out_it);
//              succ = target;
//            }
//          }
//          if (de_novo_lagrange.isEdgeForbidden(*out_it))
//            std::cout << "found forbidden: " << seqan::sourceVertex(out_it) << " -> " << seqan::targetVertex(out_it) << std::endl;
//          seqan::goNext(out_it);
//        }
//
//        if (succ != -1u && prefix.score + back_score > lower_bound)
//        {
//          cand_path.path = prefix.path;
//          cand_path.score = prefix.score + back_score;
//
//          //backtrace suffix
//          VertexDescriptor nextV = succ;
//          while (true)
//          {
//            cand_path.path.push_back(nextV);
//            if (nextV == seqan::numVertices(g_->graph)-1)
//              break;
//            
//            nextV = seqan::property(pred_backward, nextV);
//          }
//        
//          bool feasible = g_->isPathFeasible(cand_path.path);
//          std::cout << "Fast: " << (clock() - timeFast)/CLOCKS_PER_SEC << std::endl;
//          if (!feasible)
//          {
//            cand_path.path.clear();
//            cand_path.score = LagrangeProblem::MINUS_INF;
        if (candidates[i].first == FEASIBLE)
        {
          cand_path = candidates[i].second;
//          std::cout << "1" << std::endl;
//          for (Size tt = 0; tt < candidates[i].second.path.size(); ++tt)
//            std::cout << candidates[i].second.path[tt] << '\t';
//          std::cout << std::endl;
        }
        else if (candidates[i].first == INFEASIBLE) //not previously solved for deviation node i
        {
          HeuristicReturn ret = backWardHeuristic1(de_novo_lagrange, dist_backward, pred_backward, lower_bound, cand_path);
          if (ret == INFEASIBLE)
          {
            backWardHeuristic2(de_novo_lagrange, weights, graph_rev, lower_bound, i, candidates, parent_path);

            if (candidates[i].first == FEASIBLE)
            {
              cand_path = candidates[i].second;
//              std::cout << "2" << std::endl;
//              for (Size tt = 0; tt < candidates[i].second.path.size(); ++tt)
//                std::cout << candidates[i].second.path[tt] << '\t';
//              std::cout << std::endl;

            }
            else if (candidates[i].first == INFEASIBLE)//still not solved for deviation node i?
            {              
#endif
              cand_path.score = LagrangeProblem::MINUS_INF;
//              DoubleReal timeBB = clock();
              DoubleReal global_lower_bound = DeNovoLagrangeProblemBoost::MINUS_INF;
              BranchBoundDeNovo branch_bound(&global_lower_bound, de_novo_lagrange, &cand_path);
              branch_bound.solve();
//              std::cout << "BB: " << (clock() - timeBB)/CLOCKS_PER_SEC << std::endl;
#ifdef FASTSUBOPT
            }
          }
        }
        else
        {
          cand_path.score = LagrangeProblem::MINUS_INF;
        }
//        else
//        {
//          std::cout << "found path fast " << cand_path.score << std::endl;
//        }
//      }
//      else
//      {
//        std::cout << "bounding fast " << std::endl;
//      }
#endif

        //if path has better score than lowest score in heap, the replace this lowest heap element with it
        if (cand_path.score > lower_bound)
        {
          //std::cout<<"score vgl: "<<score<<" :: "<<qu_score<<std::endl;
          cand_path.parent_id = k-1;
          cand_path.deviation_node = i;

          cand_heap_.insert(cand_path);
          cand_heap_.erase(--cand_heap_.end());

//          std::cout<<"added candidate path " << i << ": " << cand_path.score << std::endl;
//          for (Size tt = 0; tt < cand_path.path.size(); ++tt)
//            std::cout << cand_path.path[tt] << '\t';

        }
        prefix.path.push_back(parent_path[i+1]);
        prefix.score += g_->getEdgeWeight(parent_path[i], parent_path[i+1]);
      }
    }    

    //std::cerr<<"TOP CAND SCORE: "<< cand_path_heap.begin()->score <<std::endl;
    //std::cerr<<"CAND SIZE: "<<cand_path_heap.size()<<std::endl;

    path = *(cand_heap_.begin());
    cand_heap_.erase(cand_heap_.begin());
  }



  //forbid edges as defined in yen paper page 714, step a
  Size YenAlgorithm::forbidEdges_(Size k, const PathSolutionList& found_paths, DeNovoLagrangeProblemBoost& dnlp)
  {
    if (k == 0)
    {
      return 0;
    }
    else
    {
      Size deviation_node = 0;
      const vector<VertexDescriptor> &last_path = found_paths[k - 1].path;

      //forbid edges from paths in found_path
      for (Size i = 0; i < k-1; ++i)
      {
        Size count = 0;
        vector<VertexDescriptor>::const_iterator last_path_it = last_path.begin();
        vector<VertexDescriptor>::const_iterator found_path_it = found_paths[i].path.begin();
          
        std::cerr << found_paths.size() << " " << last_path.size() << " " << found_paths[i].path.size() << std::endl;

        while(*last_path_it == *found_path_it)
        {
          ++last_path_it;
          ++found_path_it;
          ++count;
        }
        dnlp.forbidEdge(*(--last_path_it), *found_path_it);
#ifdef DEBUG
        cout << "forbid edge: " << *last_path_it << "::" << *found_path_it << endl;
#endif
        if (count)
          deviation_node = max(deviation_node, count - 1);
      }
      return deviation_node;
    }
  }



  int YenAlgorithm::computeLongestPaths(Size k, PathList& paths, PathScoreList& scores)
  {
    cand_heap_.clear();
    PathSolutionList found_paths;
    PathSolution tmp = { vector<VertexDescriptor>(0), -SpectrumGraphSeqan::INFINITYdist, -1u, -1u};
    std::vector<PathSolution>tmp_v(k, tmp);
    cand_heap_.insert(tmp_v.begin(), tmp_v.end());

    double tstart = clock();
    double time_total = 0;

    for (Size i = 0; i < k; i++)
    {
//      cout << "start for " << i << "-th path" << endl;
      PathSolution longest_path_tmp;
      singleIteration_(i, found_paths, longest_path_tmp);
      
      found_paths.push_back(longest_path_tmp);
      paths.push_back(longest_path_tmp.path);
      scores.push_back(longest_path_tmp.score);

      //found_path_scores.push_back(longest_path_tmp.score);
      //found_path.push_back(longest_path_tmp.path);
      cout << "found path number: " << i <<"  with score  " << found_paths.back().score << endl;
#ifdef DEBUG
      for(UInt tt = 0; tt < found_paths[i].path.size();++tt)
      {
        std::cout<<" - "<< found_paths[i].path[tt] << std::endl;
      }
#endif
      double time = (clock() - tstart) / CLOCKS_PER_SEC;
      time_total = time;
//      cout << "TIME_LAGRANGE::" << i << "::" << time << "::" << time_total << endl;
    }

    cout << " TOTAL TIME_LAGRANGE:: " <<time_total << endl;

    for (Size i = 0; i < k; i++)
    {
      vector<VertexDescriptor>::iterator found_path_it = found_paths[i].path.begin();
      for (;found_path_it != found_paths[i].path.end(); ++found_path_it)
      {
        cout << *found_path_it << endl;
      }
      cout << found_paths[i].score << endl;
    }
    return 0;
  }

  
  YenAlgorithm::HeuristicReturn YenAlgorithm::backWardHeuristic1(
                          const DeNovoLagrangeProblemBoost& de_novo_lagrange,
                          const seqan::String<DoubleReal>& dist_backward,
                          const seqan::String<VertexDescriptor>& pred_backward,
                          DoubleReal lower_bound,
                          PathSolution& cand_path                          
                          )
  {
    DoubleReal timeFast = clock();
    //get maximal weighted candidate successor using backward information
    SpectrumGraphSeqan::OutEdgeIterator out_it(g_->graph, cand_path.path.back());
    VertexDescriptor succ = -1u;
    DoubleReal back_score = SpectrumGraphSeqan::INFINITYdist / 2;
    while (!seqan::atEnd(out_it))
    {
      VertexDescriptor target = seqan::targetVertex(g_->graph, *out_it);
      if (!de_novo_lagrange.isEdgeForbidden(*out_it) && !de_novo_lagrange.isVertexForbidden(target))
      {
//        std::cout << seqan::property(dist_backward, target) << std::endl;
        if (seqan::property(dist_backward, target) + g_->getEdgeWeight(*out_it) < back_score)
        {
          back_score = seqan::property(dist_backward, target) + g_->getEdgeWeight(*out_it);
          succ = target;
        }
      }
//      if (de_novo_lagrange.isEdgeForbidden(*out_it))
//        std::cout << "found forbidden: " << seqan::sourceVertex(out_it) << " -> " << seqan::targetVertex(out_it) << std::endl;
      seqan::goNext(out_it);
    }
    
    if (succ == -1u)
    {
      cand_path.path.clear();
      cand_path.score = LagrangeProblem::MINUS_INF;
      return NOPATH;
    }

//    std::cout << "Pref: " << cand_path.score << " back score: " << back_score << std::endl;
    cand_path.score -= back_score;
    if (cand_path.score <= lower_bound)
      return PRUNED;
    
    //backtrace suffix
    VertexDescriptor nextV = succ;
    while (true)
    {
      cand_path.path.push_back(nextV);
      if (nextV == seqan::numVertices(g_->graph)-1)
        break;
      
      nextV = seqan::property(pred_backward, nextV);
    }
    
    bool feasible = g_->isPathFeasible(cand_path.path);
//    std::cout << "Fast: " << (clock() - timeFast)/CLOCKS_PER_SEC << std::endl;
    if (!feasible)
    {
      cand_path.path.clear();
      cand_path.score = LagrangeProblem::MINUS_INF;
      return INFEASIBLE;
    }
    return FEASIBLE;
  }
  
  void YenAlgorithm::backWardHeuristic2(const DeNovoLagrangeProblemBoost& de_novo_lagrange,
                                                                 const seqan::String<DoubleReal>& weights,
                                                                 SpectrumGraphSeqan::TGraph& rev_graph,
                                                                 DoubleReal lower_bound,
                                                                 Size s,
                                                                 std::vector<std::pair<HeuristicReturn, PathSolution> >& candidates,
                                                                 const std::vector<VertexDescriptor>& parent_path
                                                                 )
  {
    
    seqan::String<DoubleReal> dist_backward;
    seqan::String<VertexDescriptor> pred_backward;
    seqan::String<VertexDescriptor> order_rev = g_->getTopologicalOrdering();
    seqan::reverse(order_rev);
    
    VertexDescriptor source = order_rev[0];
    VertexDescriptor target = parent_path[s];
    
  
    
    SpectrumGraphSeqan::dagShortestPathST(rev_graph,
                                          source,
                                          target,
                                          weights,
                                          pred_backward,
                                          dist_backward,
                                          de_novo_lagrange.getForbiddenNodes(),
                                          de_novo_lagrange.getForbiddenEdges(),
                                          std::set<VertexDescriptor>(),
                                          order_rev);


    DoubleReal timeFast = clock();
    //get maximal weighted candidate successor using backward information
    
    for (Size i = s; i+1 < parent_path.size(); ++i)
    {
      if (candidates[i].first != INFEASIBLE)
        continue;
      
      SpectrumGraphSeqan::OutEdgeIterator out_it(g_->graph, parent_path[i]);
      VertexDescriptor succ = -1u;
      DoubleReal back_score = SpectrumGraphSeqan::INFINITYdist / 2;
      while (!seqan::atEnd(out_it))
      {
        VertexDescriptor target = seqan::targetVertex(g_->graph, *out_it);
        if (!de_novo_lagrange.isEdgeForbidden(*out_it) &&
            !de_novo_lagrange.isVertexForbidden(target) &&
            target != parent_path[i+1])
        {
//          std::cout << seqan::property(dist_backward, target) << std::endl;
          if (seqan::property(dist_backward, target) + g_->getEdgeWeight(*out_it) < back_score)
          {
            back_score = seqan::property(dist_backward, target) + g_->getEdgeWeight(*out_it);
            succ = target;
          }
        }
//        if (de_novo_lagrange.isEdgeForbidden(*out_it))
//          std::cout << "found forbidden: " << seqan::sourceVertex(out_it) << " -> " << seqan::targetVertex(out_it) << std::endl;
        seqan::goNext(out_it);
      }
      
      if (succ == -1u)
      {
        candidates[i].first = NOPATH;
        candidates[i].second.score = DeNovoLagrangeProblemBoost::MINUS_INF;
        //std::cout << "NOPATH " << i << std::endl;
      }
      else if (candidates[i].second.score - back_score <= lower_bound)
      {
        candidates[i].first = PRUNED;
        candidates[i].second.score = DeNovoLagrangeProblemBoost::MINUS_INF;
        //std::cout << "PRUNED " << i << std::endl;
      }
      else
      {
        candidates[i].second.path.assign(parent_path.begin(), parent_path.begin()+i+1);
        
        //backtrace suffix
        VertexDescriptor nextV = succ;
        while (true)
        {
          candidates[i].second.path.push_back(nextV);
          if (nextV == source)
            break;
          
          nextV = seqan::property(pred_backward, nextV);
        }
        
        bool feasible = g_->isPathFeasible(candidates[i].second.path);
        if (feasible)
        {
          candidates[i].first = FEASIBLE;
          candidates[i].second.score -= back_score;
          //std::cout << "FEASIBLE " << i << std::endl;
        }
      }
    }
//     std::cout << "Fast: " << (clock() - timeFast)/CLOCKS_PER_SEC << std::endl;
  }

}//namespace
