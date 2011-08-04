/***************************************************************************
 *   Copyright (C) 2005 by Gunnar W. Klau, Markus Bauer                    *
 *   gunnar@mi.fu-berlin.de, mbauer@inf.fu-berlin.de                       *
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

#include <OpenMS/ANALYSIS/DENOVO/AntilopeLagrangeProblem.h>

#include <iostream>
#include <vector>
#include <numeric>
#undef DEBUG
#define DEBUG
//using namespace std;

namespace OpenMS{

  const DoubleReal LagrangeProblem::PLUS_INF = std::numeric_limits<DoubleReal>::max();
  const DoubleReal LagrangeProblem::MINUS_INF = -1 * PLUS_INF;


  int DeNovoLagrangeProblemBoost::EvaluateProblem(
      const DVector& Dual,
      list<int>& DualIndices,
      double& DualValue,
      double& PrimalValue,
      DVector& Subgradient,
      list<int>& SubgradientIndices,
      DVector& PrimalSolution,
      DVector& PrimalFeasibleSolution
      )
  {
    computeLongestPath(DualValue, PrimalValue, Dual);
    if(DualValue == MINUS_INF)
    {
      PrimalValue = MINUS_INF;
    }
    else
    {
      computeSubgradientIndices(SubgradientIndices);
    }
    std::cout<<"Dual::Primal::LowerBound   "<<DualValue<<" :: "<<PrimalValue<<"  ::  "<<lower_bound_<<std::endl;

    if(DualValue < lower_bound_){
      PrimalValue = DualValue = MINUS_INF;
      std::cout<<"INTERRUPTED"<<std::endl;
    }
    Subgradient= subgradient_;
    return 0;
  }



  //-------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------DeNovo Specific Functions------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------

  //because this has to be done only once for each iteration in YenAlgorithm it was isolated from the path search to safe runtime
  void DeNovoLagrangeProblemBoost::forbidConflictingNodes(VertexDescriptor v)
  {
    if(v == 0)
    {
      return; //for the first node nothing needs to be done (no conflicts)
    }

    const std::vector<VertexDescriptor> & nodes = G->getConflictingNodes(v);

    for(Size i = 0; i < nodes.size(); ++i)
    {
      forbidNode(nodes[i]);
      std::cerr<<" forbid node: "<< nodes[i] <<std::endl;
    }


#ifdef DEBUG
    std::cout<<"the forbidden nodes are:"<<std::endl;

    seqan::Iterator<const seqan::String<bool> >::Type forb_it = seqan::begin(forbidden_nodes_), end_it = seqan::end(forbidden_nodes_);

    while(forb_it != end_it)
    {
      if(*forb_it)
      {
        std::cout<<"node " << forb_it - seqan::begin(forbidden_nodes_) << std::endl;
      }
      ++forb_it;
    }
#endif
  }

  ///compute the longest path from node start to node end
  void DeNovoLagrangeProblemBoost::computeLongestPath(DoubleReal &dual_score, DoubleReal &primal_score, const DVector& lambda)
  {
    std::cout << "start Pathfinding" << std::endl;
    dual_score = MINUS_INF, primal_score = MINUS_INF;

    if(prefix_.path.empty())
    {
      prefix_.path.push_back(0);
      prefix_.score = 0;
    }

    VertexDescriptor start_vertex = prefix_.path.back();

    updateWeights(lambda);

    seqan::String<VertexDescriptor> predecessors;
    seqan::String<DoubleReal> distance;

    //perform dag search
    seqan::String<VertexDescriptor> ordering;
    SpectrumGraphSeqan::VertexIterator it(G->graph);
    for(; !seqan::atEnd(it); ++it)
    {
      //TODO: remove this, as it is not necssary since order does never change. compute one central ordering in SpectrumGraph
      seqan::appendValue(ordering, *it);
    }

    seqan::dagShortestPathST(G->graph, start_vertex, seqan::back(ordering), edge_weights_, predecessors, distance, ordering, forbidden_nodes_, forbidden_edges_);

    if(seqan::back(distance) >= SpectrumGraphSeqan::INFINITYdist / 2)
    {
      std::cout << "no path found!!" << std::endl;
      return;
    }

    dual_score = seqan::back(distance);

    //backtracking and storing of the path
    std::cout << "start backtracking" << std::endl;
    std::vector<VertexDescriptor> path;
    path.reserve(30);
    subgradient_.assign(G->getNumberOfClusters(), 1);

    VertexDescriptor next_vertex = seqan::back(ordering);

    while(next_vertex != start_vertex)
    {
      path.push_back(next_vertex);
      next_vertex = seqan::getProperty(predecessors, next_vertex);
    }

    reverse(path.begin(), path.end());
    path.insert(path.begin(), prefix_.path.begin(), prefix_.path.end());


    //compute subgradients
    std::vector<Size>::const_iterator cluster_it;
    std::vector<VertexDescriptor>::const_iterator v_it = path.begin();
    for (; v_it != path.end(); ++v_it)
    {
      for(cluster_it = G->getClusters(*v_it).begin(); cluster_it != G->getClusters(*v_it).end(); ++cluster_it)
      {
        --subgradient_[*cluster_it];
      }
    }

    //Debug out
    std::cout<<"Found path:"<<std::endl;
    for(Size deb = 0; deb < path.size(); ++ deb)
    {
      std::cout<<path[deb]<<std::endl;
    }

    // add the lambda summand to the complete score
    DoubleReal lambda_update = std::accumulate(lambda.begin(), lambda.end(), 0.);

    std::cout<<"relaxed score: "<<dual_score<<std::endl;
    std::cout<<"prefix_score: "<< prefix_.score << "lambda update: "<<lambda_update<<std::endl;
    dual_score += prefix_.score - lambda_update;
    dual_score *= -1; //go from minimization back to maximization

    //If necessary adapt be primal and dual solution solution found until now
    if(!path.empty())
    {
      if(check_feasibility())
      {
        //compute the primal score of actual solution. Udate best primal solution if necessary.
        primal_score = 0;
        int prev = path[0];
        for (int i = 1; i < path.size(); ++i)
        {
          primal_score += G->getEdgeWeight(prev, path[i]);
          prev = path[i];
        }
        primal_score *= -1;
        std::cerr<<"dual score: "<<dual_score<<"  primal score: "<<primal_score<<std::endl;

        if (primal_score > best_feasible_solution_.score)
        {
          best_feasible_solution_.score = primal_score;
          best_feasible_solution_.path = path;
        }
      }
      else if(best_infeasible_solution_.score > dual_score)
      {
        best_infeasible_solution_.score = dual_score;
        best_infeasible_solution_.path = path;
      }
    }
  }


  //-------------------------------------------------------------------------------------------------------------------------
  //-----------------------------------------------------reset---------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------

  void DeNovoLagrangeProblemBoost::reset()
  {
    best_feasible_solution_.path.clear();
    best_feasible_solution_.score = MINUS_INF;
    best_infeasible_solution_.path.clear();
    best_infeasible_solution_.score = PLUS_INF;
  }



  //-------------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------check_feasibility----------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  //TODO remove this stuff. can be computed within the computeLongestPath function
  bool DeNovoLagrangeProblemBoost::check_feasibility()
  {
    bool feasible = true;
    for (int i = 0; i < subgradient_.size(); i++) {
      feasible = feasible && (subgradient_[i] > -1);
    }
    if (feasible)
      std::cout << "FEASIBLE!!" << std::endl;
    return feasible;
  }



  //-------------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------get_longest_path------------------------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------
  DeNovoLagrangeProblemBoost::path_score_pair DeNovoLagrangeProblemBoost::get_longest_path()
  {
    return best_feasible_solution_;
  }



  //-------------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------compute_subgradient_indices----------------------------------------------
  //-------------------------------------------------------------------------------------------------------------------------

  void DeNovoLagrangeProblemBoost::computeSubgradientIndices(std::list<int> &subgradient_ind)
  {
    subgradient_ind.clear();
    for (Size i = 0; i < G->getNumberOfClusters(); ++i)
    {
      if (subgradient_[i])
      {
        subgradient_ind.push_back(i);
      }
    }
  }



  void DeNovoLagrangeProblemBoost::forceNode(VertexDescriptor node)
  {
    for(Size i = 0; i < spanning_edges_[node].size(); ++i)
    {
      forbidEdge(spanning_edges_[node][i]);
    }

    //all conflicting nodes must be forbidden
    const std::vector<VertexDescriptor> &conflict_nodes = G->getConflictingNodes(node);
    for(Size i = 0; i < conflict_nodes.size(); ++ i)
    {
      forbidNode(conflict_nodes[i]);
    }
    std::cout<<"exit force node"<<std::endl;
  }



  void DeNovoLagrangeProblemBoost::forbidNode(VertexDescriptor node)
  {
    seqan::assignValue(forbidden_nodes_, node, true);
  }



  void DeNovoLagrangeProblemBoost::forbidEdge(EdgeDescriptor edge)
  {
    seqan::assignProperty(forbidden_edges_, edge, true);
  }



  void DeNovoLagrangeProblemBoost::forbidEdge(VertexDescriptor source, VertexDescriptor target)
  {
    if(EdgeDescriptor edge = seqan::findEdge(G->graph, source, target))
    {
      std::cerr<<"forbid edge: "<<source<<"   --->   " << target << "   " << edge << std::endl;
      forbidEdge(edge);
    }
  }



  void DeNovoLagrangeProblemBoost::getViolatedClusters(std::vector<Size> &clusters, const std::vector<VertexDescriptor> &path)
  {
    std::vector<int>violated_clusters(G->getNumberOfClusters(),-2);

    for(std::vector<VertexDescriptor>::const_iterator it = path.begin(); it != path.end(); ++it)
    {
      for(std::vector<Size>::const_iterator it_in = G->getClusters(*it).begin(); it_in != G->getClusters(*it).end(); ++it_in)
      {
        if(!(++violated_clusters[*it_in]) )
        {
          clusters.push_back(*it_in);
        }
      }
    }
  }



  void DeNovoLagrangeProblemBoost::updateWeights(const DVector& lambda)
  {
    edge_weights_ = edge_weights_bk_;

    SpectrumGraphSeqan::VertexIterator node_it(G->graph);

    for(; !seqan::atEnd(node_it); ++node_it)
    {
      DoubleReal lambda_sum = 0.0;
      const std::vector<Size> &clusters = G->getClusters(*node_it);
      for(std::vector<Size>::const_iterator clust_it = clusters.begin(); clust_it!=clusters.end(); ++clust_it)
      {
        lambda_sum += lambda[*clust_it];
      }

      if(lambda_sum > 0)
      {
        // adapt weight of all outgoing edges
        SpectrumGraphSeqan::OutEdgeIterator out_it(G->graph, *node_it);
        for(; !seqan::atEnd(out_it); ++out_it)
        {
          seqan::property(edge_weights_, *out_it) += lambda_sum;
        }
      }
    }
  }



  void DeNovoLagrangeProblemBoost::resetWeights(bool init)
  {
    if(init)
    {
      G->getEdgeWeights(edge_weights_bk_);
      edge_weights_ = edge_weights_bk_;

      SpectrumGraphSeqan::EdgeIterator edge_it(G->graph);
      spanning_edges_.assign(G->size(), std::vector<EdgeDescriptor>());

      while(!seqan::atEnd(edge_it))
      {
        VertexDescriptor source = seqan::sourceVertex(G->graph, *edge_it);
        VertexDescriptor target = seqan::targetVertex(G->graph, *edge_it);
        for(Size j = source + 1; j < target; ++j)
        {
          spanning_edges_[j].push_back(*edge_it);
        }
        ++edge_it;
      }
      return;
    }

    edge_weights_ = edge_weights_bk_;
  }


}//namespace
