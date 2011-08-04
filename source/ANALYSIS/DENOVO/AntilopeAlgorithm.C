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

#define DEBUG
namespace OpenMS
{  
  typedef DeNovoLagrangeProblemBoost::path_score_pair path_score_pair;

  path_score_pair YenAlgorithm::kthIteration(UInt k, SubgradientSolver * sub_solver)
  {
    path_score_pair result_path_score, prefix;
    result_path_score.score = 0;
    prefix.score = 0;
    DeNovoLagrangeProblemBoost de_novo_lagrange(G);
    forbid_edges(k, de_novo_lagrange);    

    if (k == 0)
    {
      DoubleReal global_lower_bound = DeNovoLagrangeProblemBoost::MINUS_INF;

      BranchBoundDeNovo branch_bound(&global_lower_bound, de_novo_lagrange, &result_path_score);
      branch_bound.solve();

      cand_path_heap.insert(result_path_score);
    }
    else
    {
      prefix.path.push_back(0);
      for(Size i = 1; i <= deviation_node;  ++i)
      {
        de_novo_lagrange.forbidConflictingNodes(found_path[k - 1][i]);
        prefix.path.push_back(found_path[k - 1][i]);
        prefix.score += G->getEdgeWeight(found_path[k-1][i-1], found_path[k-1][i]);
      }

      cout << "DEVIATION node: " << deviation_node << endl;
      for (int i = deviation_node; i < found_path[k - 1].size() - 1; ++i)
      {
        result_path_score.score = DeNovoLagrangeProblemBoost::MINUS_INF;
        result_path_score.path.clear();
        cout << "iteration: " << i << endl;

        de_novo_lagrange.setPrefix(prefix);
        de_novo_lagrange.forbidConflictingNodes(prefix.path.back());
        de_novo_lagrange.forbidEdge(found_path[k-1][i], found_path[k-1][i+1]);

        de_novo_lagrange.reset();
        de_novo_lagrange.resetWeights(false);
        de_novo_lagrange.set_lower_bound(cand_path_heap.rbegin()->score);

        DoubleReal global_lower_bound = DeNovoLagrangeProblemBoost::MINUS_INF;
        BranchBoundDeNovo branch_bound(&global_lower_bound, de_novo_lagrange, &result_path_score);
        branch_bound.solve();

        float score = result_path_score.score;
        float qu_score = cand_path_heap.rbegin()->score;

        //if path has better score than lowest score in heap, the replace this lowest heap element with it
        if (score > qu_score)
        {
          std::cout<<"score vgl: "<<score<<" :: "<<qu_score<<std::endl;
          cand_path_heap.insert(result_path_score);
          cand_path_heap.erase(--cand_path_heap.end());

          std::cout<<"added candidate path:"<<std::endl;
        }
        prefix.path.push_back(found_path[k - 1][i+1]);
        prefix.score += G->getEdgeWeight(found_path[k-1][i], found_path[k-1][i+1]);
      }
    }    

    std::cerr<<"TOP CAND SCORE: "<< cand_path_heap.begin()->score <<std::endl;
    std::cerr<<"CAND SIZE: "<<cand_path_heap.size()<<std::endl;

    path_score_pair k_longest = *cand_path_heap.begin();
    cand_path_heap.erase(cand_path_heap.begin());

    return k_longest;
  }



  //forbid edges as defined in yen paper page 714, step a
  int YenAlgorithm::forbid_edges(int k, DeNovoLagrangeProblemBoost &dnlp)
  {    
    deviation_node = 0;

    if (k == 0)
    {
      return 0;
    }

    else
    {
      vector<VertexDescriptor> last_path = found_path[k - 1];

      //forbid edges from paths in found_path
      for (int i = 0; i < k-1; ++i)
      {
        int count = 0;
        vector<VertexDescriptor>::iterator last_path_it = last_path.begin();
        vector<VertexDescriptor>::iterator found_path_it = found_path[i].begin();

        while(*last_path_it == *found_path_it)
        {
          ++last_path_it;
          ++found_path_it;
          ++count;
        }
        dnlp.forbidEdge(*(--last_path_it), *found_path_it);
        cout << "forbid edge: " << *last_path_it << "::" << *found_path_it << endl;
        deviation_node = max(deviation_node, count - 1);
      }      
    }
    return 0;
  }



  int YenAlgorithm::computeLongestPaths(int k)
  {
    cand_path_heap.clear();
    path_score_pair tmp = { vector<VertexDescriptor>(0), -SpectrumGraphSeqan::INFINITYdist };
    std::vector<path_score_pair>tmp_v(k, tmp);
    cand_path_heap.insert(tmp_v.begin(), tmp_v.end());

    SubgradientSolver *sub_solver = new SubgradientSolver();

    double tstart = clock();
    double time_total = 0;

    for (int i = 0; i < k; i++)
    {
      cout << "start for " << i << "-th path" << endl;
      path_score_pair longestPath_temp = kthIteration(i, sub_solver);

      found_path_scores.push_back(longestPath_temp.score);
      found_path.push_back(longestPath_temp.path);
      cout << "found path number: " << i <<"  with score  "<<found_path_scores[i]<<endl;

      for(UInt tt=0; tt<found_path[i].size();++tt)
      {
        std::cout<<" - "<<found_path[i][tt]<<std::endl;
      }

      double time = (clock() - tstart) / CLOCKS_PER_SEC;
      time_total = time;
      cout << "TIME_LAGRANGE::" << i << "::" << time << "::" << time_total << endl;

    }
    cout << " TOTAL TIME_LAGRANGE:: " <<time_total << endl;

    for (int i = 0; i < k; i++)
    {
      vector<VertexDescriptor>::iterator found_path_it = found_path[i].begin();
      for (;found_path_it !=found_path[i].end(); ++found_path_it)
      {
        cout << *found_path_it << endl;
      }
      cout << found_path_scores[i] << endl;
    }
    return 0;
  }


}//namespace
