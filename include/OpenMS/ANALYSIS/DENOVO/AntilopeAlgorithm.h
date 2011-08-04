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

#ifndef LISAK_LONGESTPATH_lagrange_H
#define LISAK_LONGESTPATH_lagrange_H

//#include <Common/CommonTypes.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeSpectrumGraph.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeLagrangeProblem.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeSubgradientSolver.h>

namespace OpenMS
{
  class YenAlgorithm
  {
    typedef SpectrumGraphSeqan::VertexDescriptor VertexDescriptor;
    typedef DeNovoLagrangeProblemBoost::path_score_pair path_score_pair;

    private:
      ///de_novoGraph
      SpectrumGraphSeqan* G;
      ///number of clusters
      int clust_num;
      ///the first node to be used in the k-th iteration
      int deviation_node;
      ///the forbidden edges (one vector for each starting nodes)
      vector<vector<bool> > forbidden_edges;
      ///the scores of the computed paths
      vector<double> found_path_scores;
      ///the computed paths as vector of nodes
      vector<vector<VertexDescriptor> > found_path;
      ///the heap for the candidate paths
      std::multiset<path_score_pair> cand_path_heap;
      ///the Lagrange Problem is used once for each starting node
      DeNovoLagrangeProblemBoost *de_novo_lagrange;

    public:
      ///constructor
      YenAlgorithm(SpectrumGraphSeqan* G_in) :
        G(G_in)
      {
        deviation_node = 0;
      }

      //Member Functions
      ///forbid edges
      int forbid_edges(int k, DeNovoLagrangeProblemBoost &dnlp);

      ///compute the k longest paths
      int computeLongestPaths(int k);

      ///the k-th iteration step
      path_score_pair kthIteration(UInt k, SubgradientSolver * sub_solver);

      std::vector<std::vector<VertexDescriptor> > get_longest_paths()
      {
        return found_path;
      }
      vector<double> get_path_scores()
      {
        return found_path_scores;
      }
  };

}
;//namespace

#endif

