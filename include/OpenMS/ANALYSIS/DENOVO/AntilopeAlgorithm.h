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
    typedef DeNovoLagrangeProblemBoost::PathSolution PathSolution;
    typedef std::vector<PathSolution> PathSolutionList;
    
    typedef std::vector<VertexDescriptor> Path;
    typedef std::vector<Path> PathList;
    typedef std::vector<DoubleReal> PathScoreList;


    private:
    
      enum HeuristicReturn
      {
        FEASIBLE,
        INFEASIBLE,
        PRUNED,
        NOPATH
      };

      ///de_novoGraph
      SpectrumGraphSeqan* g_;
    
      ///the heap for the candidate paths
      std::multiset<PathSolution> cand_heap_;
    
      seqan::String<DoubleReal> dist_backward;
      seqan::String<VertexDescriptor> pred_backward;

    
    public:
      ///constructor
      YenAlgorithm(SpectrumGraphSeqan* g_in) :
        g_(g_in)
      {
      }
    
      ///compute the k longest paths
      int computeLongestPaths(Size k, PathList& paths, PathScoreList& scores);
    
    
    protected:

      SpectrumGraphSeqan::TGraph graph_rev;

      //Member Functions
      ///forbid edges
      Size forbidEdges_(Size k, const PathSolutionList& found_paths, DeNovoLagrangeProblemBoost& dnlp);

      ///the k-th iteration step
      void singleIteration_(Size k,
                                  const PathSolutionList& found_paths,
                                  PathSolution &path);
    
    HeuristicReturn backWardHeuristic1(const DeNovoLagrangeProblemBoost& de_novo_lagrange,
                                       const seqan::String<DoubleReal>& dist_backward,
                                       const seqan::String<VertexDescriptor>& pred_backward,
                                       DoubleReal lower_bound,
                                       PathSolution& cand_path
                                       );
    
   void backWardHeuristic2(const DeNovoLagrangeProblemBoost& de_novo_lagrange,
                           const seqan::String<DoubleReal>& weights,
                           SpectrumGraphSeqan::TGraph& rev_graph,
                           DoubleReal lower_bound,
                           Size s,
                           std::vector<std::pair<HeuristicReturn, PathSolution> >& candidates,
                           const std::vector<VertexDescriptor>& parent_path
                           );


  };

}


#endif

