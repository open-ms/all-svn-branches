// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// '$'Maintainer: Sandro Andreotti '$'
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DENOVO/AntilopeBranchBound.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopeSubgradientSolver.h>

#define NDEBUG

namespace OpenMS
{
  //Custom Constructor
  BranchBoundDeNovo::BranchBoundDeNovo(DoubleReal *glb_ptr, const DeNovoLagrangeProblemBoost& dlp, path_score_pair* psp):
    glb_(glb_ptr),
    lp_(dlp),
    path_and_score_(psp)
    {}


   //Copy Constructor
  BranchBoundDeNovo::BranchBoundDeNovo(const BranchBoundDeNovo &bbd):
    glb_(bbd.glb_),
    lp_(bbd.lp_),
    path_and_score_(bbd.path_and_score_)
    {}

  //assignment operator
  BranchBoundDeNovo& BranchBoundDeNovo::operator=(const BranchBoundDeNovo &bbd)
  {
     if(this != &bbd)
     {
       glb_ = bbd.glb_;
       lp_ = bbd.lp_;
       path_and_score_ = bbd.path_and_score_;
     }
     return *this;
   }


  void BranchBoundDeNovo::solve()
  {
    SubgradientSolver sub_solver;
    lp_.set_lower_bound(max(lp_.get_lower_bound(),*glb_));
    sub_solver.InitProblem(1, lp_.num_of_duals(), DVector(lp_.num_of_duals(), 0), DVector(lp_.num_of_duals(), 100));

    SubgradientSolver::SolverProgress sp = sub_solver.Solve(lp_);

    if(sp != SubgradientSolver::SOLVED)
    {
      if(*glb_ < lp_.get_best_primal_score())
      {
        *glb_ = lp_.get_best_primal_score();
        *path_and_score_ = lp_.get_longest_path();
      }

      DeNovoLagrangeProblemBoost left_lp(lp_), right_lp(lp_);

      branch(left_lp, right_lp);
      std::cout<<"out of branch"<<std::endl;

      BranchBoundDeNovo left_node(glb_, left_lp, path_and_score_);
      left_node.solve();
      std::cout<<"left node done"<<std::endl;
      BranchBoundDeNovo right_node(glb_, right_lp, path_and_score_);
      right_node.solve();
      std::cout<<"right node done"<<std::endl;
    }

    else
    {
      if(*glb_ < lp_.get_best_primal_score())
      {
        std::cout<<"in solved if block"<<std::endl;
        *glb_ = lp_.get_best_primal_score();
        *path_and_score_ = lp_.get_longest_path();
      }
    }
  }

  void BranchBoundDeNovo::branch(DeNovoLagrangeProblemBoost &left_lp, DeNovoLagrangeProblemBoost &right_lp)
  {
    std::vector<VertexDescriptor> best_infeasible_path;
    lp_.getBestIllegalPath(best_infeasible_path);

    std::vector<Size> viol_clusters;
    lp_.getViolatedClusters(viol_clusters, best_infeasible_path);

    std::cerr<<"VIOL: "<<viol_clusters.size()<<std::endl;

    std::vector<VertexDescriptor> nodes_in_cl;
    lp_.getSpectrumGraph().getConflictingNodesByCluster(nodes_in_cl, viol_clusters.front());
    VertexDescriptor branch_node = nodes_in_cl.front();

    std::cout<<"perform branch step over node: "<<branch_node<<std::endl;

    //in the left lp we force the branch node to be deactivated
    left_lp.forbidNode(branch_node);

    //in the right lp we force the branch node to be active
    right_lp.forceNode(branch_node);
    std::cout<<"out of force node"<<std::endl;
  }





}//namespace
