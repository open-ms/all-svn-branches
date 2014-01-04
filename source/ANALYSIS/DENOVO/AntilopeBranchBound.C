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
  BranchBoundDeNovo::BranchBoundDeNovo(DoubleReal *glb_ptr, const DeNovoLagrangeProblemBoost& dlp, PathSolution* psp):
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


  void BranchBoundDeNovo::solve(const DVector &init_dual_vec)
  {
    SubgradientSolver sub_solver;
    lp_.setLowerBound(max(lp_.getLowerBound(),*glb_));
    sub_solver.initProblem(1, lp_.getDualDim(), DVector(lp_.getDualDim(), 0), DVector(lp_.getDualDim(), 100));

    if (!init_dual_vec.empty())
      sub_solver.setDualVector(init_dual_vec);

    SubgradientSolver::SolverProgress sp = sub_solver.solve(lp_);

    std::vector<VertexDescriptor> best_infeasible_path;
    
    if(sp != SubgradientSolver::SOLVED)
    {
      lp_.getBestIllegalPath(best_infeasible_path);
      if (best_infeasible_path.empty())
      {
        std::cout << "ATTENTION REGET INFEASIBLE! " << std::endl;
        //if no infeasible solution was found (is possible as we do not start with zero Lambda in branch nodes)
        //solve a single iteration with lambda zero vector.
        sub_solver.setDualVector(DVector(lp_.getDualDim(), 0));
        Param pars = sub_solver.getDefaults();
        pars.setValue("noofiterations", 1);
        sub_solver.setParameters(pars);
        sp = sub_solver.solve(lp_);
        
        lp_.getBestIllegalPath(best_infeasible_path); // now it must be non_empty or problem is feasible
        std::cout << best_infeasible_path.size() << std::endl;
//        exit(1);
      }
    }
    
    if(sp != SubgradientSolver::SOLVED)
    {
      if(*glb_ < lp_.get_best_primal_score())
      {
        *glb_ = lp_.get_best_primal_score();
        *path_and_score_ = lp_.get_longest_path();
      }

      //perform branching
      
      DeNovoLagrangeProblemBoost left_lp(lp_), right_lp(lp_);
      left_lp.resetBestFeasible();
      left_lp.resetBestInfeasible();
      
      right_lp.resetBestFeasible();
      right_lp.resetBestInfeasible();

      VertexDescriptor branch_node = lp_.getBranchNode(best_infeasible_path);
      
      std::cout<<"perform branch step over node: "<<branch_node<<std::endl;
      
      //in the left lp we force the branch node to be deactivated
      left_lp.forbidNode(branch_node);
      
      //in the right lp we force the branch node to be active
      right_lp.forceNode(branch_node);
      
//      branch(left_lp, right_lp);


      BranchBoundDeNovo left_node(glb_, left_lp, path_and_score_);
      left_node.solve(sub_solver.getDualVector());
      std::cout<<"left node done"<<std::endl;
      BranchBoundDeNovo right_node(glb_, right_lp, path_and_score_);
      right_node.solve(sub_solver.getDualVector());
      std::cout<<"right node done"<<std::endl;
    }

    else
    {
      if(*glb_ < lp_.get_best_primal_score())
      {
#ifdef DEBUG
        std::cout<<"in solved if block"<<std::endl;
#endif
        *glb_ = lp_.get_best_primal_score();
        *path_and_score_ = lp_.get_longest_path();
      }
    }
  }

//  void BranchBoundDeNovo::branch(DeNovoLagrangeProblemBoost &left_lp, DeNovoLagrangeProblemBoost &right_lp)
//  {
//    std::vector<VertexDescriptor> best_infeasible_path;
//    lp_.getBestIllegalPath(best_infeasible_path);
//
//    //std::vector<Size> viol_clusters;
//    //lp_.getViolatedClusters(viol_clusters, best_infeasible_path);
//
//    //std::cerr<<"VIOL: "<<viol_clusters.size()<<std::endl;
//
//    //const std::vector<VertexDescriptor> &nodes_in_cl = lp_.getSpectrumGraph().getCluster(viol_clusters.front());
//    
//    //VertexDescriptor branch_node = nodes_in_cl.front();
//    
//    VertexDescriptor branch_node = lp_.getBranchNode(best_infeasible_path);    
//
//    std::cout<<"perform branch step over node: "<<branch_node<<std::endl;
//
//    //in the left lp we force the branch node to be deactivated
//    left_lp.forbidNode(branch_node);
//
//    //in the right lp we force the branch node to be active
//    right_lp.forceNode(branch_node);
//    std::cout<<"out of force node"<<std::endl;
//  }





}//namespace
