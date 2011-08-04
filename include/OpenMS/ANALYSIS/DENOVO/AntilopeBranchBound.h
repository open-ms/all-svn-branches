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


#ifndef BRANCHBOUNDDENOVO_H_
#define BRANCHBOUNDDENOVO_H_

#include <OpenMS/ANALYSIS/DENOVO/AntilopeLagrangeProblem.h>

namespace OpenMS
{
  class BranchBoundDeNovo
  {
    typedef SpectrumGraphSeqan::VertexDescriptor VertexDescriptor;
    typedef SpectrumGraphSeqan::EdgeDescriptor EdgeDescriptor;
    typedef DeNovoLagrangeProblemBoost::path_score_pair path_score_pair;
    typedef std::vector<DoubleReal> DVector;

    public:
      ///Custom Constructor
      BranchBoundDeNovo(DoubleReal *glb_ptr, const DeNovoLagrangeProblemBoost& dlp, path_score_pair* psp);

      ///Copy Constructor
      BranchBoundDeNovo(const BranchBoundDeNovo &bbd);

      ///assignment operator
      BranchBoundDeNovo& operator=(const BranchBoundDeNovo &bbd);

      void solve();
      void branch(DeNovoLagrangeProblemBoost &, DeNovoLagrangeProblemBoost &);

    private:

      ///Default constructor
      BranchBoundDeNovo();

      ///Pointer to the global lower bound
      DoubleReal *glb_;

      ///the local DeNovoLagrangeProblem
      DeNovoLagrangeProblemBoost lp_;

      ///pointer to the best dual solution thus far
      path_score_pair *path_and_score_;
  };

}//namespace

#endif /* BRANCHBOUNDDENOVO_H_ */
