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
// $Maintainer: Sandro Andreotti $
// --------------------------------------------------------------------------

#ifndef ANTILOPEILP_H
#define ANTILOPEILP_H


#include <OpenMS/ANALYSIS/DENOVO/AntilopeSpectrumGraph.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>


namespace OpenMS
{

  class de_novo_ILP: public DefaultParamHandler
  {
  private:
    //convenience typedef
    typedef std::vector<DoubleReal>DRealVec;
    typedef std::vector<UInt>UIntVec;
    typedef std::vector<UIntVec>CandSetInt;
    typedef std::vector<DRealVec>CandSetReal;

  public:
    ///Constructor
    de_novo_ILP();

    ///Copy Constructor
    de_novo_ILP(const de_novo_ILP&);

    ///destructor
    virtual ~de_novo_ILP()
    {
    };

    ///assignment Operator
    de_novo_ILP & operator=(const de_novo_ILP&);

    /// build cplex model and solve - classic approach
    void computeCandidates(const SpectrumGraphSeqan &graph, CandSetInt &cands, DRealVec &scores);

  }; //de_novo_ILP
}//namespace

#endif // ANTILOPEILP_H
