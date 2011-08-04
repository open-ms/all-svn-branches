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

#include <OpenMS/ANALYSIS/DENOVO/AntilopePreprocessing.h>
#include <OpenMS/ANALYSIS/DENOVO/AntilopePMCorrect.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <cstring>
#include <set>
#include <algorithm>

namespace OpenMS {

  const  DoubleReal IdSetup::proton_m=Constants::PROTON_MASS_U;

  const  DoubleReal IdSetup::h2o_mass=EmpiricalFormula("H2O").getMonoWeight();
  const  DoubleReal IdSetup::nh3_mass=EmpiricalFormula("NH3").getMonoWeight();
  const  DoubleReal IdSetup::nh3_h2o_mass=h2o_mass+nh3_mass;

  const  DoubleReal IdSetup::BIon_offset=Residue::getInternalToFullMonoWeight()-Residue::getBIonToFullMonoWeight();
  const  DoubleReal IdSetup::BIon_h2o_offset=BIon_offset-h2o_mass;
  const  DoubleReal IdSetup::BIon_h2o_h2o_offset=BIon_offset-2*h2o_mass;
  const  DoubleReal IdSetup::BIon_nh3_offset=BIon_offset-nh3_mass;
  const  DoubleReal IdSetup::BIon_nh3_h2o_offset=BIon_offset - h2o_mass - nh3_mass;

  const  DoubleReal IdSetup::YIon_offset=Residue::getInternalToFullMonoWeight()-Residue::getYIonToFullMonoWeight();
  const  DoubleReal IdSetup::YIon_h2o_offset=YIon_offset - h2o_mass;
  const  DoubleReal IdSetup::YIon_h2o_h2o_offset=YIon_offset-2*h2o_mass;
  const  DoubleReal IdSetup::YIon_nh3_offset=YIon_offset-nh3_mass;
  const  DoubleReal IdSetup::YIon_nh3_h2o_offset=YIon_offset-h2o_mass-nh3_mass;

  const  DoubleReal IdSetup::AIon_offset=Residue::getInternalToFullMonoWeight()-Residue::getAIonToFullMonoWeight();
  const  DoubleReal IdSetup::AIon_h2o_offset=AIon_offset-h2o_mass;
  const  DoubleReal IdSetup::AIon_nh3_offset=AIon_offset-nh3_mass;

  const  DoubleReal IdSetup::CIon_offset=Residue::getInternalToFullMonoWeight()-Residue::getCIonToFullMonoWeight();

  const  DoubleReal IdSetup::ZIon_offset=Residue::getInternalToFullMonoWeight()-Residue::getZIonToFullMonoWeight();

  void IdSetup::nChooseKCombinations(int n, int k, std::vector<std::vector<int> > &combinations, std::vector<int> &tmp_comb)
  {
    if(n>k)
    {
        nChooseKCombinations(n-1,k, combinations, tmp_comb);
    }
    if(n>=k)
    {
      k--;
      tmp_comb[k]=n-1;
      if(k>0)
      {
        nChooseKCombinations(n,k, combinations, tmp_comb);
      }
      else
      {
        combinations.push_back(tmp_comb);
      }
    }
  }

  void IdSetup::nChooseKCombinations(int n, int k, std::vector<std::vector<int> > &combinations)
  {
    std::vector<int>tmp(k);
    nChooseKCombinations(n, k, combinations, tmp);
  }


int IdSetup::create_vector_A(BoolVec &A, UIntVec & A_len, AASeqVecMap & annot_map)
{
  UInt h = 400;
  UInt numelem = (UInt) (h / precision_) + 1;

  A.assign(numelem, false);
  A_len.assign(numelem, 0);

  ResidueDB *rdb=ResidueDB::getInstance();
  std::set<const Residue  *>all_residues_tmp=rdb->getResidues("Natural20");
  std::vector<String>all_amino_acids;

  //store one letter codes
  std::set<const Residue  *>::iterator it;
  for(it=all_residues_tmp.begin();it!=all_residues_tmp.end(); ++it)
  {
    String amino_acid=(*it)->getOneLetterCode();
    if(!((join_isobarics_ && amino_acid.hasPrefix("I")) || (join_Q_K_ &&  amino_acid.hasPrefix("Q"))))
    {
      all_amino_acids.push_back((*it)->getOneLetterCode());
    }
  }

  UInt number_aa=all_amino_acids.size();
  AASequence tag;
  std::vector<std::vector<int> >aa_tag_ids;

  for(Size len=1; len<=max_span_; ++len)
  {
    aa_tag_ids.clear();
    nChooseKCombinations(number_aa, len, aa_tag_ids);
    for(Size comb=0; comb<aa_tag_ids.size(); ++comb)
    {
      //get the AAs
      tag="";
      for(Size pos=0; pos<len;++pos)
      {
        tag+=all_amino_acids[aa_tag_ids[comb][pos]];
      }      
      int mass_index=(int)(tag.getMonoWeight(Residue::Internal)/precision_ + 0.5);
      for(Int offset=Int(-delta_/precision_); offset<=Int(delta_/precision_);++offset)
      {
        Size mass_index_off=mass_index+offset;
        if (mass_index_off<numelem)
        {
          A[mass_index_off]=true;
          annot_map[mass_index_off].push_back(tag);
          if(A_len[mass_index_off]!=1)
          {
            A_len[mass_index_off]=len;
          }
        }
      }
    }
  }
  return 0;
}

}//namespace


