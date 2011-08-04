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

#include <stdio.h>

#include <OpenMS/ANALYSIS/DENOVO/AntilopePMCorrect.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{

    //Default Constructor
    ParentMassCorrection::ParentMassCorrection():
            delta_(0.5),
            step_size_(0.1),
            range_(0){}

    ///Copy Constructor
    ParentMassCorrection::ParentMassCorrection(const ParentMassCorrection &pmc):
      delta_(pmc.delta_),
      step_size_(pmc.step_size_),
      range_(pmc.range_){}

    ///assignment operator
    ParentMassCorrection& ParentMassCorrection::operator=(const ParentMassCorrection &pmc)
    {
      if(this==&pmc) return *this;
			else
			{
        delta_=pmc.delta_;
        step_size_=pmc.step_size_;
        range_ = pmc.range_;
			}
			return *this;
		}



    //operator () calling the Correction function, returning the corrected parent mass. Also the parent mass of the inpute is changed
    DoubleReal ParentMassCorrection::operator () (PeakSpectrum &input_spec, DoubleReal intervall, DoubleReal step_size, DoubleReal delta)
    {
      range_ = intervall;
      step_size_ = step_size;
      delta_ = delta;

      return computeCorrection_(input_spec);
    }


    ///the actual heart of this class performing the correction
    DoubleReal ParentMassCorrection::computeCorrection_(PeakSpectrum &orig_spec)
    {
      //get the original parent mass
      DoubleReal orig_parent_mass=orig_spec.getPrecursors()[0].getPosition()[0];
      Size charge =orig_spec.getPrecursors()[0].getCharge();
      orig_parent_mass=(orig_parent_mass*charge)-((charge-1)*Constants::PROTON_MASS_U);

      Int number_of_steps=Int(range_/step_size_);

      std::vector<ShiftedSimilarity_> similarities;

      //for(Int i=-1*number_of_steps; i<number_of_steps;++i)
      for(Int i = -1*number_of_steps; i <= 0; ++i)
      {
        UInt loc_match_count=0;
        DoubleReal loc_mse=0;
        DoubleReal last_hit_mz=0.0;
        DoubleReal loc_match_intens=0.0;

        PeakSpectrum::ConstIterator forw_it=orig_spec.begin();

        PeakSpectrum::ConstReverseIterator rev_it=orig_spec.rbegin();

        while(forw_it!=orig_spec.end() && rev_it !=orig_spec.rend())
        {
          //std::cout<<"in while "<<i<<std::endl;
          DoubleReal mirror_mz=(1+orig_parent_mass + i*step_size_) - rev_it->getMZ();

          if(forw_it->getMZ() - mirror_mz < -delta_)
          {
            ++forw_it;
          }

          else if (forw_it->getMZ() - mirror_mz > delta_)
          {
            ++rev_it;
          }

          else
          {
            if(forw_it->getMZ() - last_hit_mz  >56.0)
            {
              ++loc_match_count;
              loc_match_intens+=forw_it->getIntensity() * rev_it->getIntensity();
              //std::cout<<"matching peaks: "<<forw_it->getMZ()<<"  "<<rev_it->getMZ()<<"   "<<orig_parent_mass + i*step_size_<<std::endl;
              last_hit_mz=forw_it->getMZ();
              loc_mse+=(forw_it->getMZ() - mirror_mz)*(forw_it->getMZ() - mirror_mz);
            }
            ++forw_it;
            ++rev_it;
          }
        }
        //std::cout<<"orig parent mass: "<<orig_parent_mass<<"  i "<<i<<"  step size  "<<step_size_<<std::endl; 
        similarities.push_back(ShiftedSimilarity_((orig_parent_mass + i*step_size_), loc_match_count, loc_match_intens, (loc_mse/loc_match_count)));
      }

      //debuggin
      
      for(UInt i = 0; i < similarities.size(); ++i)
      {
        std::cout<<"DEBUG: "<<similarities[i].shifted_parent_mass<<"  "<<similarities[i].shared_peaks_count<<"  "<<similarities[i].shared_peaks_intens<<"   "<<similarities[i].mse<<std::endl;
      }
      

      sort(similarities.begin(),similarities.end());

      corrected_mass_ = (similarities.front().shifted_parent_mass + charge - 1)/charge;

      return (corrected_mass_);

    }//computeCorrection_




}//namespace
