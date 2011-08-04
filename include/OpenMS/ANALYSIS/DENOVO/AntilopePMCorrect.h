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

#ifndef _PARENTMASSCORRECTION_H
#define	_PARENTMASSCORRECTION_H


#include <OpenMS/KERNEL/StandardTypes.h>


namespace OpenMS
{


  class ParentMassCorrection
  {

  public:

    ///Default Constructor
    ParentMassCorrection();

    ///Copy Constructor
    ParentMassCorrection(const ParentMassCorrection &p_m_c);

    ///assignment operator
    ParentMassCorrection& operator=(const ParentMassCorrection &pmc_i);

    ///() Operator to call the computation corrected parent manss
    DoubleReal operator () (PeakSpectrum &input_spec, DoubleReal intervall=2, DoubleReal step_size=0.1, DoubleReal delta=0.5);

    ///virtual destructor
    virtual ~ParentMassCorrection(){};

    //getter & setter

    ///return the step_size
    DoubleReal getStepSize(){return step_size_;};

    ///set the step_size
    void setStepSize(DoubleReal step_in){step_size_=step_in;};

    ///return delta
    DoubleReal getDelta(){return delta_;};

    ///set delta
    void setDelta(DoubleReal delta_in){delta_=delta_in;};

    ///get the corrected mass
    DoubleReal getCorrectedMass(){return corrected_mass_;};




  protected:

    //similarty typedef
    struct ShiftedSimilarity_
    {
      DoubleReal shifted_parent_mass;
      UInt shared_peaks_count;
      DoubleReal shared_peaks_intens;
      DoubleReal mse;

      //constructor
      ShiftedSimilarity_(DoubleReal p_mass, UInt spc,DoubleReal intens_in, DoubleReal mse_in):
              shifted_parent_mass(p_mass),
              shared_peaks_count(spc),
              shared_peaks_intens(intens_in),
              mse(mse_in){};

      bool operator<(const ShiftedSimilarity_ &s2)const
      {
        if(shared_peaks_intens != s2.shared_peaks_intens)
        {
          return shared_peaks_intens > s2.shared_peaks_intens;
        }
        else
        {
          return mse < s2.mse;
        }
      }
    };

    ///the actual hear of this class performing the correction
    DoubleReal computeCorrection_(PeakSpectrum  &spectrum);

    DoubleReal delta_;

    DoubleReal step_size_;

    DoubleReal range_;

    DoubleReal corrected_mass_;








  };


}

#endif	/* _PARENTMASSCORRECTION_H */

