// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
#ifndef OPENMS_ANALYSIS_TARGETED_RTPROBABILITY_H
#define OPENMS_ANALYSIS_TARGETED_RTPROBABILITY_H


#include <OpenMS/KERNEL/FeatureMap.h>



namespace OpenMS 
{
  /**
    @brief RTProbability calculates probabilities 

    @htmlinclude OpenMS_RTProbability.parameters

    @ingroup Analysis_TARGETED
  */
  class OPENMS_DLLAPI RTProbability : public DefaultParamHandler
  {
  public:
    
    /// Default constructor
    RTProbability();
    
    /// Copy constructor
    RTProbability(const RTProbability& rhs);
    
    /// Desctructor
    virtual ~RTProbability();
    
    /// assignment operator
    RTProbability& operator = (const RTProbability& rhs);
        
    void learnGaussian(FeatureMap<>& features,String rt_model_path,DoubleReal min_score=0.5);
    
    DoubleReal getRTProbability(DoubleReal min_obs_rt,DoubleReal max_obs_rt, DoubleReal pred_rt);
		void setGaussianParameters(DoubleReal mu, DoubleReal sigma);

    DoubleReal getGaussSigma()
    {
      return sigma_;
    }
    DoubleReal getGaussMu()
    {
      return mu_;
    }
    
  protected:

    Int getScanNumber_(DoubleReal rt);

    void generateDistributionImage_(const std::vector<DoubleReal>& ids,
                                    const String& formula, const String& filename,DoubleReal min, DoubleReal diff);
    
    DoubleReal sigma_;
    DoubleReal mu_;

    
  };
 
} // namespace OpenMS

#endif//#ifndef OPENMS_ANALYSIS_TARGETED_RTPROBABILITY_H
