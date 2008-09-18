// -*- C++: make; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_RANKCORRELATION_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_RANKCORRELATION_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseQuality.h>

#include "gsl/gsl_cdf.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

#include <math.h>
#include<limits>

namespace OpenMS
{
	/** 
		@brief Measures the quality of a modelfit to some realworld data. 
	  
 		Implementation of class BaseQuality. The quality is measured as
 		the (squared) Spearman correlation coefficent between data and model.
		
		The Spearman correlation coefficient is a non-parametric measure of
		correlation.
 		
 		@ingroup FeatureFinder
 	*/
  class RankCorrelation 
    : public BaseQuality
  {

  typedef BaseModel<2>::IntensityType IntensityType;
  typedef BaseModel<2>::CoordinateType CoordinateType;
  typedef std::vector<IntensityType> IntensityVector;
  
  public:
    /// Default constructor
    RankCorrelation();

    /// destructor 
    virtual ~RankCorrelation();

    /// evaluates the quality of the fit of @p model to @p set
    double evaluate(const IndexSet& set, const BaseModel<2>& model);
    
    /// evaluates the quality of the fit of @p model to @p set along dimension @p dim
   	double evaluate(const IndexSet& set, const BaseModel<1>& model, UInt dim);
		
		/// evaluates the quality of the fit of @p model to @p lint along dimension @p dim
		double evaluate(const Math::LinearInterpolation<double,double>& lint, const BaseModel<1>& model);

		/// creates instance of this class (function is called by factory).
   	static BaseQuality* create()
   	{ 
   		return new RankCorrelation();
   	}

		/// name 
	  static const String getProductName()
	  { 
	  	return "RankCorrelation"; 
	  }
		
  }; //class
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_RANKCORRELATION_H
