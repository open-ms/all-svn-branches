// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EuclideanDistance.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

	EuclideanDistance::EuclideanDistance()
		:	BaseQuality()
	{
		setName(getProductName());
		check_defaults_ = false;
	}

	EuclideanDistance::~EuclideanDistance()
	{
	}

	double EuclideanDistance::evaluate(const IndexSet& set, const BaseModel<2>& model)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		
		typedef BaseModel<2>::IntensityType Intensity;
		
		Intensity temp_diff = 0;		// differences between different coordinate vectors
		Intensity sum_diff   = 0;		// sum of differences		
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			temp_diff = (model.getIntensity(traits_->getPeakPos(*it)) - traits_->getPeakIntensity(*it));
			sum_diff += temp_diff * temp_diff;
		}
		
		return (0 - sqrt(sum_diff));
	}
	
	double EuclideanDistance::evaluate(const IndexSet& set, const BaseModel<1>& model, UInt dim)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		
		typedef BaseModel<2>::IntensityType Intensity;
		
		Intensity temp_diff = 0;		// differences between different coordinate vectors
		Intensity sum_diff   = 0;		// sum of differences		
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			temp_diff = (model.getIntensity(traits_->getPeakPos(*it)[dim]) - traits_->getPeakIntensity(*it));			
			sum_diff += temp_diff * temp_diff;
		}
		
		return (0-sqrt(sum_diff));
	}

	double EuclideanDistance::evaluate(const Math::LinearInterpolation<double,double>& lint, const BaseModel<1>& iso_model)
	{
			// compute average intensity and intensity sum of spectrum
			double mz_data_sum = 0.0;
			for (UInt i=0;i<lint.getData().size();++i)
			{
				mz_data_sum += lint.getData()[i];
			}		
		
			// compute model intensity sum 
			double mz_model_sum = 0.0;		
			for (UInt i=0;i<lint.getData().size();++i)
			{
				mz_model_sum += iso_model.getIntensity(lint.index2key(i) );
			} 		
		
			double temp_diff = 0;		// differences between different coordinate vectors
			double sum_diff   = 0;		// sum of differences		
			for (UInt i=0;i<lint.getData().size();++i)
			{
					temp_diff = (iso_model.getIntensity(lint.index2key(i) ) / mz_model_sum) - (lint.getData()[i] / mz_data_sum);
					sum_diff +=  temp_diff * temp_diff;
			}
					
			return (0-sqrt(sum_diff));
	}
	
}
