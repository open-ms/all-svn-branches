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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>

#include <cmath>
#include <cfloat>
using namespace std;

namespace OpenMS
{
  SpectraSTSimilarityScore::SpectraSTSimilarityScore()
    : PeakSpectrumCompareFunctor()
  {
		setName(SpectraSTSimilarityScore::getProductName());
  }

  SpectraSTSimilarityScore::SpectraSTSimilarityScore(const SpectraSTSimilarityScore& source)
    : PeakSpectrumCompareFunctor(source)
  {
  }

  SpectraSTSimilarityScore::~SpectraSTSimilarityScore()
  {
  }

  SpectraSTSimilarityScore& SpectraSTSimilarityScore::operator = (const SpectraSTSimilarityScore& source)
  {
		if (this != &source)
		{
    	PeakSpectrumCompareFunctor::operator = (source);
		}
    return *this;
  }

	double SpectraSTSimilarityScore::operator () (const PeakSpectrum& spec) const
	{
		return operator () (spec, spec);
	}
	
  double SpectraSTSimilarityScore::operator () (const PeakSpectrum& s1, const PeakSpectrum& s2) const
  {
		double score(0);
		BinnedSpectrum bin1(1,1,s1);
		BinnedSpectrum bin2(1,1,s2);
		
		//normalize bins
		
		//magnitute of the spectral vector
		Real magnitude1(0);
		Real magnitude2(0);
		for(SparseVector<Real>::SparseVectorIterator iter1 = bin1.getBins().begin(); iter1 < bin1.getBins().end(); ++iter1)
		{
			magnitude1 += pow((double)*iter1,2);
		}
		magnitude1 = sqrt(magnitude1);
		//normalize bins of bin1
		for(SparseVector<Real>::SparseVectorIterator iter1 = bin1.getBins().begin(); iter1 < bin1.getBins().end(); ++iter1)
		{
			*iter1 = (Real)*iter1/magnitude1;
		}
		
		for(SparseVector<Real>::SparseVectorIterator iter2 = bin2.getBins().begin(); iter2 < bin1.getBins().end(); ++iter2)
		{
			magnitude2 += pow((double)*iter2,2);
		}
		magnitude2 = sqrt(magnitude2);
		//normalize bins of bin1
		for(SparseVector<Real>::SparseVectorIterator iter2 = bin2.getBins().begin(); iter2 < bin2.getBins().end(); ++iter2)
		{
			*iter2 = (Real)*iter2/magnitude2;
		}		
		
		Size shared_bins = min(bin1.getBinNumber(),bin2.getBinNumber());
		for(Size s = 0; s < shared_bins; ++s)
		{
			if((double)bin1.getBins()[s] >0.0 && (double)bin2.getBins()[s]>0.0)
			{
				score += ((double)bin1.getBins()[s]*(double)bin2.getBins()[s]);
			}
		}	
		
    return score;
	
	}
	
	double SpectraSTSimilarityScore::operator() (const BinnedSpectrum& bin1,const BinnedSpectrum& bin2)	const
	{
		double score(0);
		
		Size shared_bins = min(bin1.getBinNumber(),bin2.getBinNumber());
		for(Size s = 0; s < shared_bins; ++s)
		{
			if(bin1.getBins()[s] >0 && bin2.getBins()[s]>0)
			{
				score += (bin1.getBins()[s]*bin2.getBins()[s]);
			}
		}	
		
    return score;	
	}
	
	bool SpectraSTSimilarityScore::preprocess(PeakSpectrum& spec, Real remove_peak_intensity_threshold, Real filter_all_peaks_below_mz,Real min_peak_number)
	{
		UInt total_intensity(0), intensity_below_mz(0);
		for(PeakSpectrum::iterator k = spec.begin(); k < spec.end(); ++k)
		{
			if(k->getIntensity() >  remove_peak_intensity_threshold)
			{
				total_intensity +=  k->getIntensity();
				k->setIntensity(sqrt(k->getIntensity()));
				if (k->getMZ() < filter_all_peaks_below_mz) 
				{
     			intensity_below_mz += k->getIntensity();
    		}
			}
			else
			{
				k = spec.erase(k);
			}
		}
		//if not enough peaks in the specturm pass that one out
		//Remove spectra with almost no peaks above a certain m/z value. 
		//All query spectra with 95%+ of the total intensity below <m/z> will be removed.
		if(spec.size() >= min_peak_number && intensity_below_mz < 0.95 * total_intensity)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	
	BinnedSpectrum SpectraSTSimilarityScore::transform(const PeakSpectrum& spec)
	{
		BinnedSpectrum bin(1,1,spec);
		Real magnitude(0);
		for(SparseVector<Real>::SparseVectorIterator iter = bin.getBins().begin(); iter < bin.getBins().end(); ++iter)
		{
			magnitude += pow((double)*iter,2);
		}
		magnitude = sqrt(magnitude);
		//normalize bins
		for(SparseVector<Real>::SparseVectorIterator iter = bin.getBins().begin(); iter < bin.getBins().end(); ++iter)
		{
			*iter = (Real)*iter/magnitude;
		}	
		return bin;
	}
	
	double SpectraSTSimilarityScore::dot_bias(const BinnedSpectrum& bin1, const BinnedSpectrum& bin2, double dot_product) const
	{
		double numerator(0);

		Size shared_bins = min(bin1.getBinNumber(),bin2.getBinNumber());
		for(Size s = 0; s < shared_bins; ++s)
		{
			if(bin1.getBins()[s] >0 && bin2.getBins()[s]>0)
			{
				numerator += (pow(bin1.getBins()[s],2)*pow(bin2.getBins()[s],2));
			}
		}			
		numerator = sqrt(numerator);
		
		if(dot_product)
		{
			return (double)numerator/dot_product;
		}
		else
		{
			return (double)numerator/(*this)(bin1,bin2);
		}
	}
	
	double SpectraSTSimilarityScore::delta_D(double top_hit, double runner_up)
	{
		return (double)(top_hit-runner_up)/top_hit;
	}
	
	double SpectraSTSimilarityScore::compute_F(double dot_product, double delta_D,double dot_bias)
	{
		double b(0);
		if(dot_bias < 0.1 || ( 0.35 < dot_bias && dot_bias <= 0.4))
		{
			b = 0.12;
		}
		else if( 0.4 < dot_bias && dot_bias <= 0.45)
		{
			b = 0.18;
		}
		else if(dot_bias > 0.45)
		{
			b  = 0.24;
		}
		return 0.6*dot_product + 0.4*delta_D - b;
	}

	
}
