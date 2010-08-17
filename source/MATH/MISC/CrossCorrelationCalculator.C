// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/CrossCorrelationCalculator.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
 #include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_randist.h>
#include <cmath>

namespace OpenMS
{
CrossCorrelationCalculator::CrossCorrelationCalculator(DoubleReal stepwidth,DoubleReal gauss_mean, DoubleReal gauss_sigma) : stepwidth_(stepwidth), gauss_mean_(gauss_mean), gauss_sigma_(gauss_sigma) {
	// TODO Auto-generated constructor stub

}

CrossCorrelationCalculator::~CrossCorrelationCalculator() {
	// TODO Auto-generated destructor stub
}

void CrossCorrelationCalculator::calculate(const MSExperiment<Peak1D>& input, MSExperiment<Peak1D>& output,Size spectrum_selection_id)
{
	mz_min=input.getMinMZ();
	mz_max=input.getMaxMZ();
	// make sure that output is clear
	output.clear(true);

	// copy experimental settings
	static_cast<ExperimentalSettings&>(output) = input;

	// resize output with respect to input
	output.resize(input.size()-mz_min);

	Size progress = 0;

	startProgress(0,input.size(),"Calculating Autocorrelation");
	for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
	{
		if (input[scan_idx].getMSLevel() != 1)
		{
			output[scan_idx] = input[scan_idx];
		}
		else if (scan_idx==spectrum_selection_id)
		{
			analyzeSpectrum(input[scan_idx], output[scan_idx],true);
		}
		else
		{
			analyzeSpectrum(input[scan_idx], output[scan_idx]);
		}
		setProgress(++progress);
	}
	endProgress();
}

void CrossCorrelationCalculator::analyzeSpectrum(const MSSpectrum<Peak1D>& input, MSSpectrum<Peak1D>& output, bool gauss_fitting)
{
	//Copy spectrum settings
	output.clear(true);
	output.SpectrumSettings::operator=(input);
	output.MetaInfoInterface::operator=(input);
	output.setRT(input.getRT());
	output.setMSLevel(input.getMSLevel());
	output.setName(input.getName());
	output.setType(SpectrumSettings::PEAKS);

	std::vector<DoubleReal> mz_vec;
	std::vector<DoubleReal> intensity_vec;
	mz_min=input.begin()->getMZ();

	//Fill intensity and m/z vector for interpolation
	for (MSSpectrum<Peak1D>::ConstIterator mz_it=input.begin(); mz_it!=input.end(); ++mz_it)
	{
		mz_vec.push_back(mz_it->getMZ());
		intensity_vec.push_back(mz_it->getIntensity());
	}

	//Interpolate
	gsl_interp_accel* acc_spl = gsl_interp_accel_alloc();
	gsl_spline* spline_spl = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
	gsl_spline_init(spline_spl, &(*mz_vec.begin()), &(*intensity_vec.begin()), mz_vec.size());

	DoubleReal max_distance=mz_max-mz_min;

	//n must be a power of two for gsl Fourier transformation; take next higher size for n, which is a power of two
	Size s = pow(2,(ceil(log2(max_distance/stepwidth_))));

	//Create a vector containing interpolated values with a spacing of "stepwidth"
	std::vector<DoubleReal> data(s,0);
	Size i=0;
	for (DoubleReal x=mz_min;x<=mz_max;x+=stepwidth_)
	{
		data[i] = gsl_spline_eval (spline_spl, x, acc_spl);
		++i;
	}


	//Fast fourier transform the values
	gsl_fft_real_radix2_transform (&(*data.begin()), 1, s);

	std::vector<DoubleReal> data1;


	if (gauss_fitting)
	{
		std::vector<DoubleReal> gauss_fitted_data(s,0);
		Size j=0;
		for (DoubleReal x=0.0-(gauss_mean_-mz_min);x<=0.0+(mz_max-gauss_mean_);x+=stepwidth_)
		{
			gauss_fitted_data[j] = gsl_spline_eval (spline_spl, gauss_mean_+x, acc_spl)*gsl_ran_gaussian_pdf (x, 0.5*gauss_sigma_);
//			std::cout << j << "\t" << gauss_mean_+x << "\t" << x << "\t" << gsl_spline_eval (spline_spl, gauss_mean_+x, acc_spl) << "\t" << gsl_ran_gaussian_pdf (x, 0.5*gauss_sigma_) << std::endl;
			++j;
		}
//		for (j=0;j<gauss_fitted_data.size();++j)
//		{
//			std::cout << j << "\t" << gauss_fitted_data[j] << std::endl;
//		}
		gsl_fft_real_radix2_transform (&(*gauss_fitted_data.begin()), 1, s);
		data1.insert(data1.end(),gauss_fitted_data.begin(),gauss_fitted_data.end());
	}
	else
	{
		data1.insert(data1.end(),data.begin(),data.end());
	}
//	if (gauss_fitting)
//	{
//		for (Size j=0;j<f.size();++j)
//		{
//			std::cout << f[j] << "\t" << data[j] << std::endl;
//		}
//	}

	//Multiply the fast fourier transformed complex values with the complex conjugate
	//Have a look at the GNU Scientific Library reference manual for a description of the data structure
	data[0]=data[0]*data1[0];
	data[s/2]=data[s/2]*data1[s/2];
	for (i = 1; i <= s/2; ++i)
	{
		data[i]=data[i]*data1[i]+data[s-i]*data1[s-i];
		data[s-i]=0.0;
	}

	//Compute inverse fast fourier transformation
	gsl_fft_halfcomplex_radix2_inverse (&(*data.begin()), 1, s);

	//Create output
	DoubleReal shift=0.0;
	for (i = 0; i <= s/2 ; ++i)
	{
		if (data[i] <= 0 || (gauss_fitting && gauss_mean_+shift > mz_max) )
		{
			shift+=stepwidth_;
			continue;
		}
		Peak1D peak;
		peak.setMZ(shift);
		peak.setIntensity(data[i]);
		output.push_back(peak);
		shift+=stepwidth_;
	}
}
}

