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
#include <OpenMS/CONCEPT/Constants.h>
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

std::vector<DoubleReal> CrossCorrelationCalculator::calculate(const MSExperiment<Peak1D>& input, MSExperiment<Peak1D>& output,Size spectrum_selection_id)
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
	std::vector<DoubleReal> data;
	for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
	{
		if (input[scan_idx].getMSLevel() != 1)
		{
			output[scan_idx] = input[scan_idx];
		}
		else if (scan_idx==spectrum_selection_id)
		{
			std::vector<DoubleReal> act_data=analyzeSpectrum(input[scan_idx], output[scan_idx]);
			data.insert(data.end(),act_data.begin(),act_data.end());
		}
		else
		{
			analyzeSpectrum(input[scan_idx], output[scan_idx]);
		}
		setProgress(++progress);
	}
	endProgress();
	return data;
}

std::vector<DoubleReal> CrossCorrelationCalculator::analyzeSpectrum(const MSSpectrum<Peak1D>& input, MSSpectrum<Peak1D>& output, bool gauss_fitting)
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

	DoubleReal window_size=10;
	DoubleReal stepwidth=stepwidth_;
	DoubleReal act_mz=gauss_mean_;

	//n must be a power of two for gsl Fourier transformation; take next higher size for n, which is a power of two
	Size vector_size = pow(2,(ceil(log2(window_size/stepwidth))));

	//Create a vector containing interpolated values with a spacing of "stepwidth"
	std::vector<DoubleReal> data(vector_size,0);
	Size i=0;
	DoubleReal starting_offset=std::min(act_mz-mz_min,3*gauss_sigma_);
	for (DoubleReal x=act_mz-starting_offset;x<=act_mz+window_size;x+=stepwidth)
	{
		data[i] = gsl_spline_eval (spline_spl, x, acc_spl);
		++i;
	}

	std::vector<DoubleReal> gauss_fitted_data(vector_size,0);
	DoubleReal gauss_normalization_factor=gauss_sigma_*sqrt(2*Constants::PI);
	i=0;
	for (DoubleReal x=-starting_offset;x<=window_size;x+=stepwidth)
	{
		gauss_fitted_data[i] = data[i]*gsl_ran_gaussian_pdf (x, gauss_sigma_)*gauss_normalization_factor;
		++i;
	}

	//Fourier transformation of the values
	gsl_fft_real_radix2_transform (&(*data.begin()), 1, vector_size);
	gsl_fft_real_radix2_transform (&(*gauss_fitted_data.begin()), 1, vector_size);

	//Multiply the fourier transformed complex values with the complex conjugate
	//Have a look at the GNU Scientific Library reference manual for a description of the data structure
	data[0]=data[0]*gauss_fitted_data[0];
	data[vector_size/2]=data[vector_size/2]*gauss_fitted_data[vector_size/2];
	for (i = 1; i <= vector_size/2; ++i)
	{
		data[i]=data[i]*gauss_fitted_data[i]+data[vector_size-i]*gauss_fitted_data[vector_size-i];
		data[vector_size-i]=0.0;
	}

	//Compute inverse fourier transformation
	gsl_fft_halfcomplex_radix2_inverse (&(*data.begin()), 1, vector_size);

	DoubleReal intensity_normalization_factor=gsl_spline_eval_integ (spline_spl,act_mz-starting_offset,act_mz+window_size, acc_spl);

	//Create output
	DoubleReal shift=0.0;
	for (i = 0; i <= vector_size/2 ; ++i)
	{
		if (data[i] <= 0)
		{
			shift+=stepwidth;
			continue;
		}
		Peak1D peak;
		peak.setMZ(shift);
		peak.setIntensity(data[i]/intensity_normalization_factor);
		output.push_back(peak);
		shift+=stepwidth;
	}
	gsl_spline_free(spline_spl);
	gsl_interp_accel_free(acc_spl);
	return data;
}

std::vector<DoubleReal> CrossCorrelationCalculator::getExactPositions(const std::vector<DoubleReal>& data,std::vector<DoubleReal> positions, DoubleReal tolerance)
{
	std::vector<DoubleReal> x(data.size()/2,0);
	std::vector<DoubleReal> y(data.size()/2,0);
	DoubleReal shift=0.0;
	for (Size i=0;i<data.size()/2;++i)
	{
		x[i]=shift;
		y[i]=data[i];
		shift+=stepwidth_;
	}
	gsl_interp_accel* acc_spl = gsl_interp_accel_alloc();
	gsl_spline* spline_spl = gsl_spline_alloc(gsl_interp_cspline, x.size());
	gsl_spline_init(spline_spl, &(*x.begin()), &(*y.begin()), x.size());

	std::vector<DoubleReal> exact_positions;

	for (std::vector<DoubleReal>::iterator vec_it=positions.begin();vec_it!=positions.end();++vec_it)
	{
		DoubleReal position=*vec_it;
		DoubleReal exact_position=0.0;
		DoubleReal max_intensity=0.0;
		for (DoubleReal act_position=position-tolerance;act_position<=position+tolerance;act_position+=stepwidth_)
		{
			DoubleReal act_intensity=gsl_spline_eval (spline_spl, act_position, acc_spl);
			if (act_intensity > max_intensity)
			{
				exact_position=act_position;
				max_intensity=act_intensity;
			}
		}
		exact_positions.push_back(exact_position);
	}
	return exact_positions;
}

}

