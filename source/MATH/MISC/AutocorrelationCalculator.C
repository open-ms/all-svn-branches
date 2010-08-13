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


#include <OpenMS/MATH/MISC/AutocorrelationCalculator.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
 #include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>



namespace OpenMS
{
AutocorrelationCalculator::AutocorrelationCalculator (DoubleReal stepwidth) : stepwidth_(stepwidth)
{

}

AutocorrelationCalculator::~AutocorrelationCalculator(){
	// TODO Auto-generated destructor stub
}

void AutocorrelationCalculator::calculate(const MSExperiment<Peak1D>& input, MSExperiment<Peak1D>& output)
{
	mz_min=input.getMinMZ();
	mz_max=input.getMaxMZ();
	max_intensity=input.getMaxInt();
	// make sure that output is clear
	output.clear(true);

	// copy experimental settings
	static_cast<ExperimentalSettings&>(output) = input;

	// resize output with respect to input
	output.resize(input.size()-mz_min);

	Size progress = 0;

	startProgress(0,input.size(),"Autocorrelation calculation");
	for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
	{
		if (input[scan_idx].getMSLevel() != 1)
		{
			output[scan_idx] = input[scan_idx];
		}
		else
		{
			analyzeSpectrum(input[scan_idx], output[scan_idx]);
		}
		setProgress(++progress);
	}
	endProgress();
}

struct f_params{ gsl_spline* interpolation; gsl_interp_accel* acc_spl; DoubleReal shift; DoubleReal max_intensity; DoubleReal max_mz;};

DoubleReal f (DoubleReal x, void * params)
{
	struct f_params* pList = (struct f_params *) params;
	gsl_spline* interpolation = pList->interpolation;
	gsl_interp_accel* acc_spl = pList->acc_spl;
	DoubleReal shift = pList->shift;
//	DoubleReal max_intensity = pList->max_intensity;
	DoubleReal max_mz = pList->max_mz;
	DoubleReal f_x=gsl_spline_eval (interpolation, x, acc_spl);
	DoubleReal f_x_shifted=gsl_spline_eval (interpolation, x+shift, acc_spl);
	DoubleReal f=0.0;
	if (x+shift <= max_mz)
		f =  f_x*f_x_shifted;
	return (f);
}


void AutocorrelationCalculator::analyzeSpectrum(const MSSpectrum<Peak1D>& input, MSSpectrum<Peak1D>& output)
{
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
//	DoubleReal last_mz=input.begin()->getMZ();
	for (MSSpectrum<Peak1D>::ConstIterator mz_it=input.begin(); mz_it!=input.end(); ++mz_it)
	{
		mz_vec.push_back(mz_it->getMZ());
		intensity_vec.push_back(mz_it->getIntensity());
	}

	gsl_interp_accel* acc_spl = gsl_interp_accel_alloc();
	gsl_spline* spline_spl = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
	gsl_spline_init(spline_spl, &(*mz_vec.begin()), &(*intensity_vec.begin()), mz_vec.size());

	if (input.size() < 6) return;

	DoubleReal max_distance=mz_max-mz_min;
//	std::cout << mz_min << "\t" << mz_max << "\t" << max_distance << std::endl;

	for (DoubleReal shift=0.0;shift<=max_distance;shift+=stepwidth_)
	{
		gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

		gsl_integration_workspace* w = gsl_integration_workspace_alloc (mz_vec.size());
		struct f_params list= {spline_spl,acc_spl,shift,max_intensity,mz_max};

		gsl_function F;
		F.function = &f;
		F.params = &list;

		DoubleReal result, error;

		DoubleReal relerr=0.01;
		int status=1;
		while(status)
		{
			status=gsl_integration_qags (&F, mz_min, mz_max, 0, relerr, mz_vec.size(), w, &result, &error);
			relerr*=1.1;
		}
		Peak1D peak;
		peak.setMZ(shift);
		peak.setIntensity(result);
		if (result > 0)
		output.push_back(peak);
		gsl_set_error_handler(old_handler);
		gsl_integration_workspace_free (w);
	}
	gsl_spline_free(spline_spl);
	gsl_interp_accel_free(acc_spl);
}

}

