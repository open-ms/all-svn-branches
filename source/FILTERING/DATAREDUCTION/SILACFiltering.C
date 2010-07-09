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
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <iostream>

namespace OpenMS
{
DoubleReal SILACFiltering::mz_stepwidth=0;
DoubleReal SILACFiltering::intensity_cutoff=0;
gsl_interp_accel* SILACFiltering::acc_lin=0;
gsl_interp_accel* SILACFiltering::acc_spl=0;
gsl_spline* SILACFiltering::spline_lin=0;
gsl_spline* SILACFiltering::spline_spl=0;
Int SILACFiltering::feature_id=0;


SILACFiltering::SILACFiltering(MSExperiment<Peak1D>& exp_,DoubleReal mz_stepwidth_,DoubleReal intensity_cutoff_) : exp(exp_)
{
	mz_stepwidth=mz_stepwidth_;
	intensity_cutoff=intensity_cutoff_;
}

void SILACFiltering::addFilter(SILACFilter& filter) {
	filters.insert(&filter);
}

SILACFiltering::~SILACFiltering() {
//	for (std::set<SILACFilter*>::iterator it=filters.begin();it!=filters.end();++it)
//	{
//		delete(*it);
//	}
}

void SILACFiltering::blockValue(DoubleReal value,SILACFilter* source)
{
	for (std::set<SILACFilter*>::iterator it=filters.begin();it!=filters.end();++it)
	{
		if (*it!=source)
			(*it)->blockValue(value);
	}
}

bool SILACFiltering::filterPtrCompare::operator ()(SILACFilter* a, SILACFilter* b) const
{
	return (a->getEnvelopeDistanceLightHeavy() < b->getEnvelopeDistanceLightHeavy());
}

void SILACFiltering::filterDataPoints()
{
	startProgress(0,exp.size(),"filtering raw data");



	std::vector<DataPoint> data;

	for (MSExperiment<Peak1D>::Iterator rt_it=exp.begin(); rt_it!=exp.end();++rt_it)
	{
		setProgress(rt_it-exp.begin());
		Size number_data_points = rt_it->size();
//		std::cout << peak_it->getRT() << " " << rt_it->getRT() << std::endl;
		// spectra with less than 10 data points are being ignored
		if (number_data_points>=10) { //filter MS1 spectra (
			// read one OpenMS spectrum into GSL structure
			std::vector<DoubleReal> mz_vec;
			std::vector<DoubleReal> intensity_vec;
			Int j = 0;
			DoubleReal last_mz=rt_it->begin()->getMZ();
			for (MSSpectrum<>::Iterator mz_it=rt_it->begin(); mz_it!=rt_it->end(); ++mz_it)
			{
				if (mz_it->getMZ() > last_mz+2*mz_stepwidth)
				{
					DoubleReal stepwidth=(mz_it->getMZ()-last_mz)/2;
					for (DoubleReal act_mz=last_mz+2*mz_stepwidth; act_mz < mz_it->getMZ()-2*mz_stepwidth; act_mz+=stepwidth)
					{
						mz_vec.push_back(act_mz);
						intensity_vec.push_back(0.0);
					}
				}
				mz_vec.push_back(mz_it->getMZ());
				intensity_vec.push_back(mz_it->getIntensity());
				last_mz=mz_it->getMZ();
			}
//			for (std::vector<DoubleReal>::iterator it=mz_vec.begin();it!=mz_vec.end();++it)
//			{
//				std::cout << *it << std::endl;
//			}

			acc_lin = gsl_interp_accel_alloc();
			spline_lin = gsl_spline_alloc(gsl_interp_linear, mz_vec.size());
			gsl_spline_init(spline_lin, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());
			// spline interpolation
			// used for exact ratio calculation (more accurate when real peak pairs are present)
			acc_spl = gsl_interp_accel_alloc();
			spline_spl = gsl_spline_alloc(gsl_interp_cspline, mz_vec.size());
			gsl_spline_init(spline_spl, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

			for (std::set<SILACFilter*>::iterator filter_it=filters.begin();filter_it!=filters.end();++filter_it)
			{
				SILACFilter* filter_ptr=*filter_it;
				filter_ptr->reset();
			}
			DoubleReal mz_min = mz_vec[0];
			DoubleReal mz_max = mz_vec[mz_vec.size()-9];
			DoubleReal rt=rt_it->getRT();
			for (DoubleReal mz=mz_min; mz<mz_max; mz+=mz_stepwidth)
			{
				if (gsl_spline_eval (spline_spl, mz, acc_spl) < 0.0)
					continue;
				for (std::set<SILACFilter*>::iterator filter_it=filters.begin();filter_it!=filters.end();++filter_it)
				{
					SILACFilter* filter_ptr=*filter_it;
					{
						if (filter_ptr->isFeature(rt,mz))
						{
							const std::vector<DoubleReal>& peak_values=filter_ptr->getPeakValues();
							for (std::vector<DoubleReal>::const_iterator peak_it=peak_values.begin();peak_it!=peak_values.end();++peak_it)
							{
								blockValue(mz,filter_ptr);
							}
							++feature_id;
						}
					}
				}
			}
			gsl_spline_free(spline_lin);
			gsl_interp_accel_free(acc_lin);
			gsl_spline_free(spline_spl);
			gsl_interp_accel_free(acc_spl);
		}
	}
	endProgress();
}

}
