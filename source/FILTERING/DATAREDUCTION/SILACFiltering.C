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
#include <OpenMS/FORMAT/MzMLFile.h>
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
DoubleReal SILACFiltering::mz_min=0;


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
	if (a->envelope_distances.size()!=b->envelope_distances.size())
		return a->envelope_distances.size() > b->envelope_distances.size();

		if (a->charge!=b->charge)
			return a->charge > b->charge;
	std::set<DoubleReal>::iterator envelope_distances_it_a=a->envelope_distances.begin();
	std::set<DoubleReal>::iterator envelope_distances_it_b=b->envelope_distances.begin();
	while(*envelope_distances_it_a == *envelope_distances_it_b)
	{
		++envelope_distances_it_a;
		++envelope_distances_it_b;
	}
	return (*envelope_distances_it_a > *envelope_distances_it_b);
}


void SILACFiltering::filterDataPoints()
{
	startProgress(0,exp.size(),"filtering raw data");

	std::vector<DataPoint> data;

	MSExperiment<Peak1D> exp_out;
	MSExperiment<Peak1D> ratios1;
	MSExperiment<Peak1D> ratios2;
	mz_min=exp.getMinMZ();
	static_cast<ExperimentalSettings&>(exp_out) = exp;
	exp_out.resize(exp.size());
	ratios1.resize(exp.size());
	ratios2.resize(exp.size());

	Size scan_idx=0;
	for (MSExperiment<Peak1D>::Iterator rt_it=exp.begin(); rt_it!=exp.end();++rt_it)
	{
		MSSpectrum<Peak1D>& input=exp[scan_idx];
		MSSpectrum<Peak1D>& output=exp_out[scan_idx];
		output.clear(true);
		output.SpectrumSettings::operator=(input);
		output.MetaInfoInterface::operator=(input);
		output.setRT(input.getRT());
		output.setMSLevel(input.getMSLevel());
		output.setName(input.getName());
		output.setType(SpectrumSettings::PEAKS);

		MSSpectrum<Peak1D>& ratio1=ratios1[scan_idx];
		ratio1.clear(true);
		ratio1.SpectrumSettings::operator=(input);
		ratio1.MetaInfoInterface::operator=(input);
		ratio1.setRT(input.getRT());
		ratio1.setMSLevel(input.getMSLevel());
		ratio1.setName(input.getName());
		ratio1.setType(SpectrumSettings::PEAKS);

		MSSpectrum<Peak1D>& ratio2=ratios2[scan_idx];
		ratio2.clear(true);
		ratio2.SpectrumSettings::operator=(input);
		ratio2.MetaInfoInterface::operator=(input);
		ratio2.setRT(input.getRT());
		ratio2.setMSLevel(input.getMSLevel());
		ratio2.setName(input.getName());
		ratio2.setType(SpectrumSettings::PEAKS);

		setProgress(rt_it-exp.begin());
		Size number_data_points = rt_it->size();
//		std::cout << peak_it->getRT() << " " << rt_it->getRT() << std::endl;
		// spectra with less than 10 data points are being ignored
		if (number_data_points>=10) { //filter MS1 spectra (
			// read one OpenMS spectrum into GSL structure
			std::vector<DoubleReal> mz_vec;
			std::vector<DoubleReal> intensity_vec;
			mz_min=rt_it->begin()->getMZ();
			DoubleReal last_mz=rt_it->begin()->getMZ();
			for (MSSpectrum<>::Iterator mz_it=rt_it->begin(); mz_it!=rt_it->end(); ++mz_it)
			{
				if (mz_it->getMZ() > last_mz+2*mz_stepwidth)
				{
					for (DoubleReal act_mz=last_mz+2*mz_stepwidth; act_mz < mz_it->getMZ()-2*mz_stepwidth; act_mz+=mz_stepwidth)
					{
						mz_vec.push_back(act_mz);
						intensity_vec.push_back(0.0);
					}
				}
				mz_vec.push_back(mz_it->getMZ());
				intensity_vec.push_back(mz_it->getIntensity());
				last_mz=mz_it->getMZ();
			}
			acc_lin = gsl_interp_accel_alloc();
			spline_lin = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
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
			MSSpectrum<>::Iterator mz_it=rt_it->begin();
			last_mz=mz_it->getMZ();
			++mz_it;
			for ( ;mz_it!=rt_it->end(); ++mz_it)
			{

				for (DoubleReal mz=last_mz; mz<mz_it->getMZ() && std::abs(last_mz-mz_it->getMZ())<3*mz_stepwidth ;mz+=((last_mz+mz_it->getMZ())/2)-last_mz)
				{
					if (gsl_spline_eval (spline_lin, mz, acc_lin) <= 0.0)
						continue;
					for (std::set<SILACFilter*>::iterator filter_it=filters.begin();filter_it!=filters.end();++filter_it)
					{
						SILACFilter* filter_ptr=*filter_it;
						if (filter_ptr->isFeature(rt,mz))
						{
							const std::vector<DoubleReal>& peak_positions=filter_ptr->getPeakPositions();
							for (std::vector<DoubleReal>::const_iterator peak_it=peak_positions.begin();peak_it!=peak_positions.end();++peak_it)
							{
								blockValue(*peak_it,filter_ptr);
							}
//							if (filter_it!=filters.begin())
								output.push_back(filter_ptr->peak);
								ratio1.push_back(filter_ptr->peak_ratio1);
								ratio2.push_back(filter_ptr->peak_ratio2);
							++feature_id;
						}
//						std::cout << "---------" << std::endl;
					}
				}
				last_mz=mz_it->getMZ();
			}
			gsl_spline_free(spline_lin);
			gsl_interp_accel_free(acc_lin);
			gsl_spline_free(spline_spl);
			gsl_interp_accel_free(acc_spl);
		}
		exp_out[scan_idx]=output;
		ratios1[scan_idx]=ratio1;
		ratios2[scan_idx]=ratio2;
		++scan_idx;
	}
	MzMLFile file1;
	file1.store("/home/steffen/Studium/Master/demos/nijmegen/temp/out.mzML",exp_out);
	MzMLFile file2;
	file2.store("/home/steffen/Studium/Master/demos/nijmegen/temp/ratios1_2.mzML",ratios1);
	MzMLFile file3;
	file3.store("/home/steffen/Studium/Master/demos/nijmegen/temp/ratios2_3.mzML",ratios2);
	endProgress();
}

}
