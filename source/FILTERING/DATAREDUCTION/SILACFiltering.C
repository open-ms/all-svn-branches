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

}

void SILACFiltering::blockPositions(const std::vector<DoubleReal>& peak_positions,SILACFilter* source)
{
	//Add values to the filter blacklists
	for (std::set<SILACFilter*>::iterator it=filters.begin();it!=filters.end();++it)
	{
		//Iterate over all values
		for (std::vector<DoubleReal>::const_iterator peak_it=peak_positions.begin();peak_it!=peak_positions.end();++peak_it)
		{
			//If the values originate from the cluster itself, the monoisotopic peak is skipped
			if (*it==source && peak_it==peak_positions.begin())
				++peak_it;
			//Add value to the blacklist
			(*it)->blockValue(*peak_it);
		}
	}
}

bool SILACFiltering::filterPtrCompare::operator ()(SILACFilter* a, SILACFilter* b) const
{
	//The filters are ordered by a certain hierarchy
	//First: amount of mass shifts (quadruplets, triplets, doublets, singlets); highest amount first
	if (a->envelope_distances.size()!=b->envelope_distances.size())
		return a->envelope_distances.size() > b->envelope_distances.size();

		//Second: charge; highest charge first
		if (a->charge!=b->charge)
			return a->charge > b->charge;

			//Third: size of mass shifts; lowest mass shift first
			//Mass shifts of both compared filters are iterated until they differ
			std::set<DoubleReal>::iterator envelope_distances_it_a=a->envelope_distances.begin();
			std::set<DoubleReal>::iterator envelope_distances_it_b=b->envelope_distances.begin();
			while(*envelope_distances_it_a == *envelope_distances_it_b)
			{
				++envelope_distances_it_a;
				++envelope_distances_it_b;
			}
			//the different mass shift is compared
			return (*envelope_distances_it_a < *envelope_distances_it_b);
}


void SILACFiltering::filterDataPoints()
{
	startProgress(0,exp.size(),"filtering raw data");

	std::vector<DataPoint> data;
	//Find out lowest m/z value
	mz_min=exp.getMinMZ();
	
	//Iterate over all spectra of the experiment
	for (MSExperiment<Peak1D>::Iterator rt_it=exp.begin(); rt_it!=exp.end();++rt_it)
	{
		
		setProgress(rt_it-exp.begin());
		Size number_data_points = rt_it->size();
		// spectra with less than 10 data points are being ignored
		if (number_data_points>=10) { //filter MS1 spectra (
			// read one OpenMS spectrum into GSL structure
			std::vector<DoubleReal> mz_vec;
			std::vector<DoubleReal> intensity_vec;
			mz_min=rt_it->begin()->getMZ();
			DoubleReal last_mz=rt_it->begin()->getMZ();
			//Fill intensity and m/z vector for interpolation. Add zeros in the area with no data points to improve cubic spline fit
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
			
			//akima interpolation; returns 0 in regions with no raw data points
			acc_lin = gsl_interp_accel_alloc();
			spline_lin = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
			gsl_spline_init(spline_lin, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

			// spline interpolation
			// used for exact ratio calculation (more accurate when real peak pairs are present)
			acc_spl = gsl_interp_accel_alloc();
			spline_spl = gsl_spline_alloc(gsl_interp_cspline, mz_vec.size());
			gsl_spline_init(spline_spl, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

			DoubleReal rt=rt_it->getRT();

			//Iterate over all filters
			for (std::set<SILACFilter*>::iterator filter_it=filters.begin();filter_it!=filters.end();++filter_it)
			{
				// current spectrum: rt_it
				/*std::cout << "RT = " << rt_it->getRT();
				std::cout << "   min MZ = " << rt_it->getMin()[0];
				std::cout << "   max MZ = " << rt_it->getMax()[0] << std::endl;*/
				
				//Extract current spectrum
				MSSpectrum<>::Iterator mz_it=rt_it->begin();
				
				last_mz=mz_it->getMZ();
				++mz_it;
				//Iterate over the spectrum with a step width that is oriented on the raw data point positions				
				for ( ;mz_it!=rt_it->end(); ++mz_it) // iteration correct
				{
					std::cout << "extern: m/z=" << mz_it->getMZ() << std::endl;
					
					// loop by Steffen Sass; point of the two-tier-loop: Iterate where raw data points are, not in empty space => (1) better run time (2) less noise, spline fit has 'viel Phantasie' in regions without data points
					//Choose half of the data point distances as stepwidth to take interpolated intensities between the data points into account
					//for (DoubleReal mz=last_mz; mz < mz_it->getMZ() && std::abs(last_mz-mz_it->getMZ()) < 3 * mz_stepwidth ;mz+=((last_mz+mz_it->getMZ())/2)-last_mz)
					
					// We do not move with mz_stepwidth over the spline fit, but with about a third of the local mz differences
					for (DoubleReal mz=last_mz; mz < mz_it->getMZ(); mz+=(std::abs(mz_it->getMZ() - last_mz))/3)
					{
						std::cout << "    intern m/z = " << mz << std::endl;
						
						if (gsl_spline_eval (spline_lin, mz, acc_lin) <= 0.0)
							continue;

						SILACFilter* filter_ptr=*filter_it;
						//Check if the filter at the given position is a SILAC pair
						if (filter_ptr->isPair(rt,mz))
						{
							//Retrieve peak positions for blacklisting
							const std::vector<DoubleReal>& peak_positions=filter_ptr->getPeakPositions();
							blockPositions(peak_positions,filter_ptr);
							++feature_id;
						}
					}
					last_mz=mz_it->getMZ();
				}

			}
			//Clear the interpolations
			gsl_spline_free(spline_lin);
			gsl_interp_accel_free(acc_lin);
			gsl_spline_free(spline_spl);
			gsl_interp_accel_free(acc_spl);
		}
		//Reset the filters
		for (std::set<SILACFilter*>::iterator filter_it=filters.begin();filter_it!=filters.end();++filter_it)
		{
			SILACFilter* filter_ptr=*filter_it;
			filter_ptr->reset();
		}
	}
	endProgress();
}

}
