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
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut $
// --------------------------------------------------------------------------
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#ifdef _OPENMP
#ifdef OPENMS_WINDOWSPLATFORM
#include <omp.h>
#endif
#endif

#include <iostream>

using namespace std;

namespace OpenMS
{
DoubleReal SILACFiltering::mz_stepwidth=0;
DoubleReal SILACFiltering::intensity_cutoff=0;
gsl_interp_accel* SILACFiltering::current_lin=0;
gsl_interp_accel* SILACFiltering::current_spl=0;
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
	filters.push_back(&filter);
}

SILACFiltering::~SILACFiltering() {

}

void SILACFiltering::filterDataPoints()
{
	startProgress(0,exp.size(),"filtering raw data");

	std::vector<DataPoint> data;
	std::list<BlacklistEntry> blacklist;
	//Find out lowest m/z value
	mz_min=exp.getMinMZ();
	
  //Iterate over all spectra of the experiment
  for (MSExperiment<Peak1D>::Iterator rt_it=exp.begin(); rt_it!=exp.end();++rt_it)
	{
    setProgress(rt_it-exp.begin());
	Size number_data_points = rt_it->size();    // number of (m/z, intensity) data points in this spectrum
		
	DoubleReal rt=rt_it->getRT();    // retention time of this spectrum
		
	//kill archaic blacklist entries
	std::list<BlacklistEntry> tempBlacklist;     // instead of erasing entries in blacklist, benerate new tempBlacklist and copy it then back to blacklist
	for (std::list<BlacklistEntry>::iterator blacklist_it=blacklist.begin(); blacklist_it!=blacklist.end(); ++blacklist_it)
	{
		// if the blacklisting is older than 10 s drop it from the list.
		if ( abs(rt - (*blacklist_it).rtInitial) < 10) tempBlacklist.push_back(*blacklist_it);
	}
	blacklist.resize(tempBlacklist.size());
	copy(tempBlacklist.begin(), tempBlacklist.end(), blacklist.begin());		
		
	// spectra with less than 10 data points are being ignored		
    if (number_data_points>=10)
    {
			// filter MS1 spectra
			// read one spectrum into GSL structure
			std::vector<DoubleReal> mz_vec;
			std::vector<DoubleReal> intensity_vec;
			mz_min=rt_it->begin()->getMZ();
			DoubleReal last_mz=rt_it->begin()->getMZ();
		
			// INTERPOLATION (Akima and Spline interpolation in order to have intensities at any m/z.)
			//Fill intensity and m/z vector for interpolation. Add zeros in the area with no data points to improve cubic spline fit
			for (MSSpectrum<>::Iterator mz_it=rt_it->begin(); mz_it!=rt_it->end(); ++mz_it)
			{
				if (mz_it->getMZ() > last_mz+2*mz_stepwidth) // If the mz gap is rather larger, fill in zeros. These addtional Stützstellen improve interpolation where no signal (i.e. data points) is.
				{
					for (DoubleReal current_mz=last_mz+2*mz_stepwidth; current_mz < mz_it->getMZ()-2*mz_stepwidth; current_mz+=mz_stepwidth)
					{
						mz_vec.push_back(current_mz);
						intensity_vec.push_back(0.0);
					}
				}
				mz_vec.push_back(mz_it->getMZ());
				intensity_vec.push_back(mz_it->getIntensity());
				last_mz=mz_it->getMZ();
			}
			
			//akima interpolation,    returns 0 in regions with no raw data points
			current_lin = gsl_interp_accel_alloc();
			spline_lin = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
			gsl_spline_init(spline_lin, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

			// spline interpolation,    used for exact ratio calculation (more accurate when real peak pairs are present)
			current_spl = gsl_interp_accel_alloc();
			spline_spl = gsl_spline_alloc(gsl_interp_cspline, mz_vec.size());
			gsl_spline_init(spline_spl, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());


			//Iterate over all filters
			for (std::list<SILACFilter*>::iterator filter_it=filters.begin(); filter_it!=filters.end(); ++filter_it)
			{
				// debug output: Are all filters in the right order?
				//std::cout << "filter: isotopes per peptide = " << (*filter_it)->getIsotopesPerPeptide() << " charge = " << (*filter_it)->getCharge() << std::endl;
				
				MSSpectrum<>::Iterator mz_it=rt_it->begin();
				
				last_mz=mz_it->getMZ();
				++mz_it;
				//Iterate over the spectrum with a step width that is oriented on the raw data point positions				
				for ( ;mz_it!=rt_it->end(); ++mz_it) // iteration correct
				{
					// loop by Steffen Sass; point of the two-tier-loop: Iterate where raw data points are, not in empty space => (1) better run time (2) less noise, spline fit has 'viel Phantasie' in regions without data points
					//Choose half of the data point distances as stepwidth to take interpolated intensities between the data points into account
					//for (DoubleReal mz=last_mz; mz < mz_it->getMZ() && std::abs(last_mz-mz_it->getMZ()) < 3 * mz_stepwidth ;mz+=((last_mz+mz_it->getMZ())/2)-last_mz)

					// We do not move with mz_stepwidth over the spline fit, but with about a third of the local mz differences
					for (DoubleReal mz=last_mz; mz < mz_it->getMZ(); mz+=(std::abs(mz_it->getMZ() - last_mz))/3)
					{
						if (gsl_spline_eval (spline_lin, mz, current_lin) <= 0.0)
							continue;

						//Check if m/z position is blacklisted
						bool isBlacklisted = false;
						for (std::list<BlacklistEntry>::iterator blacklist_it = blacklist.begin(); blacklist_it!=blacklist.end(); ++blacklist_it)
						{
							if ((mz > blacklist_it->mzBlack_min) && (mz < blacklist_it->mzBlack_max) && blacklist_it->generatingFilter != (*filter_it))
							{
								isBlacklisted = true;
								break;
							}
						}
						
						if (isBlacklisted == false)   //Check the other filters only if m/z is not blacklisted
						{
							if ((*filter_it)->isPair(rt,mz))   //Check if the mz at the given position is a SILAC pair
							{
								//Retrieve peak positions for blacklisting
								const std::vector<DoubleReal>& peak_positions=(*filter_it)->getPeakPositions();
								std::vector<DoubleReal>::const_iterator peak_positions_it = peak_positions.begin();
								DoubleReal peak_width = SILACFilter::getPeakWidth(*peak_positions_it);

								// filling the blacklist
								BlacklistEntry next_entry;
								
								next_entry.mzBlack_min = mz - 0.8 * peak_width;
								next_entry.mzBlack_max = mz + 0.8 * peak_width;
								next_entry.rtInitial = rt;
								next_entry.generatingFilter = (*filter_it);  // Filter generating the blacklisting.
								blacklist.push_back(next_entry);
								
								peak_positions_it++;
						
								for (; peak_positions_it != peak_positions.end(); peak_positions_it++)
								{
									
									DoubleReal peak_width = SILACFilter::getPeakWidth(*peak_positions_it);
									
									next_entry.mzBlack_min = *peak_positions_it - 0.8 * peak_width;
									next_entry.mzBlack_max = *peak_positions_it + 0.8 * peak_width;
									next_entry.rtInitial = rt;
									next_entry.generatingFilter = NULL;
									blacklist.push_back(next_entry);
								}

								++feature_id;
							}
						}	
					}			
					last_mz=mz_it->getMZ();
				}

			}
			//Clear the interpolations
			gsl_spline_free(spline_lin);
			gsl_interp_accel_free(current_lin);
			gsl_spline_free(spline_spl);
			gsl_interp_accel_free(current_spl);
		}

  }
	endProgress();
}



}
