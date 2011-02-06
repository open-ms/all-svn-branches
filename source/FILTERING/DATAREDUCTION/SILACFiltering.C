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

#include <iostream>

using namespace std;

namespace OpenMS
{
  DoubleReal SILACFiltering::mz_stepwidth = 0;
  DoubleReal SILACFiltering::intensity_cutoff = 0;
  DoubleReal SILACFiltering::intensity_correlation = 0;
  bool SILACFiltering::allow_missing_peaks = true;
  gsl_interp_accel* SILACFiltering::current_aki = 0;
  gsl_interp_accel* SILACFiltering::current_spl = 0;
  gsl_spline* SILACFiltering::spline_aki = 0;
  gsl_spline* SILACFiltering::spline_spl = 0;
  Int SILACFiltering::feature_id = 0;
  DoubleReal SILACFiltering::mz_min = 0;

	SILACFiltering::SILACFiltering(MSExperiment<Peak1D>& exp_, DoubleReal mz_stepwidth_, DoubleReal intensity_cutoff_, DoubleReal intensity_correlation_, bool allow_missing_peaks_) : exp(exp_)
	{
		mz_stepwidth = mz_stepwidth_;
		intensity_cutoff = intensity_cutoff_;
		intensity_correlation = intensity_correlation_;
		allow_missing_peaks = allow_missing_peaks_;
	}

	void SILACFiltering::addFilter(SILACFilter& filter)
	{
		filters.push_back(&filter);
	}

	SILACFiltering::~SILACFiltering()
	{

	}

	void SILACFiltering::filterDataPoints()
	{
    startProgress(0, exp.size(), "filtering raw data");

    vector<DataPoint> data;
    vector<BlacklistEntry> blacklist;     // create global blacklist

    mz_min = exp.getMinMZ();      // get lowest m/z value of the experiment

    // Iterate over all filters
    for (list<SILACFilter*>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
    {
      // Iterate over all spectra of the experiment (iterate over rt)
      for (MSExperiment<Peak1D>::Iterator rt_it = exp.begin(); rt_it != exp.end(); ++rt_it)
      {
        previous_entries.clear();     // clear vector of blacklist entries from previous m/z position

        // set progress
        // calculate with progress for the current rt run and progress for the filter run, each scaled by total numbers of filters
        setProgress((rt_it - exp.begin()) / filters.size() + distance(filters.begin(), filter_it) * exp.size() / filters.size());

        Size number_data_points = rt_it->size();    // number of (m/z, intensity) data points in this spectrum

        DoubleReal rt = rt_it->getRT();    // retention time of this spectrum

        // spectra with less than 10 data points are being ignored
        if (number_data_points >= 10)
        {
          // filter MS1 spectra
          // read one spectrum into GSL structure
          vector<DoubleReal> mz_vec;
          vector<DoubleReal> intensity_vec;
          mz_min = rt_it->begin()->getMZ();
          DoubleReal last_mz = rt_it->begin()->getMZ();

          // INTERPOLATION (Akima and Spline interpolation in order to have intensities at any m/z.)
          // Fill intensity and m/z vector for interpolation. Add zeros in the area with no data points to improve cubic spline fit
          for (MSSpectrum<>::Iterator mz_it = rt_it->begin(); mz_it != rt_it->end(); ++mz_it)
          {
            if (mz_it->getMZ() > last_mz + 2 * mz_stepwidth) // If the mz gap is rather larger, fill in zeros. These addtional Stützstellen improve interpolation where no signal (i.e. data points) is.
            {
              for (DoubleReal current_mz = last_mz + 2 * mz_stepwidth; current_mz < mz_it->getMZ() - 2 * mz_stepwidth; current_mz += mz_stepwidth)
              {
                mz_vec.push_back(current_mz);
                intensity_vec.push_back(0.0);
              }
            }
            mz_vec.push_back(mz_it->getMZ());
            intensity_vec.push_back(mz_it->getIntensity());
            last_mz = mz_it->getMZ();
          }

          // akima interpolation, returns 0 in regions with no raw data points
          current_aki = gsl_interp_accel_alloc();
          spline_aki = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
          gsl_spline_init(spline_aki, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

          // spline interpolation, used for exact ratio calculation (more accurate when real peak pairs are present)
          current_spl = gsl_interp_accel_alloc();
          spline_spl = gsl_spline_alloc(gsl_interp_cspline, mz_vec.size());
          gsl_spline_init(spline_spl, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

          MSSpectrum<>::Iterator mz_it = rt_it->begin();

          last_mz = mz_it->getMZ();
					++mz_it;

          // Iterate over the spectrum with a step width that is oriented on the raw data point positions (iterate over mz)
          for ( ; mz_it != rt_it->end(); ++mz_it) // iteration correct
					{
            // We do not move with mz_stepwidth over the spline fit, but with about a third of the local mz differences
            for (DoubleReal mz = last_mz; mz < mz_it->getMZ(); mz += (abs(mz_it->getMZ() - last_mz)) / 3)
						{
              //---------------------------------------------------------------
              // BLUNT INTENSITY FILTER (Just check that intensity at current m/z position is above the intensity cutoff)
              //---------------------------------------------------------------

              if (gsl_spline_eval (spline_aki, mz, current_aki) < intensity_cutoff)
							{
                continue;
							}


              //--------------------------------------------------
              // BLACKLIST FILTER (check if current m/z and rt position is blacklisted)
              //--------------------------------------------------

              bool isBlacklisted = false;
              bool isFriend = false;

              // iterate over the blacklist
              for (vector<BlacklistEntry>::iterator blacklist_it = blacklist.begin(); blacklist_it != blacklist.end(); ++blacklist_it)
              {
                if (isBlacklisted == true || isFriend == true)
                  break;     // if there is already a blacklisted position or current and generating filter are friends

                // check if current and potential corresponding positions are blacklisted
                // (i.e. if positions are inside m/z and rt range of one BlacklistEntry)
                else
                {
                  Int numberOfPeptides = (*filter_it)->getNumberOfPeptides();     // number of labelled peptides + 1 for cuurent filter
                  Int isotopes_per_peptide = (*filter_it)->getIsotopesPerPeptide();     // get number of isotopic peaks per peptide for current filter
                  DoubleReal isotope_distance = (*filter_it)->getIsotopeDistance();     // get distance between isotopic peaks of a peptide in [Th] for current filter

                  // iterate over peptides
                  for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
                  {
                    if (isBlacklisted == true || isFriend == true)
                      break;          // if there is already a blacklisted position or current and generating filter are friends

                    else
                    {
                      DoubleReal mz_peptide_separation = (*filter_it)->getMzPeptideSeparations()[peptide];      // get m/z shift for next potential peak

                      // iterate over isotopes
                      for (Int isotope = 0; isotope < isotopes_per_peptide; isotope++)
                      {
                        if (isBlacklisted == true || isFriend == true)
                          break;      // if there is already a blacklisted position or current and generating filter are friends

                        else
                        {
                          DoubleReal current_mz_position = mz + mz_peptide_separation + isotope * isotope_distance;     // calculate m/z position for next potential peak

                          // perform check
                          if ((blacklist_it->range).encloses(current_mz_position, rt))
                          {
                            // check for blacklisted potential monoisotopic peak if current and generating filter are friends
                            // (i.e. if current and generating filter only differ in number of isotopic peaks per peptide)
                            if (peptide == 0 && isotope == 0)
                            {
                              // check if current an generating flter are equal in charge and number of mass shifts
                              if (blacklist_it->generatingFilter != NULL && (*filter_it)->getCharge() == (blacklist_it->generatingFilter)->getCharge() && (*filter_it)->getMassSeparationsSize() == (blacklist_it->generatingFilter)->getMassSeparationsSize())
                              {
                                vector<DoubleReal> current_filter_mass_separations = (*filter_it)->getMassSeparations();      // get mass shifts for current filter
                                vector<DoubleReal> generating_filter_mass_separations = (blacklist_it->generatingFilter)->getMassSeparations();     // get mass shifts for generating filter

                                // iterate over mass shifts
                                for (unsigned int i = 0; i < current_filter_mass_separations.size(); i++)
                                {
                                  // check mass shifts
                                  if (current_filter_mass_separations[i] != generating_filter_mass_separations[i])
                                  {
                                    isBlacklisted = true;     // current and generating filter differ in mass shift
                                    break;
                                  }
                                  else
                                  {
                                    isFriend = true;      // current and generating filter are friends
                                    break;
                                  }
                                }
                              }
                            }

                            else
                            {
                              isBlacklisted = true;      // one of the positions is found in the blacklist
                              break;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }


              // Check the other filters only if current m/z and rt position is not blacklisted or if current and generating filter are friends
              if (isBlacklisted == false || isFriend == true)
							{
                if ((*filter_it)->isSILACPattern(rt, mz))      // Check if the mz at the given position is a SILAC pair
								{
                  //--------------------------------------------------
                  // blacklisting
                  //--------------------------------------------------

                  BlacklistEntry next_entry;      // create blacklist entry for current m/z and rt position
                  vector<BlacklistEntry> new_entries;     // vector of blacklist entries for current m/z and rt position and following peaks
                  DRange<2> range;      // create range for current m/z and rt position
                  DRange<2> range_united;   // create range that combines cuurent and previous range
                  bool united = false;
                  vector<BlacklistEntry>::iterator previous_it = previous_entries.begin();      // iterator over blacklist entries for previous m/z position

                  // retrieve peak positions for blacklisting
                  const vector<DoubleReal>& peak_positions = (*filter_it)->getPeakPositions();

                  // iterate over the blacklist
                  for (vector<DoubleReal>::const_iterator peak_positions_it = peak_positions.begin(); peak_positions_it != peak_positions.end(); ++peak_positions_it)
                  {                    
                    // get peak width for corresponding peak
                    DoubleReal peak_width = SILACFilter::getPeakWidth(*peak_positions_it);

                    range.setMinX(*peak_positions_it - 0.8 * peak_width);     // set min m/z position of blacklisted range
                    range.setMaxX(*peak_positions_it + 0.8 * peak_width);     // set max m/z position of blacklisted range
                    range.setMinY(rt - 10);     // set min rt position of blacklisted range
                    range.setMaxY(rt + 10);     // set max rt position of blacklisted range

                    // check if range for current m/z position intersects with range from previous m/z position
                    if (previous_entries.size() != 0 && range.isIntersected(previous_it->range))
                    {
                      range_united = range.united(previous_it->range);      // create new minimal range containing current and previous range
                      next_entry.range = range_united;      // add united range to blacklist entry for current m/z and rt position
                      united = true;      // current and previous range intersect and have been combined
                    }
                    else
                      next_entry.range = range;     // add current range to blacklist entry for current m/z and rt position

                    // set generating filter to current filter for monoisotopic peak and to NULL for following peaks
                    if (peak_positions_it == peak_positions.begin())
                      next_entry.generatingFilter = (*filter_it);     // add generating filter pointer to blacklist entry for current m/z and rt position

                    else
                      next_entry.generatingFilter = NULL;     // add NULL pointer to blacklist entry for following peaks

                    new_entries.push_back(next_entry);     // add pointer of current BlacklistEntry to "new_entries"
                    previous_it++;      // get next previous blacklist entry
                  }

                  // erase blacklist entries from previous m/z position if previous and current range intersect and have been combined
                  if (united == true)
                    blacklist.resize(blacklist.size() - previous_entries.size());

                  blacklist.insert(blacklist.end(), new_entries.begin(), new_entries.end());      // insert blacklist entries for current m/z and rt position and following peaks
                  previous_entries.clear();     // clear vector of  blacklist entries for previous m/z position
                  previous_entries.swap(new_entries);     // set blacklist entries for current m/z and rt position to blacklist entries for previous m/z position
                  new_entries.clear();      // clear vector of blacklist entries for current m/z and rt position

									++feature_id;
								}
							}	
            }

            last_mz = mz_it->getMZ();
					}
				}

        // free the interpolation objects
        gsl_spline_free(spline_aki);
        gsl_interp_accel_free(current_aki);
				gsl_spline_free(spline_spl);
				gsl_interp_accel_free(current_spl);
      }
	  }

    endProgress();
	}
}
