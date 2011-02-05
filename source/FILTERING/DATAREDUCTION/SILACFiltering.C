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
  gsl_interp_accel* SILACFiltering::current_lin = 0;
  gsl_interp_accel* SILACFiltering::current_spl = 0;
  gsl_spline* SILACFiltering::spline_lin = 0;
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

    Int old_blacklist_size = 0;

    vector<DataPoint> data;
    vector<BlacklistEntry> blacklist;

    mz_min = exp.getMinMZ();      // find out lowest m/z value

    // Iterate over all filters
    for (list<SILACFilter*>::iterator filter_it = filters.begin(); filter_it != filters.end(); ++filter_it)
    {
      // Iterate over all spectra of the experiment (iterate over rt)
      for (MSExperiment<Peak1D>::Iterator rt_it = exp.begin(); rt_it != exp.end(); ++rt_it)
      {
        previous_entries.clear();     // vector of pointers of previous BlacklistEntry
//        cout << "size of previous after initializing: " << previous_entries.size() << endl;

        bool temp = false;

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
          current_lin = gsl_interp_accel_alloc();
          spline_lin = gsl_spline_alloc(gsl_interp_akima, mz_vec.size());
          gsl_spline_init(spline_lin, &*mz_vec.begin(), &*intensity_vec.begin(), mz_vec.size());

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

              if (gsl_spline_eval (spline_lin, mz, current_lin) < intensity_cutoff)
							{
                continue;
							}


              //--------------------------------------------------
              // BLACKLIST FILTER (check if current m/z and rt position is blacklisted)
              //--------------------------------------------------

              bool isBlacklisted = false;
              bool following_isBlacklisted = false;

              // iterate over the blacklist
              for (vector<BlacklistEntry>::iterator blacklist_it = blacklist.begin(); blacklist_it != blacklist.end(); ++blacklist_it)
              {
                // check if position is blacklisted and if so break
                if (isBlacklisted == true)
                {
                  break;
                }

                // check if position is blacklisted (i.e. position is inside m/z and rt range and if current and generating filter are different)
                else if ((blacklist_it->range).encloses(mz, rt) && (*filter_it) != blacklist_it->generatingFilter)
                {
                  // check if current filter is friend of generating filter
                  // i.e. current and generating filter differ in isotopes_per_peptide
                  // but are equal in charge state and mass shifts

                  // check if genearting filter is not NULL and if size of mass shifts for current and generating filter is equal
                  if (blacklist_it->generatingFilter != NULL && (*filter_it)->getMassSeparationsSize() == (blacklist_it->generatingFilter)->getMassSeparationsSize() && (*filter_it)->getCharge() == (blacklist_it->generatingFilter)->getCharge())
                  {
                    vector<DoubleReal> current_filter_mass_separations = (*filter_it)->getMassSeparations();
                    vector<DoubleReal> generating_filter_mass_separations = (blacklist_it->generatingFilter)->getMassSeparations();

                    // check if mass shifts of current and generating filter are equal
                    for (unsigned int i = 0; i < current_filter_mass_separations.size(); i++)
                    {
                      if (current_filter_mass_separations[i] != generating_filter_mass_separations[i])
                      {
                        isBlacklisted = true;     // if mass shifts of current and generating filter differ in only one entry set "isBlacklisted" to true
                        break;
                      }
                    }
                  }
                  else
                  {
                    isBlacklisted = true;
                    break;
                  }
                }
                else
                {
                  // check if m/z positions for potential following peaks are blacklisted
                  // calculate expected positions of pontential following peaks
                  Int numberOfPeptides = (*filter_it)->getNumberOfPeptides();     // number of labelled peptides + 1 for cuurent filter
                  Int isotopes_per_peptide = (*filter_it)->getIsotopesPerPeptide();     // get number of isotopic peaks per peptide for current filter
                  DoubleReal isotope_distance = (*filter_it)->getIsotopeDistance();     // get distance between isotopic peaks of a peptide in [Th] for current filter

                  // iterate over number of peptides
                  for (Int peptide = 0; peptide <= numberOfPeptides; peptide++)
                  {
                    if (following_isBlacklisted == true)
                    {
                      break;    // check if there is already a blacklisted m/z position for potential following peaks
                    }
                    else
                    {
                      // iterate over number of isotopic peaks per peptide
                      for (Int isotope = 0; isotope < isotopes_per_peptide - 1; isotope++)
                      {
                        if (peptide == 0 && isotope == 0)
                        {
                          continue;   // do not check again if current m/z position is blacklisted
                        }
                        else
                        {
                          DoubleReal mz_peptide_separation = (*filter_it)->getMzPeptideSeparations()[peptide];      // get m/z shift for next potential peak
                          DoubleReal mz_following_peak = mz + mz_peptide_separation + isotope * isotope_distance;     // calculate m/z position for next potential peak

                          // check if next potential peak is blacklisted
                          if ((blacklist_it->range).encloses(mz_following_peak, rt))
                          {
                            isBlacklisted = true;
                            following_isBlacklisted = true;
                            break;
                          }
                        }
                      }
                    }
                  }
                }
              }


              // Check the other filters only if current m/z and rt position is not blacklisted
              if (isBlacklisted == false)
							{
                if ((*filter_it)->isSILACPattern(rt, mz))      // Check if the mz at the given position is a SILAC pair
								{
                  //--------------------------------------------------
                  // blacklisting
                  //--------------------------------------------------

                  BlacklistEntry next_entry;
                  DRange<2> range;
                  DRange<2> range_united;
                  vector<BlacklistEntry> new_entries;

                  vector<BlacklistEntry>::iterator previous_it = previous_entries.begin();

                  bool united = false;

                  // Retrieve peak positions for blacklisting
                  const vector<DoubleReal>& peak_positions = (*filter_it)->getPeakPositions();
//                  cout << "start of blacklisting - size of peak positions: " << peak_positions.size() << " size of previous: " << previous_entries.size() << " size of blacklist: " << blacklist.size() << endl;

                  // filling the blacklist
                  for (vector<DoubleReal>::const_iterator peak_positions_it = peak_positions.begin(); peak_positions_it != peak_positions.end(); ++peak_positions_it)
                  {                    
                    // get peak width for corresponding peak
                    DoubleReal peak_width = SILACFilter::getPeakWidth(*peak_positions_it);

                    range.setMinX(*peak_positions_it - 0.8 * peak_width);
                    range.setMaxX(*peak_positions_it + 0.8 * peak_width);
                    range.setMinY(rt - 10);
                    range.setMaxY(rt + 10);

/*                    if (previous_entries.size() != 0)
                    {
                      temp = range.isIntersected(previous_it->range);

                      cout << "0 before if intersected - size of previous: " << previous_entries.size() << " temp: " << temp << endl;
                      cout << "range previous - x_min: " << (previous_it->range).minX() << " x_max: " << (previous_it->range).maxX() << " y_min: " << (previous_it->range).minY() << " y_max: " << (previous_it->range).maxY() << endl;
                      cout << "range current - x_min: " << range.minX() << " x_max: " << range.maxX() << " y_min: " << range.minY() << " y_max: " << range.maxY() << endl;
                    }
*/
                    if (previous_entries.size() != 0 && range.isIntersected(previous_it->range))
                    {
                      range_united = range.united(previous_it->range);
                      next_entry.range = range_united;
                      united = true;

//                      cout << "if intersected: " << temp << endl;
                    }
                    else
                    {
                      next_entry.range = range;

//                      cout << "else intersected: " << temp << endl;
                    }

                    // set generating filter to current filter for monoisotopic peak and to NULL for following peaks
                    if (peak_positions_it == peak_positions.begin())
                    {
 //                     cout << "if generating: " << temp << endl;
                      next_entry.generatingFilter = (*filter_it);     // filter generating the  blacklist entry
                    }
                    else
                    {
//                      cout << "else generating: " << temp << endl;
                      next_entry.generatingFilter = NULL;
                    }

                    new_entries.push_back(next_entry);     // add pointer of current BlacklistEntry to "new_entries"
                    previous_it++;
                  }

                  if (united == true)
                  {
                    cout << "before resize - size of blacklist: " << blacklist.size() << endl;
                    blacklist.resize(blacklist.size() - previous_entries.size());
                  }

                  old_blacklist_size += 6;

                  //cout << "after resize - size of blacklist: " << blacklist.size() << " size of blacklist_temp: " << blacklist_temp.size() << endl;
                  blacklist.insert(blacklist.end(), new_entries.begin(), new_entries.end());
                  cout << "after insert - size of united blacklist: " << blacklist.size() << " size \"old\" blacklist: " << old_blacklist_size << endl;



                  //cout << "before swap size of previous: " << previous_entries.size() << " size of new: " << new_entries.size() << endl;
                  previous_entries.clear();
                  previous_entries.swap(new_entries);
                  new_entries.clear();

//                  cout << "after swap size of previous: " << previous_entries.size() << " size of new: " << new_entries.size() << endl;
//                  cout << "size of blacklist: " << blacklist.size() << endl;


									++feature_id;
								}
							}	
						}			
            last_mz = mz_it->getMZ();
					}
				}

        // Clear the interpolations
				gsl_spline_free(spline_lin);
				gsl_interp_accel_free(current_lin);
				gsl_spline_free(spline_spl);
				gsl_interp_accel_free(current_spl);
			}
	  }

		endProgress();
	}
}
