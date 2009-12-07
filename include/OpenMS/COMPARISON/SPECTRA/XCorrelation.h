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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_XCORRELATION_H
#define OPENMS_COMPARISON_SPECTRA_XCORRELATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <vector>
#include <list>
#include <limits>
#include <utility>
#include <algorithm>

#define XCORR_DEBUG
#undef  XCORR_DEBUG

namespace OpenMS
{
	class Peak1D;

	/**
		@brief compares the correlation of two spectra

		@htmlinclude OpenMS_XCorrelation.parameters

		@ingroup SpectraComparison
	*/
	template <typename PeakT = Peak1D>
	class XCorrelation
		: public DefaultParamHandler
	{

		public:

		typedef PeakT PeakType;
		typedef MSSpectrum<PeakType> SpectrumType;
		typedef typename MSSpectrum<PeakType>::ConstIterator ConstSpectrumIterator;
		typedef typename MSSpectrum<PeakType>::Iterator SpectrumIterator;

		// @name Constructors and Destructors
		// @{
		/// default constructor
		XCorrelation()
			: DefaultParamHandler("XCorrelation")
		{
			defaults_.setValue("peak_tolerance", 0.3, "Defines the absolut (in Da) peak tolerance");
			defaults_.setValue("parentmass_tolerance", 3.0, "Defines the absolut (in Da) parent mass tolerance");
			defaults_.setValue("min_shift", 0.0, "Defines the minimal absolut (in Da) shift between the two spectra");
			defaults_.setValue("min_dist", 57.0214637230, "Defines the minimal distance (in Da) between the two peaks of one sort that may be connected in a sparse matching");
			defaults_.setValue("correlation_scoring", "intensity", "If intensity, correlation scores on basis of matched intensity values, if matchnumber correlation scores solely on basis of the number of matches");
			defaults_.setValidStrings("correlation_scoring", StringList::create("intensity,matchnumber"));
			defaultsToParam_();
		}

		/// copy constructor
		XCorrelation(const XCorrelation& source)
			: DefaultParamHandler(source)
		{
		}

		/// destructor
		~XCorrelation()
		{
		}

		/// assignment operator
		XCorrelation& operator = (const XCorrelation& source)
		{
			if (this != &source)
			{
				DefaultParamHandler::operator = (source);
			}
			return *this;
		}

		// @}


		/**
			@brief Method to calculate best sum of intensities from matching by mapping s1 onto s2 (bigger is better)

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		DoubleReal bestMatchIntensity(MSSpectrum<PeakType>& s1, MSSpectrum<PeakType>& s2) const
		{
			DoubleReal peak_tolerance = (double)param_.getValue("peak_tolerance");
			DoubleReal shiftwindow = (double)param_.getValue("parentmass_tolerance");

			//~ DoubleReal stepsize = (double)param_.getValue("scan_resolution");
			DoubleReal stepsize = 2*peak_tolerance;
			/// @improvement for prms score initialization with 0 is no good use numeric_limits<int>::min() o.s.
			/// @improvement depend stepsize from resolution( stepsize = resolution??)


			DoubleReal best_score(0);
			for(DoubleReal shift = -shiftwindow; shift <= shiftwindow; shift+=stepsize)
			{
				DoubleReal cur_score(0);

				std::vector<std::pair<Size, Size> > matches_unshift_all;
				for(Size i = 0; i < s1.size(); ++i) /// @improvement not that i and i+1 matches still overlap ...
				{
					DoubleReal mz = s1[i].getMZ() + shift;
					ConstSpectrumIterator start(s2.MZBegin(mz-peak_tolerance));
					ConstSpectrumIterator end = (s2.MZEnd(mz+peak_tolerance));
					DoubleReal best_pair_score(0);
					for(ConstSpectrumIterator it = start; it != end; ++it)
					{
						Size j = it - s2.begin();
						DoubleReal cur_pair_score(s1[i].getIntensity()+s2[j].getIntensity());
						if(cur_pair_score>best_pair_score)
						{
							best_pair_score = cur_pair_score;
						}
					}
					cur_score += best_pair_score;
				}

				if(cur_score>best_score)
				{
					best_score = cur_score;
				}
			}
			return best_score;
		}
		///

		/**
			@brief Method to calculate best sum of distances from matching by mapping s1 onto s2 (smaller is better)

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		DoubleReal bestMatchDistance(MSSpectrum<PeakType>& s1, MSSpectrum<PeakType>& s2) const
		{
			DoubleReal peak_tolerance = (double)param_.getValue("peak_tolerance");
			DoubleReal shiftwindow = (double)param_.getValue("parentmass_tolerance");
			DoubleReal stepsize = 2*peak_tolerance;

			DoubleReal best_score(0);
			for(DoubleReal shift = -shiftwindow; shift <= shiftwindow; shift+=stepsize)
			{
				DoubleReal cur_score(0);

				std::vector<std::pair<Size, Size> > matches_unshift_all;
				ConstSpectrumIterator pivot = s2.begin();
				for(Size i = 0; i < s1.size(); ++i) /// @improvement not that i and i+1 matches still overlap ...
				{
					DoubleReal mz = s1[i].getMZ() + shift;
					ConstSpectrumIterator start(s2.MZBegin(pivot,mz-peak_tolerance,s2.end()));
					pivot = start;
					ConstSpectrumIterator end = (s2.MZEnd(pivot,mz+peak_tolerance,s2.end()));
					DoubleReal best_pair_score(0);
					for(ConstSpectrumIterator it = start; it != end; ++it)
					{
						Size j = it - s2.begin();
						DoubleReal cur_pair_score(fabs(s1[i].getMZ()-s2[j].getMZ()));
						if(cur_pair_score<best_pair_score)
						{
							best_pair_score = cur_pair_score;
						}
					}
					cur_score += best_pair_score;
				}

				if(cur_score<best_score)
				{
					best_score = cur_score;
				}
			}
			return best_score;
		}
		///

		/**
			@brief Method to find the maximal long set of sparse matchings by mapping s1 onto s2

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		void maxSparseMatch(const SpectrumType& s1, const SpectrumType& s2, std::list<std::pair<Size,Size> >& best_matches, DoubleReal& shift) const
		{
			DoubleReal peak_tolerance = (double)param_.getValue("peak_tolerance");
			DoubleReal min_dist = (double)param_.getValue("min_dist");

			//~ find best with DP from all matches(mapping s1 to s2)
			std::vector<std::pair<Size, Size> > all_matches;
			ConstSpectrumIterator pivot = s2.begin();
			for(Size i = 0; i < s1.size(); ++i)
			{
				DoubleReal mz = s1[i].getMZ() + shift;
				ConstSpectrumIterator start = (s2.MZBegin(pivot,mz-peak_tolerance,s2.end()));
				pivot = start;
				ConstSpectrumIterator end = (s2.MZEnd(pivot,mz+peak_tolerance/* +0.00001 */,s2.end()));
				for(ConstSpectrumIterator it = start; it != end; ++it)
				{
					Size j = it - s2.begin();
					all_matches.push_back(std::pair<Size,Size>(i,j));
				}
			}

			//~ initialization:
			/// @attention indices here are indices in the dp-table and to get from these to indices in all_matches you have to decrease once;
			std::pair<Size,DoubleReal> best_dp_col(0.0,0.0); //will hold index of the last in the best path so far and the corresponding score
			std::vector< std::pair<Size,DoubleReal> > dp_table(all_matches.size()+1, best_dp_col); // [predecessor, predecessor score], size+1 due to DP initialization in 0
			Size next_too_close = 1;		// Index for the match (in dp_table and all_matches) that is next BEHIND i and too close
			Size max_index = 0;					// Index for best match BEHIND next_too_close with the highest score (already cumulated) in the dp_table

			//~ "recursion":
			for(Size i=1; i< dp_table.size(); ++i)
			{
				while(next_too_close < i and
							s1[all_matches[next_too_close-1].first].getMZ() <= (s1[all_matches[i-1].first].getMZ()-min_dist+2*peak_tolerance) and
							s2[all_matches[next_too_close-1].second].getMZ() <= (s2[all_matches[i-1].second].getMZ()-min_dist+2*peak_tolerance))
				/// @improvement find out if not better only to add one time the peak_tolerance
				{
					// This is executed only when next_too_close is advanced (then there may be a new max_index?)
					if (dp_table[next_too_close].second > dp_table[max_index].second)
					{
						max_index = next_too_close;
					}
					++next_too_close;
				}
				dp_table[i].first = max_index;
				dp_table[i].second = dp_table[max_index].second + s1[all_matches[i-1].first].getIntensity() + s2[all_matches[i-1].second].getIntensity(); //path score with intensities!
				/*debug std::cout << "dp table (" << i << "): " << dp_table[i].first << " | " << dp_table[i].second;*/

				//~ for traceback: best path ends in best_dp_col.first
				if (dp_table[i].second>best_dp_col.second)
				{
					best_dp_col.second=dp_table[i].second;
					best_dp_col.first=i;
					/*debug std::cout << " - new best ";*/
				}
				/*debug std::cout << std::endl;*/
			}

			//~ traceback:
			std::list<Size> best_path;
			while (best_dp_col.first>0)
			{
				best_path.push_back(best_dp_col.first-1);
				best_dp_col.first=dp_table[best_dp_col.first].first;
			}
			while(best_path.size()>0)
			{
				best_matches.push_back(all_matches[best_path.back()]);
				best_path.pop_back();
			}
		}
		///

		/**
			@brief Method to find best correlation of s1 and s2 by shifting around

			@param ...
			@param pm_diff_shift boolean indicates if s2 is considered shifted by pm difference or not
			@return ...

			...
			@see ...
		*/
		void getXCorrelation(SpectrumType& s1, SpectrumType& s2, DoubleReal& best_score1 , DoubleReal& best_score2, DoubleReal& best_shift, std::list<std::pair<Size,Size> >& best_matches, bool pm_diff_shift = false) const
		{
			DoubleReal peak_tolerance = (double)param_.getValue("peak_tolerance");
			DoubleReal parentmass_tolerance = (double)param_.getValue("parentmass_tolerance");
			DoubleReal shift_step = 2 * peak_tolerance;

			//~ can also deal with neg. shifts!
			//~ if(s2.getPrecursors().front().getMZ() < s1.getPrecursors().front().getMZ())
			//~ {
				//~ throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "prerequisite is the s1-precursor mz is not greater than s2-precursor mz");
			//~ }

			DoubleReal pm_s1 = s1.getPrecursors().front().getMZ();
			int c_s1 = s1.getPrecursors().front().getCharge();
			DoubleReal pm_s2 = s2.getPrecursors().front().getMZ();
			int c_s2 = s2.getPrecursors().front().getCharge();
			/// @attention if the precursor charge is unknown, i.e. 0 best guess is its doubly charged
			(c_s1==0)?c_s1=2:c_s1=c_s1;
			(c_s2==0)?c_s2=2:c_s2=c_s2;
			/// @attention singly charged mass difference!
			DoubleReal pm_diff = (pm_s2*c_s2 - (c_s2-1)*Constants::PROTON_MASS_U)-(pm_s1*c_s1 - (c_s1-1)*Constants::PROTON_MASS_U);
			/* debug std::cout << pm_diff << std::endl; */
			if(!pm_diff_shift)
			{
				pm_diff = 0.0;
			}

			if(s1.getPrecursors().front().getCharge()>2 or s2.getPrecursors().front().getCharge()>2)
			{
				throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "prerequisite are spectra from chargestate 2 or lower");
			}

			// reset the correlation values
			best_score1 = 0;
			best_score2 = 0;
			best_shift = std::numeric_limits<double>::min();
			best_matches.clear();
			s1.sortByPosition();
			s2.sortByPosition();

			if((String)param_.getValue("correlation_scoring")=="intensity")
			{
				//~ find matches in unshifted s1 to s2 max sparses matches with DP
				std::list<std::pair<Size, Size> > matches_unshift;
				maxSparseMatch(s1,s2,matches_unshift, pm_diff);
				/*debug std::cout << "matches_unshift.size(): " << matches_unshift.size() << std::endl;*/

				best_matches = matches_unshift;
				//~ calc score
				for(std::list<std::pair<Size,Size> >::iterator it = best_matches.begin(); it != best_matches.end(); ++it)
				{
					best_score1 += s1[it->first].getIntensity();
					best_score2 += s2[it->second].getIntensity();
				}

				//~ shifts only neccessary if spectra reasonable apart or explicitly called:
				if(fabs(pm_diff)>parentmass_tolerance or !pm_diff_shift)
				{
					DoubleReal shift(pm_diff-parentmass_tolerance);
					DoubleReal range(pm_diff+parentmass_tolerance);
					for(; shift<=range; shift+=shift_step)
					{
						DoubleReal close_to_unshift(fabs(fabs(pm_diff) - fabs(shift))+0.00001);
						/// @improvement some situation 0.3 < 0.3 = true pops up if not +0.00001 added
						if(close_to_unshift<peak_tolerance)
						{
							continue;
						}

						/// @improvement add a correlation_scoring method matches that optimizes only the matchnumber (new max_sparse_overlap)
						//~ find max sparse matches in shifted s1 to s2 with DP
						std::list<std::pair<Size, Size> > matches_shift;
						maxSparseMatch(s1,s2,matches_shift,shift);
						/*debug std::cout << "matches_shift.size(): " << matches_shift.size() << " shift: " << shift << std::endl;*/

						DoubleReal score1 = std::numeric_limits<double>::min();
						DoubleReal score2 = std::numeric_limits<double>::min();
						std::list<std::pair<Size,Size> > matches;

						if(matches_shift.size() > 0)
						{
							/*debug for(std::vector<std::pair<Size,Size> >::iterator it = matches_unshift_all.begin(); it != matches_unshift_all.end(); ++it)
											{
												std::cout << it->first << " | " << it->second << std::endl;
											}*/
							/*debug for(std::vector<std::pair<Size,Size> >::iterator it = matches_shift_all.begin(); it != matches_shift_all.end(); ++it)
											{
												std::cout << it->first << " / " << it->second << std::endl;
											}*/
							/*debug for(std::list<std::pair<Size,Size> >::iterator it = matches_unshift.begin(); it != matches_unshift.end(); ++it)
											{
												std::cout << it->first << " * " << it->second << std::endl;
											}*/
							/*debug for(std::list<std::pair<Size,Size> >::iterator it = matches_shift.begin(); it != matches_shift.end(); ++it)
											{
												std::cout << it->first << " # " << it->second << std::endl;
											}*/

							//~ unite matches
							if(matches_unshift.size() > 0)
							{
								//~ both matches are > 0
								if(shift<pm_diff)
								{
									//~ prevent matches_unshift from being emptied in merge
									matches = matches_unshift;
									matches_shift.merge(matches,PairComparatorFirstElement< std::pair<Size,Size> >());
									matches = matches_shift;
								}
								else
								{
									matches = matches_unshift;
									matches.merge(matches_shift,PairComparatorFirstElement< std::pair<Size,Size> >());
								}
							}
							else
							{
								matches = matches_shift;
							}
							matches.unique(PairMatcherFirstElement< std::pair<Size,Size> >()); // dirty: takes only the match from the one where the other has been merged into - overlap HERE is highly unlikely
							//~ also this behaviour is partly wanted as a unshifted match is always?(maybe only in -shift) to be considered better than a shifted one

							//~ calc score
							for(std::list<std::pair<Size,Size> >::iterator it = matches.begin(); it != matches.end(); ++it)
							{
								score1 += s1[it->first].getIntensity();
								score2 += s2[it->second].getIntensity();
							}
							if((score1+score2 > best_score1 + best_score2) || (score1+score2 == best_score1+best_score2 && fabs(shift)<fabs(best_shift)))
							{
								best_score1=score1;
								best_score2=score2;
								best_shift=shift;
								best_matches = matches;
							}
						}
					}
				}
				if(best_shift==std::numeric_limits<double>::min())
				{
					best_shift=pm_diff;
				}
			}
		}
		///

	};
}
#endif //OPENMS_COMPARISON_SPECTRA_XCORRELATION_H
