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
			defaults_.setValue("max_shift", 150.0, "Defines the maximal absolut (in Da) shift between the two spectra");
			defaults_.setValue("min_dist", 57.0, "Defines the minimal distance (in Da) between the two peaks of one sort that may be connected in a sparse matching");
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

		void getXCorrelation(SpectrumType& s1, SpectrumType& s2, DoubleReal& best_score1 , DoubleReal& best_score2, DoubleReal& best_shift, std::list<std::pair<Size,Size> >& best_matches) const
		{
			DoubleReal peak_tolerance = (double)param_.getValue("peak_tolerance");
			DoubleReal parentmass_tolerance = (double)param_.getValue("parentmass_tolerance");
			DoubleReal max_shift = (double)param_.getValue("max_shift");
			DoubleReal shift_step = 2 * peak_tolerance;
			DoubleReal shift = 0;

			// reset the correlation values
			best_score1 = 0;
			best_score2 = 0;
			best_shift = 0;
			best_matches.clear();
			s1.sortByPosition();
			s2.sortByPosition();

			if(s2.getPrecursors().front().getMZ() < s1.getPrecursors().front().getMZ())
			{
				/// @improvement throw error? else maybe no range to cover because pm_diff is positive
			}
			DoubleReal pm_diff(s1.getPrecursors().front().getMZ() - s2.getPrecursors().front().getMZ());

			if(pm_diff>max_shift || s1.getPrecursors().front().getCharge()>2 || s2.getPrecursors().front().getCharge()>2)
			{
				//~ the above criteria disqualify s1 and s2 as (useful) spectral pairs from the start
				/// @improvement throw error?
				return;
			}


			shift = -parentmass_tolerance;
			DoubleReal range(parentmass_tolerance);
			for(; shift<=range; shift+=shift_step)
			{
				DoubleReal score1 = std::numeric_limits<double>::min();
				DoubleReal score2 = std::numeric_limits<double>::min();
				std::list<std::pair<Size,Size> > matches;

				/// @improvement add a correlation_scoring method matches that optimizes only the matchnumber (new max_sparse_overlap)
				if((String)param_.getValue("correlation_scoring")=="intensity")
				{
					//~ find matches in unshifted s1 to s2
					std::vector<std::pair<Size, Size> > matches_unshift_all;
					for(Size i = 0; i < s1.size(); ++i)
					{
						DoubleReal mz = s1[i].getMZ() + shift;
						ConstSpectrumIterator start(s2.MZBegin(mz-peak_tolerance));
						ConstSpectrumIterator end = (s2.MZEnd(mz+peak_tolerance));
						for(ConstSpectrumIterator it = start; it != end; ++it)
						{
							Size j = it - s2.begin();
							matches_unshift_all.push_back(std::pair<Size,Size>(i,j));
						}
					}
					//~ max sparses matches with DP
					std::list<std::pair<Size, Size> > matches_unshift;
					maxSparseMatches(s1,s2,matches_unshift_all,matches_unshift);


					if(abs(pm_diff)>parentmass_tolerance)
					{
						shift = pm_diff-parentmass_tolerance;
						range = pm_diff+parentmass_tolerance;
						//~ find matches in shifted s1 to s2
						std::vector<std::pair<Size, Size> > matches_shift_all;
						for(Size i = 0; i < s1.size(); ++i)
						{
							DoubleReal mz = s1[i].getMZ() + shift;
							ConstSpectrumIterator start(s2.MZBegin(mz-peak_tolerance));
							ConstSpectrumIterator end = (s2.MZEnd(mz+peak_tolerance));
							for(ConstSpectrumIterator it = start; it != end; ++it)
							{
								Size j = it - s2.begin();
								matches_shift_all.push_back(std::pair<Size,Size>(i,j));
							}
						}
						//~ max sparses matches with DP
						std::list<std::pair<Size, Size> > matches_shift;
						maxSparseMatches(s1,s2,matches_shift_all,matches_shift);
						//~ unite matches
						set_union(matches_unshift.begin(),matches_unshift.end(),matches_shift.begin(),matches_shift.end(),matches.begin());	//thanks to Size it will sort correctly
						matches.unique(matchPredicate); // dirty: takes only the match with higher index in s2 - but a overlap HERE is highly unlikely
					}
					else
					{
						matches = matches_unshift;
					}

					//~ calc score
					for(std::list<std::pair<Size,Size> >::iterator it = matches.begin(); it != matches.end(); ++it)
					{
						score1 += s1[it->first].getIntensity();
						score2 += s2[it->second].getIntensity();
					}
					if((score1+score2 > best_score1 + best_score2) || (score1+score2 == best_score1+best_score2 && abs(shift)<abs(best_shift)))
		      {
						best_score1=score1;
						best_score2=score2;
						best_shift=shift;
						best_matches = matches;
					}
				}
			}
		}

		void maxSparseMatches(const SpectrumType& s1, const SpectrumType& s2, std::vector<std::pair<Size,Size> >& all_matches, std::list<std::pair<Size,Size> >& best_matches) const
		{
			DoubleReal peak_tolerance = (double)param_.getValue("peak_tolerance");
			DoubleReal min_dist = (double)param_.getValue("min_dist");

			DoubleReal twocol [2] = { 0,0 };
			std::vector< DoubleReal[] > dp_table(all_matches.size()+1, twocol); // [predecessor, predecessor score], size+1 due to DP initialization in 0
			Size next_too_close;		// Index for the match (in dp_table and all_matches) that is next BEHIND i and too close
			Size max_index;					// Index for best match BEHIND next_too_close with the highest score (already cumulated) in the dp_table

			float best_match;     // Best path score so far
			int best_match_index;  // Index of the last in the best path so far

			max_index = 0;       next_too_close = 1;
			best_match = 0;    best_match_index = 0;
			for (Size i=1; i< dp_table.size(); ++i)
			{
				while(next_too_close < i and
							s1[all_matches[next_too_close-1].first].getMZ() <= (s1[all_matches[i-1].first].getMZ()-min_dist+2*peak_tolerance) and
							s2[all_matches[next_too_close-1].second].getMZ() <= (s2[all_matches[i-1].second].getMZ()-min_dist+2*peak_tolerance))
				{
					// This is executed only when next_too_close is advanced (then there may be a new max_index?)
					if (dp_table[next_too_close][1] > dp_table[max_index][1])
					{
						max_index = next_too_close;
					}
					next_too_close++;
				}
				dp_table[i][0] = max_index;
				dp_table[i][1] = dp_table[max_index][1] + s1[all_matches[i-1].first][1] + s2[all_matches[i-1].second][1]; //path score with intensities!

				// for traceback: best path ends in best_match_index
				if (dp_table[i][1]>best_match)
				{
					best_match=dp_table[i][1];
					best_match_index=i;
				}
			}

			std::list<Size> best_path;
			while (best_match_index>0)
			{
				best_path.push_back(best_match_index-1);
				best_match_index=(Size)dp_table[best_match_index][0];
			}

			Size best_path_size = best_path.size();
			//~ best_matches.reserve(best_path_size); not for lists

			for(Size i=0; i < best_path.size(); ++i)
			{
				best_matches.push_back(std::pair<Size,Size> (all_matches[best_path.back()].first,all_matches[best_path.back()].second));
				best_path.pop_back();
			}
		}

		//~ struct compareByFirst
  //~ : std::binary_function<const std::pair<Size,Size>,const std::pair<Size,Size>,bool>
//~ {
  //~ bool operator() (IntRealString left, IntRealString right ) const { return left.i_ < right.i_; }
//~ };
		// a binary predicate (as function) for the unique function of a list:
		bool matchPredicate (const std::pair<Size,Size>& first, const std::pair<Size,Size>& second) const
		{
			return ( first.first == second.first );
		}

	};
}
#endif //OPENMS_COMPARISON_SPECTRA_XCORRELATION_H
