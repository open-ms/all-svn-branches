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
#ifndef OPENMS_COMPARISON_SPECTRA_ANTISYMMETRICALIGNMENT_H
#define OPENMS_COMPARISON_SPECTRA_ANTISYMMETRICALIGNMENT_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

#include <vector>
#include <map>
#include <utility>
#include <algorithm>

#define ALIGNMENT_DEBUG
#undef  ALIGNMENT_DEBUG

namespace OpenMS
{

	//~ warning: unused parameter 'parent_mass2'
	//~ warning: unused parameter 'common_scores'
	//~ warning: unused parameter 'common_shifted_scores'

struct ASA_DP_Tables
{
	DoubleReal pm_s1, pm_s2;

	std::vector<int> common;
	std::vector<int> common_shifted;
	std::vector<double> common_scores;
	std::vector<double> common_shifted_scores;
	std::vector<int> prev;
	std::vector<int> next;
	std::vector<int> prev_shifted;
	std::vector<int> next_shifted;
	std::vector<std::vector<int> > left_jumps;
	std::vector<std::vector<int> > left_jumps_shifted;
	std::vector<std::vector<int> > right_jumps;
	std::vector<std::vector<int> > right_jumps_shifted;
	std::vector<std::vector<int> > left_neighbors;
	std::vector<std::vector<int> > right_neighbors;

	std::vector<DoubleReal> peaks;
	std::vector<DoubleReal> peaks2;
	float delta;
	float delta2;

	std::vector<double> prefix1;
	std::vector<std::vector<double> > suffix1;
	std::vector<std::vector<double> > prefix2;
	std::vector<double> suffix2;
	std::vector<double> prefix1_L;
	std::vector<double> suffix1max;
	std::vector<double> suffix1_R;
	std::vector<double> prefix2max;
	std::vector<double> prefix2_L;
	std::vector<double> suffix2_R;
	std::vector<std::vector<double> > D1;
	std::vector<std::vector<std::vector<double> > > D2;
	std::vector<std::vector<std::vector<double> > > D3;
	std::vector<std::vector<double> > D2max;
	std::vector<std::vector<double> > D3max;
	std::vector<std::vector<double> > M1_L;
	std::vector<std::vector<double> > M1_R;
	std::vector<std::vector<double> > M2_L;
	std::vector<std::vector<std::vector<double> > > M2_R;
	std::vector<std::vector<std::vector<double> > > M3_L;
	std::vector<std::vector<double> > M3_R;

	enum celltype {cell_prefix1, cell_suffix1, cell_prefix2, cell_suffix2,
				cell_prefix1_L, cell_suffix1_R, cell_prefix2_L, cell_suffix2_R,
				cell_D1, cell_D2, cell_D3,
				cell_M1_L, cell_M1_R, cell_M2_L, cell_M2_R,
				cell_M3_L, cell_M3_R, cell_Invalid};

	double best_score;
	int best_i;
	int best_j;
	int best_s;
	celltype best_t;

	ASA_DP_Tables(Size n, Size m, Size l, DoubleReal pm1, DoubleReal pm2) // n = res_1.size();  m = common.size()/2; l = res_2.size();
				:pm_s1(pm1), pm_s2(pm2),best_score(0.0),best_i(-1),best_j(-1),best_s(-1),
				peaks(n),peaks2(l),
				common(n,-1), common_shifted(n,-1), common_scores(n,0), common_shifted_scores(n,0),
				prev(n, -1), next(n,-1), prev_shifted(n,-1), next_shifted(n,-1),
				left_jumps(n), left_jumps_shifted(n), right_jumps(n),right_jumps_shifted(n), left_neighbors(n), right_neighbors(n),
				prefix1(m, std::numeric_limits<DoubleReal>::min()), suffix1(m),prefix2(m), suffix2(m, std::numeric_limits<DoubleReal>::min()), prefix1_L(m, std::numeric_limits<DoubleReal>::min()), suffix1max(m, std::numeric_limits<DoubleReal>::min()), suffix1_R(m, std::numeric_limits<DoubleReal>::min()), prefix2max(m, std::numeric_limits<DoubleReal>::min()), prefix2_L(m, std::numeric_limits<DoubleReal>::min()), suffix2_R(m, std::numeric_limits<DoubleReal>::min()),
				D1(m,std::vector<double>(m, std::numeric_limits<DoubleReal>::min())),
				D2(m,std::vector<std::vector<double> >(m)),
				D3(m,std::vector<std::vector<double> >(m)),
				D2max(m,std::vector<double>(m, std::numeric_limits<DoubleReal>::min())),
				D3max(m,std::vector<double>(m, std::numeric_limits<DoubleReal>::min())),
				M1_L(m,std::vector<double>(m, std::numeric_limits<DoubleReal>::min())),
				M1_R(m,std::vector<double>(m, std::numeric_limits<DoubleReal>::min())),
				M2_L(m,std::vector<double>(m, std::numeric_limits<DoubleReal>::min())),
				M2_R(m,std::vector<std::vector<double> >(m)),
				M3_L(m,std::vector<std::vector<double> >(m)),
				M3_R(m,std::vector<double>(m, std::numeric_limits<DoubleReal>::min())),best_t(cell_Invalid)
	{
		//~ celltype best_t = cell_Invalid;
	}
};
///

	/**
		@brief Aligns the peaks of two spectra

		@htmlinclude OpenMS_AntisymmetricAlignment.parameters

		@ingroup SpectraComparison
	*/
	template <typename PeakT = Peak1D>
	class AntisymmetricAlignment
		: public DefaultParamHandler
	{
		typedef MSSpectrum<PeakT> SpectrumType;

	public:

		// @name Constructors and Destructors
		// @{
		/// default constructor
		AntisymmetricAlignment()
		: DefaultParamHandler("AntisymmetricAlignment")
		{
				defaults_.setValue("peak_tolerance", 0.3, "Defines the absolut (in Da) peak tolerance");
				defaults_.setValue("parentmass_tolerance", 3.0, "Defines the absolut (in Da) parent mass tolerance");
				defaults_.setValue("min_dist", 57.0214637230 , "Defines the minimal distance (in Da) between the two peaks of one sort that may be connected in a sparse matching");
				defaults_.setValue("sv_penalty", -5.0 , "Defines the penalty for being in the same 'vertex' in DP");
				defaults_.setValue("dif_penalty", -5.0 , "Defines the penalty for jumping the parentmass difference in DP");
				defaultsToParam_();
		}

		/// copy constructor
		AntisymmetricAlignment(const AntisymmetricAlignment& source)
			: DefaultParamHandler(source)
		{
		}

		/// destructor
		virtual ~AntisymmetricAlignment()
		{
		}

		/// assignment operator
		AntisymmetricAlignment& operator = (const AntisymmetricAlignment& source)
		{
			if (this != &source)
			{
				DefaultParamHandler::operator = (source);
			}
			return *this;
		}

		// @}

		/**
			@brief Will find valid (aa difference separated) consecutive peaks to the left and right of a given peak in a given spectrum

			@param index to given peak
			@param s1sym the given spectrum
			@param jump_masses the valid mass jumps for the matches
			@param left_matches will contain the valid peak indices to the left of index
			@param right_matches will contain the valid peak indices to the right of index
		*/
		void getAASteps(Size& index, SpectrumType& s1sym, std::vector<Real>& jump_masses,
						std::vector<int>& leftMatches, std::vector<int>& rightMatches) const
		{
			Real peak_tolerance = (Real)param_.getValue("peak_tolerance");
			Real pt_plus = peak_tolerance +.00001f;
			Size match_counter = 0;
			DoubleReal current_diff;
			DoubleReal max_jump = jump_masses[jump_masses.size()-1]+peak_tolerance;

			//~ if we start leftmost this can be skipped as there is nothing to jump left to
			if(index>0)
			{
				//~ max |index| left jumps possible
				leftMatches.resize(index);
				//~ but so far zero
				match_counter=0;

				//~ typename SpectrumType::Iterator it_lo1 = s1sym.begin()+index-1;
				//~ typename SpectrumType::Iterator it_lo2 = s1sym.MZBegin(s1sym[index].getMZ()-max_jump-.00001);

				std::reverse_iterator< typename SpectrumType::Iterator > it_lo1(s1sym.begin()+index); //points to before 'index'
				std::reverse_iterator< typename SpectrumType::Iterator > it_lo2(s1sym.MZBegin(s1sym[index].getMZ()-max_jump-.00001)); //points to before MZBegin... maybe before s1sym.begin() -std::min(s1sym[0].getMZ(), MZBegin...)?
				std::vector<Real>::iterator jm_left_pivot = jump_masses.begin();

				/// @improvement posssible speedup or slowdown by comp if it_lo2.getMZ() <= it_lo1.getMZ() and only go on then ...
				//~ iterate from lo1 to lo2(not included!) and find if mz diff between lo1(iterated) and index concurs with a entry in jump_masses +-peak_tolerance
				for(;it_lo1 != it_lo2; ++it_lo1)
				{
					current_diff = s1sym[index].getMZ()-it_lo1->getMZ();
					//~ difference will only grow -> advance the pivot
					jm_left_pivot = std::lower_bound(jm_left_pivot, jump_masses.end(), (Real)current_diff-pt_plus);

					if(jm_left_pivot!=jump_masses.end() and EqualInTolerance<Real>(pt_plus)((Real)current_diff,*jm_left_pivot))
					{
						leftMatches[match_counter++]=(s1sym.size()-1)-(it_lo1-s1sym.rbegin());
					}
				}
				leftMatches.resize(match_counter);
			}

			if(index<(s1sym.size()-1))
			{
				//~ max |s1sym|-|index| right jumps possible
				rightMatches.resize(s1sym.size()-index);
				//~ but so far zero
				match_counter=0;

				typename SpectrumType::Iterator it_hi1 = s1sym.begin()+index+1;
				typename SpectrumType::Iterator it_hi2 = s1sym.MZEnd(s1sym[index].getMZ()+max_jump+.00001);
				std::vector<Real>::iterator jm_right_pivot = jump_masses.begin();

				//~ iterate from hi1 to <hi2 and find if mz diff between it_hi1(iterated) and index concurs with a entry in jump_masses +-peak_tolerance
				for(;it_hi1 != it_hi2; ++it_hi1)
				{
					current_diff = it_hi1->getMZ()-s1sym[index].getMZ();
					//~ difference will only grow -> advance the pivot
					jm_right_pivot = std::lower_bound(jm_right_pivot, jump_masses.end(), (Real)current_diff-pt_plus);

					if(jm_right_pivot!=jump_masses.end() and EqualInTolerance<Real>(pt_plus)((Real)current_diff,*jm_right_pivot))
					{
						rightMatches[match_counter++]=it_hi1-s1sym.begin();
					}
				}
				rightMatches.resize(match_counter);
			}
		}
		///

		/**
			@brief Will yield a symmetric spectrum

			@param spectrum the peak spectrum to be made symmetric

			@return a symmetric version of spectrum with a integer MetaDataArray, where x[i]>0 indicates a synthetic peak.

			result is a spectrum with property sym[i].getMZ + sym[n-i-1].getMZ = parentMass. Synthetic peaks will be added with a inntensity of zero and in addition remarked in a metadataarray.
		@important previous sortByPosition() of spectrum
		@important is the 0 and pm peak addition neccessary!?
		*/
		MSSpectrum<PeakT> getSymetricSpectrum(const MSSpectrum<PeakT>& spectrum) const
		{
			SpectrumType sym(spectrum);
			if(sym.empty())
			{
				return sym;
			}

			typename SpectrumType::IntegerDataArray ida_symmetric; ida_symmetric.resize(sym.size(),0);
			sym.getIntegerDataArrays().push_back(ida_symmetric);
			sym.getIntegerDataArrays().push_back(ida_symmetric);

			Real peak_tolerance = (Real)param_.getValue("peak_tolerance");

			/// @important uncharged mass difference!
			DoubleReal pm = sym.getPrecursors().front().getUnchargedMass();

			Size sym_pasttheend_index = sym.size();
			for(Size sym_index = 0; sym_index < sym_pasttheend_index; ++sym_index)
			{
				if(sym.getIntegerDataArrays().front()[sym_index]!=0)
				{
					continue;
				}
				DoubleReal target = pm-sym[sym_index].getMZ()+2*Constants::PROTON_MASS_U;

				typename SpectrumType::Iterator lo = sym.MZBegin(target-peak_tolerance);
				typename SpectrumType::Iterator hi = sym.MZEnd(lo,target+peak_tolerance,sym.begin()+sym_pasttheend_index);

				//~ pairs get metadataarray entries -n/n where n is number of the pair (count from 1) and the lower mz peak has the negative pair number
				int sign_it_mz = -1;
				//~ int sign_target = sign_it_mz * -1
				((target-sym[sym_index].getMZ())<0)?sign_it_mz*=(-1):sign_it_mz=sign_it_mz;

				//~ empty range - add a peak
				if(lo==hi)
				{
					/// @improvement catch target-is-negative case
					PeakT tmp;
					tmp.setPosition(target);
					sym.push_back(tmp);
					sym.getIntegerDataArrays().front().push_back( sign_it_mz*(-1)* sym_index+1 );
					sym.getIntegerDataArrays().front()[sym_index]= sign_it_mz * sym_index+1;
					sym.getIntegerDataArrays()[1].push_back(1);
				}
				else
				{
					//~ get best match in tolerance window
					for(typename SpectrumType::Iterator it_mirror_mz = lo; it_mirror_mz != hi; ++it_mirror_mz)
					{
						if( fabs(it_mirror_mz->getMZ()-sym[sym_index].getMZ()) < fabs(lo->getMZ()-sym[sym_index].getMZ()) )
						{
							lo=it_mirror_mz;
						}
					}
					sym.getIntegerDataArrays().front()[sym_index]= sign_it_mz*(sym_index+1);
					sym.getIntegerDataArrays().front()[(lo-sym.begin())]= sign_it_mz*(-1)*(sym_index+1);
				}
			}

			std::swap(sym.getIntegerDataArrays()[0],sym.getIntegerDataArrays()[1]);
			sym.getIntegerDataArrays()[0].setName("synthetic peaks");
			sym.getIntegerDataArrays().resize(1);
			sym.sortByPosition();

			return sym;
		}
		///

		/**
			@brief Will prune the spectrum s into a version without non-singly charged peaks

			@param spectrum a spectrum to be deconvoluted

			@return a version of spectrum with charge deconvoluted peaks to charge one

		@important works only if id is present
		*/
		void filterHighCharges(MSSpectrum<PeakT>& spectrum,Size c_lo,Size c_hi) const
		{
			Real peak_tolerance = param_.getValue("peak_tolerance");

			if(spectrum.empty() or spectrum.getPeptideIdentifications().empty() or spectrum.getPeptideIdentifications().front().getHits().empty())
			{
				return;
			}
			///@improvement without id like denovo over isotopes?

			TheoreticalSpectrumGenerator tsg;
			RichPeakSpectrum ts;
			tsg.getSpectrum(ts, spectrum.getPeptideIdentifications().front().getHits().front().getSequence(), c_lo, c_hi);
			ts.sortByPosition();

			typename MSSpectrum<PeakT>::Iterator lo = spectrum.begin();
			typename MSSpectrum<PeakT>::Iterator hi = spectrum.end();
			for(Size i = 0; i < ts.size(); ++i)
			{
				DoubleReal target(ts[i].getMZ());
				lo = spectrum.MZBegin(lo,target-peak_tolerance-0.00001,spectrum.end());
				hi = spectrum.MZEnd(lo,target+peak_tolerance+0.00001,spectrum.end());
				for(; lo!=hi; ++lo)
				{
					lo->setIntensity(0);
				}
			}

			lo = spectrum.begin();
			while(lo!=spectrum.end())
			{
				if(lo->getIntensity() == 0)
				{
					lo = spectrum.erase(lo);
					///@improvement is singly charged version yet contained? else this could be added synthetically
				}
				else
				{
					++lo;
				}
			}

			return;
		}
		///

		/**
			@brief Will prune the first spectrum with respect to matching/shift-matching peaks of the second spectrum

			@param spectrum to be pruned
			@param spectrum to be matched

			metadata will contain matching information!
		*/
		void matchPrunePeaks(MSSpectrum<PeakT>& spec_1,MSSpectrum<PeakT>& spec_2, DoubleReal pm_diff) const
		{
			Real peak_tolerance = param_.getValue("peak_tolerance");

			typename SpectrumType::IntegerDataArray ida_matches; ida_matches.resize(spec_1.size(),-1);
			spec_1.getIntegerDataArrays().insert(spec_1.getIntegerDataArrays().begin()+1, ida_matches);
			spec_1.getIntegerDataArrays()[1].setName("normal match indices");
			//~ indices in spec_2 for pm_diff shift matches
			spec_1.getIntegerDataArrays().insert(spec_1.getIntegerDataArrays().begin()+2, ida_matches);
			spec_1.getIntegerDataArrays()[2].setName("shift match indices");

			//~ remove peaks without a match or complement peak match: zero intensity in spec_1 for unmatching peaks (if both of a pair dont fit to any in spec_2)
			for(typename SpectrumType::Iterator it_res = spec_1.begin(); it_res != spec_1.end(); ++it_res)
			{
				bool no_match_normal = false;
				typename SpectrumType::ConstIterator lo = spec_2.MZBegin(it_res->getMZ()-peak_tolerance);
				typename SpectrumType::ConstIterator hi = spec_2.MZEnd(it_res->getMZ()+peak_tolerance);
				if(lo==hi)
				{
					no_match_normal = true;
				}
				else
				{
					//~ get best match in tolerance window
					for(typename SpectrumType::ConstIterator it_match_mz = lo; it_match_mz != hi; ++it_match_mz)
					{
						if( fabs(it_match_mz->getMZ()-it_res->getMZ()) < fabs(lo->getMZ()-it_res->getMZ()) )
						{
							lo=it_match_mz;
						}
					}
					spec_1.getIntegerDataArrays()[1][(it_res-spec_1.begin())] = lo-spec_2.begin();
				}

				lo = spec_2.MZBegin(it_res->getMZ()+pm_diff-peak_tolerance);
				hi = spec_2.MZEnd(it_res->getMZ()+pm_diff+peak_tolerance);
				if(lo==hi)
				{
					if(no_match_normal)
					{
						it_res->setIntensity(0);
					}
				}
				else
				{
					//~ get best match in tolerance window
					for(typename SpectrumType::ConstIterator it_match_mz = lo; it_match_mz != hi; ++it_match_mz)
					{
						if( fabs(it_match_mz->getMZ()-it_res->getMZ()) < fabs(lo->getMZ()-it_res->getMZ()) )
						{
							lo=it_match_mz;
						}
					}
					spec_1.getIntegerDataArrays()[2][(it_res-spec_1.begin())] = lo-spec_2.begin();
				}
			}

			//~ number of complementing peaks in spec_1 that both do not match the right peaks in spec_2
			Size removed_pairs = 0;
			//~ remove the peak pairs mentioned above
			Size spec_1_index = 0;
			while(spec_1_index < spec_1.size() and (Size)spec_1_index < (Size)((spec_1.size()/2) - removed_pairs) )
			{
				if(spec_1[spec_1_index].getIntensity() == 0)
				{
					Size res_mirror_index = spec_1.size() - 1 - spec_1_index;
					if(spec_1[res_mirror_index].getIntensity() == 0 )
					{
						for(Size i = 0; i < 3; ++i)
						{
							typename SpectrumType::IntegerDataArray::iterator it_mirror_ida = spec_1.getIntegerDataArrays()[i].begin()+res_mirror_index;
							spec_1.getIntegerDataArrays()[i].erase(it_mirror_ida);
						}
						for(Size i = 0; i < 3; ++i)
						{
							typename SpectrumType::IntegerDataArray::iterator it_ida = spec_1.getIntegerDataArrays()[i].begin()+spec_1_index;
							spec_1.getIntegerDataArrays()[i].erase(it_ida);
						}
						typename SpectrumType::Iterator res_mirror_it = spec_1.begin() + res_mirror_index;
						spec_1.erase(res_mirror_it);
						typename SpectrumType::Iterator spec_1_it = spec_1.begin() + spec_1_index;
						spec_1.erase(spec_1_it);
						++removed_pairs;
					}
					else
					{
						++spec_1_index;
					}
				}
				else
				{
					++spec_1_index;
				}
			}
		}
		///


		/**
			@brief will traceback the antisymmetric path of a found alignment

			@param dp_tables the filled dynamic programming tables
			@param dif_penalty the penalty score for using a peak in the alignment that might origin from the other ion species
			@param score of the alignment

			@return a pair of vectors of indices containing the path of aligned indices in a row

			...
		*/
		std::pair< std::vector<int>,std::vector<int> > traceback(ASA_DP_Tables& dp_tables, DoubleReal& dif_penalty, DoubleReal& score) const
		{
			int n0 = dp_tables.prev.size();
			int N = n0-1;

			std::vector<int> path1;
			std::vector<int> path2;

			if (dp_tables.best_i == -1)
			{
				score = dp_tables.best_score/dp_tables.pm_s1;
				return std::pair<std::vector<int>,std::vector<int> >(path1, path2);
			}

			int i = dp_tables.best_i;
			int j = dp_tables.best_j;
			int s = dp_tables.best_s;
			ASA_DP_Tables::celltype t;
			t = dp_tables.best_t;
			while(1)
			{
				if(t == ASA_DP_Tables::cell_prefix1)
				{
					path1.insert(path1.begin(), i);
					double prev_score = 0.0;
					int index = 0;
					int next_i = 0;
					for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
					{
						int i2 = *it;
						if(dp_tables.prefix1[i2] > prev_score)
						{
							prev_score = dp_tables.prefix1[i2];
							index = 1;
							next_i = i2;
						}
					}
					int i2 = dp_tables.prev[i];
					if(i2 != -1 && dp_tables.prefix1_L[i2] > prev_score)
					{
						prev_score = dp_tables.prefix1_L[i2];
						index = 2;
						next_i = i2;
					}

					i = next_i;
					if(index == 0)
					{
						break;
					}
					else if (index == 1)
					{
				// t is unchanged
					}
					else
					{
						t = ASA_DP_Tables::cell_prefix1_L;
					}
				}
				else if(t == ASA_DP_Tables::cell_suffix1)
				{
					path1.push_back(N-j);
					if (s == 0)
					{
						double next_score = 0.0;
						int index = 0;
						int next_j = 0;
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
						{
							int j2 = *it;
							if(dp_tables.suffix2[N-j2]+dif_penalty > next_score)
							{
								next_score = dp_tables.suffix2[N-j2]+dif_penalty;
								index = 1;
								next_j = N-j2;
							}
						}
						int j2 = dp_tables.next[N-j];
						if(j2 != -1 && dp_tables.suffix2_R[N-j2]+dif_penalty > next_score)
						{
							next_score = dp_tables.suffix2_R[N-j2]+dif_penalty;
							index = 2;
							next_j = N-j2;
						}
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps_shifted[N-j].begin(); it != dp_tables.right_jumps_shifted[N-j].end(); ++it)
						{
							int j2 = *it;
							if(dp_tables.suffix1max[N-j2] > next_score)
							{
								next_score = dp_tables.suffix1max[N-j2];
								index = 3;
								next_j = N-j2;
							}
						}
						j2 = dp_tables.next_shifted[N-j];
						if(j2 != -1 && dp_tables.suffix1_R[N-j2] > next_score)
						{
							next_score = dp_tables.suffix1_R[N-j2];
							index = 4;
							next_j = N-j2;
						}

						j = next_j;
						if(index == 0)
						{
							break;
						}
						else if(index == 1)
						{
							t = ASA_DP_Tables::cell_suffix2;
						}
						else if(index == 2)
						{
							t = ASA_DP_Tables::cell_suffix2_R;
						}
						else if(index == 3)
						{
							// t is unchanged
							s = (std::max_element(dp_tables.suffix1[j].begin(),dp_tables.suffix1[j].end())-dp_tables.suffix1[j].begin());
						}
						else
						{
							t = ASA_DP_Tables::cell_suffix1_R;
						}
					}
					else
					{
						// t is unchanged
						j = dp_tables.right_neighbors[N-j][s-1];
						j = N-j;
						s = (std::max_element(dp_tables.suffix1[j].begin(),dp_tables.suffix1[j].end())-dp_tables.suffix1[j].begin());
					}
				}
				else if(t == ASA_DP_Tables::cell_prefix2)
				{
					path2.insert(path2.begin(), i);
					if(s == 0)
					{
						double prev_score = 0.0;
						int index = 0;
						int next_i = 0;
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
						{
							int i2 = *it;
							if(dp_tables.prefix1[i2]+dif_penalty > prev_score)
							{
								prev_score = dp_tables.prefix1[i2]+dif_penalty;
								index = 1;
								next_i = i2;
							}
						}
						int i2 = dp_tables.prev[i];
						if(i2 != -1 && dp_tables.prefix1_L[i2]+dif_penalty > prev_score)
						{
							prev_score = dp_tables.prefix1_L[i2]+dif_penalty;
							index = 2;
							next_i = i2;
						}
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps_shifted[i].begin(); it != dp_tables.left_jumps_shifted[i].end(); ++it)
						{
							int i2 = *it;
							if(dp_tables.prefix2max[i2] > prev_score)
							{
								prev_score = dp_tables.prefix2max[i2];
								index = 3;
								next_i = i2;
							}
						}
						i2 = dp_tables.prev_shifted[i];
						if(i2 != -1 && dp_tables.prefix2_L[i2] > prev_score)
						{
							prev_score = dp_tables.prefix2_L[i2];
							index = 4;
							next_i = i2;
						}

						i = next_i;
						if(index == 0)
						{
							break;
						}
						else if(index == 1)
						{
							t = ASA_DP_Tables::cell_prefix1;
						}
						else if(index == 2)
						{
							t = ASA_DP_Tables::cell_prefix1_L;
						}
						else if(index == 3)
						{
							// t is unchanged
							s = (std::max_element(dp_tables.prefix2[i].begin(),dp_tables.prefix2[i].end())-dp_tables.prefix2[i].begin());
						}
						else
						{
							t = ASA_DP_Tables::cell_prefix2_L;
						}
					}
					else
					{
						// t is unchanged
						i = dp_tables.left_neighbors[i][s-1];
						s = (std::max_element(dp_tables.prefix2[i].begin(),dp_tables.prefix2[i].end())-dp_tables.prefix2[i].begin());
					}
				}
				else if(t == ASA_DP_Tables::cell_suffix2)
				{
					path2.push_back(N-j);
					double next_score = 0.0;
					int index = 0;
					int next_j = 0;
					for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
					{
						int j2 = *it;
						if (dp_tables.suffix2[N-j2] > next_score)
						{
							next_score = dp_tables.suffix2[N-j2];
							index = 1;
							next_j = N-j2;
						}
					}
					int j2 = dp_tables.next[N-j];
					if (j2 != -1 && dp_tables.suffix2_R[N-j2] > next_score) {
					next_score = dp_tables.suffix2_R[N-j2];
					index = 2;
					next_j = N-j2;
				}

				j = next_j;
				if(index == 0)
				{
					break;
				}
				else if(index == 1)
				{
					// t is unchanged
				}
				else
				{
					t = ASA_DP_Tables::cell_suffix2_R;
				}
			}
				else if(t == ASA_DP_Tables::cell_prefix1_L)
				{
					if(dp_tables.prefix1_L[i] == dp_tables.prefix1[i])
					{
						t = ASA_DP_Tables::cell_prefix1;
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_prefix2_L)
				{
					if(dp_tables.prefix2_L[i] == dp_tables.prefix2max[i])
					{
						t = ASA_DP_Tables::cell_prefix2;
						s = (std::max_element(dp_tables.prefix2[i].begin(),dp_tables.prefix2[i].end())-dp_tables.prefix2[i].begin());
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_suffix1_R)
				{
					if (dp_tables.suffix1_R[j] == dp_tables.suffix1max[j])
					{
						t = ASA_DP_Tables::cell_suffix1;
						s = (std::max_element(dp_tables.suffix1[j].begin(),dp_tables.suffix1[j].end())-dp_tables.suffix1[j].begin());
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_suffix2_R)
				{
					if(dp_tables.suffix2_R[j] == dp_tables.suffix2[j])
					{
						t = ASA_DP_Tables::cell_suffix2;
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_D1)
				{
					if(i >= j)
					{
						path1.insert(path1.begin(), i);
						double prev_score = dp_tables.suffix2[j]+dif_penalty;
						int index = 0;
						int next_i = 0;
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
						{
							int i2 = *it;
							if(dp_tables.D1[i2][j] > prev_score)
							{
								prev_score = dp_tables.D1[i2][j];
								index = 1;
								next_i = i2;
							}
						}
						int i2 = dp_tables.prev[i];
						if(i2 != -1 && dp_tables.M1_L[i2][j] > prev_score)
						{
							prev_score = dp_tables.M1_L[i2][j];
							index = 2;
							next_i = i2;
						}

						i = next_i;
						// note that if index=0, then the value of i is irrelevant.
						if(index == 0)
						{
							t = ASA_DP_Tables::cell_suffix2;
						}
						else if(index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = ASA_DP_Tables::cell_M1_L;
						}
					}
					else
					{
						path2.push_back(N-j);
						double next_score = dp_tables.prefix1[i]+dif_penalty;
						int index = 0;
						int next_j = 0;
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
						{
							int j2 = *it;
							if (dp_tables.D1[i][N-j2] > next_score)
							{
								next_score = dp_tables.D1[i][N-j2];
								index = 1;
								next_j = N-j2;
							}
						}
						int j2 = dp_tables.next[N-j];
						if(j2 != -1 && dp_tables.M1_R[i][N-j2] > next_score)
						{
							next_score = dp_tables.M1_R[i][N-j2];
							index = 2;
							next_j = N-j2;
						}

						j = next_j;
						if(index == 0)
						{
							t = ASA_DP_Tables::cell_prefix1;
						}
						else if(index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = ASA_DP_Tables::cell_M1_R;
						}
					}
				}
				else if(t == ASA_DP_Tables::cell_D2)
				{
					if (i > j)
					{
						path2.insert(path2.begin(), i);
						if(s == 0)
						{
							double prev_score = dp_tables.suffix2[j];
							int index = 0;
							int next_i = 0;
							for (std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
							{
								int i2 = *it;
								if(dp_tables.D1[i2][j] > prev_score)
								{
									prev_score = dp_tables.D1[i2][j];
									index = 1;
									next_i = i2;
								}
							}
							int i2 = dp_tables.prev[i];
							if(i2 != -1 && dp_tables.M1_L[i2][j] > prev_score)
							{
								prev_score = dp_tables.M1_L[i2][j];
								index = 2;
								next_i = i2;
							}
							for(std::vector<int>::const_iterator it = dp_tables.left_jumps_shifted[i].begin(); it != dp_tables.left_jumps_shifted[i].end(); ++it)
							{
								int i2 = *it;
								if(dp_tables.D2max[i2][j] > prev_score)
								{
									prev_score = dp_tables.D2max[i2][j];
									index = 3;
									next_i = i2;
								}
							}
							i2 = dp_tables.prev_shifted[i];
							if(i2 != -1 && dp_tables.M2_L[i2][j] > prev_score)
							{
								prev_score = dp_tables.M2_L[i2][j];
								index = 4;
								next_i = i2;
							}

							i = next_i;
							if(index == 0)
							{
								t = ASA_DP_Tables::cell_suffix2;
							}
							else if(index == 1)
							{
								t = ASA_DP_Tables::cell_D1;
							}
							else if(index == 2)
							{
								t = ASA_DP_Tables::cell_M1_L;
							}
							else if(index == 3)
							{
								// t is unchanged
								s = (std::max_element(dp_tables.D2[i][j].begin(),dp_tables.D2[i][j].end())-dp_tables.D2[i][j].begin());
							}
							else
							{
								t = ASA_DP_Tables::cell_M2_L;
							}
						}
						else
						{
							// t is unchanged
							i = dp_tables.left_neighbors[i][s-1];
							s = (std::max_element(dp_tables.D2[i][j].begin(),dp_tables.D2[i][j].end())-dp_tables.D2[i][j].begin());
						}
					}
					else
					{ // i <= j
						path2.push_back(N-j);
						double next_score = dp_tables.prefix2[i][s];
						int index = 0;
						int next_j = 0;
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
						{
							int j2 = *it;
							if(dp_tables.D2[i][N-j2][s] > next_score)
							{
								next_score = dp_tables.D2[i][N-j2][s];
								index = 1;
								next_j = N-j2;
							}
						}
						int j2 = dp_tables.next[N-j];
						if(j2 != -1 && dp_tables.M2_R[i][N-j2][s] > next_score)
						{
							index = 2;
							next_j = N-j2;
						}

						j = next_j;
						if(index == 0)
						{
							t = ASA_DP_Tables::cell_prefix2;
						}
						else if(index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = ASA_DP_Tables::cell_M2_R;
						}
					}
				}
				else if(t == ASA_DP_Tables::cell_D3)
				{
					if(i >= j)
					{
						path1.insert(path1.begin(), i);
						double prev_score = dp_tables.suffix1[j][s];
						int index = 0;
						int next_i = 0;
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
						{
							int i2 = *it;
							if(dp_tables.D3[i2][j][s] > prev_score)
							{
								prev_score = dp_tables.D3[i2][j][s];
								index = 1;
								next_i = i2;
							}
						}
						int i2 = dp_tables.prev[i];
						if(i2 != -1 && dp_tables.M3_L[i2][j][s] > prev_score)
						{
							index = 2;
							next_i = i2;
						}

						i = next_i;
						if(index == 0)
						{
							t = ASA_DP_Tables::cell_suffix1;
						}
						else if (index == 1)
						{
							// t is unchanged
						}
						else
						{
							t = ASA_DP_Tables::cell_M3_L;
						}
					}
					else
					{ // i < j
						path1.push_back(N-j);
						if(s == 0)
						{
							double next_score = dp_tables.prefix1[i];
							int index = 0;
							int next_j = 0;
							for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
							{
								int j2 = *it;
								if(dp_tables.D1[i][N-j2] > next_score)
								{
									next_score = dp_tables.D1[i][N-j2];
									index = 1;
									next_j = N-j2;
								}
							}
							int j2 = dp_tables.next[N-j];
							if(j2 != -1 && dp_tables.M1_R[i][N-j2] > next_score)
							{
								next_score = dp_tables.M1_R[i][N-j2];
								index = 2;
								next_j = N-j2;
							}
							for(std::vector<int>::const_iterator it = dp_tables.right_jumps_shifted[N-j].begin(); it != dp_tables.right_jumps_shifted[N-j].end(); ++it)
							{
								int j2 = *it;
								if(dp_tables.D3max[i][N-j2] > next_score)
								{
									next_score = dp_tables.D3max[i][N-j2];
									index = 3;
									next_j = N-j2;
								}
							}
							j2 = dp_tables.next_shifted[N-j];
							if(j2 != -1 && dp_tables.M3_R[i][N-j2] > next_score)
							{
								next_score = dp_tables.M3_R[i][N-j2];
								index = 4;
								next_j = N-j2;
							}

							j = next_j;
							if(index == 0)
							{
								t = ASA_DP_Tables::cell_prefix1;
							}
							else if(index == 1)
							{
								t = ASA_DP_Tables::cell_D1;
							}
							else if(index == 2)
							{
								t = ASA_DP_Tables::cell_M1_R;
							}
							else if( index == 3)
							{
								// t is unchanged
								s = (std::max_element(dp_tables.D3[i][j].begin(),dp_tables.D3[i][j].end())-dp_tables.D3[i][j].begin());
							}
							else
							{
								t = ASA_DP_Tables::cell_M3_R;
							}
						}
						else
						{
							// t is unchanged
							j = dp_tables.right_neighbors[N-j][s-1];
							j = N-j;
							s = (std::max_element(dp_tables.D3[i][j].begin(),dp_tables.D3[i][j].end())-dp_tables.D3[i][j].begin());
						}
					}
				}
				else if(t == ASA_DP_Tables::cell_M1_L)
				{
					if(dp_tables.M1_L[i][j] == dp_tables.D1[i][j])
					{
						t = ASA_DP_Tables::cell_D1;
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_M2_L)
				{
					if (dp_tables.M2_L[i][j] == dp_tables.D2max[i][j])
					{
						t = ASA_DP_Tables::cell_D2;
						s = (std::max_element(dp_tables.D2[i][j].begin(),dp_tables.D2[i][j].end())-dp_tables.D2[i][j].begin());
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_M3_L)
				{
					if (dp_tables.M3_L[i][j][s] == dp_tables.D3[i][j][s])
					{
						t = ASA_DP_Tables::cell_D3;
					}
					else
					{
						i -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_M1_R)
				{
					if(dp_tables.M1_R[i][j] == dp_tables.D1[i][j])
					{
						t = ASA_DP_Tables::cell_D1;
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_M2_R)
				{
					if(dp_tables.M2_R[i][j][s] == dp_tables.D2[i][j][s])
					{
						t = ASA_DP_Tables::cell_D2;
					}
					else
					{
						j -= 1;
					}
				}
				else if(t == ASA_DP_Tables::cell_M3_R)
				{
					if(dp_tables.M3_R[i][j] == dp_tables.D3max[i][j])
					{
						t = ASA_DP_Tables::cell_D3;
						s = (std::max_element(dp_tables.D3[i][j].begin(),dp_tables.D3[i][j].end())-dp_tables.D3[i][j].begin());
					}
					else
					{
						j -= 1;
					}
				}
			}

			score = dp_tables.best_score/dp_tables.pm_s1;
			return std::pair<std::vector<int>,std::vector<int> >(path1, path2);
		}
		///

		/**
			@brief ...

			@param

			...
		*/
		void calculateAlignementPath(ASA_DP_Tables& dp_tables, Real peak_tolerance, DoubleReal sv_penalty, DoubleReal dif_penalty) const
		{
			int n0 = dp_tables.common.size();
			int N = n0-1;
			int n = (n0+1)/2;

			// fill prefix1
			for(int i = 0; i < n; ++i)
			{
				if(dp_tables.common[i] != -1)
				{
					double prev_score = 0.0;
					for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
					{
						prev_score = std::max(prev_score, dp_tables.prefix1[*it]);
					}
					int i_prev = dp_tables.prev[i];
					if (i_prev != -1)
					{
						prev_score = std::max(prev_score, dp_tables.prefix1_L[i_prev]);
					}
					dp_tables.prefix1[i] = dp_tables.common_scores[i] + prev_score;
				}

				if(i != 0)
				{
					dp_tables.prefix1_L[i] = std::max(dp_tables.prefix1_L[i-1], dp_tables.prefix1[i]);
				}
				else
				{
					dp_tables.prefix1_L[i] = dp_tables.prefix1[i];
				}
			}

			// fill suffix1
			for(int j = 0; j < n; ++j)
			{
				int l = dp_tables.right_neighbors[N-j].size();
				dp_tables.suffix1[j].resize(l+1, std::numeric_limits<DoubleReal>::min());

				if (dp_tables.common[N-j] != -1)
				{
					double next_score = 0.0;
					for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin();
					it != dp_tables.right_jumps[N-j].end(); ++it)
					{
						next_score = std::max(next_score, dp_tables.suffix2[N-*it]+dif_penalty);
					}
					int j_next = dp_tables.next[N-j];
					if(j_next != -1)
					{
						next_score = std::max(next_score, dp_tables.suffix2_R[N-j_next]+dif_penalty);
					}
					for(std::vector<int>::const_iterator it = dp_tables.right_jumps_shifted[N-j].begin(); it != dp_tables.right_jumps_shifted[N-j].end(); ++it)
					{
						next_score = std::max(next_score, dp_tables.suffix1max[N-*it]);
					}
					j_next = dp_tables.next_shifted[N-j];
					if(j_next != -1)
					{
						next_score = std::max(next_score, dp_tables.suffix1_R[N-j_next]);
					}
					dp_tables.suffix1[j][0] = dp_tables.common_scores[N-j] + next_score;

					for(int s = 0; s < l; ++s)
					{
						int j_next = dp_tables.right_neighbors[N-j][s];
						if (dp_tables.common[j_next] == -1)
						{
							continue;
						}
						double penalty = 0.0;
						double y = dp_tables.peaks2[dp_tables.common[N-j]]+dp_tables.peaks2[dp_tables.common[j_next]];
						if(fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
						{
							penalty = sv_penalty;
						}
						double next_score = dp_tables.suffix1max[N-j_next];
						dp_tables.suffix1[j][s+1] = dp_tables.common_scores[N-j] + next_score + penalty;
					}
				}

				dp_tables.suffix1max[j] = *std::max_element(dp_tables.suffix1[j].begin(), dp_tables.suffix1[j].end());
				if(j != 0)
				{
					dp_tables.suffix1_R[j] = std::max(dp_tables.suffix1_R[j-1], dp_tables.suffix1max[j]);
				}
				else
				{
					dp_tables.suffix1_R[j] = dp_tables.suffix1max[j];
				}
			}

			// fill prefix2
			for(int i = 0; i < n; ++i)
			{
				int l = dp_tables.left_neighbors[i].size();
				dp_tables.prefix2[i].resize(l+1 ,std::numeric_limits<DoubleReal>::min());

				if (dp_tables.common_shifted[i] != -1)
				{// s = 0
					double prev_score = 0.0;
					for (std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin();it != dp_tables.left_jumps[i].end(); ++it)
					{
						prev_score = std::max(prev_score, dp_tables.prefix1[*it]+dif_penalty);
					}
					int i2 = dp_tables.prev[i];
					if (i2 != -1)
					{
						prev_score = std::max(prev_score, dp_tables.prefix1_L[i2]+dif_penalty);
					}
					for(std::vector<int>::const_iterator it = dp_tables.left_jumps_shifted[i].begin(); it != dp_tables.left_jumps_shifted[i].end(); ++it)
					{
						prev_score = std::max(prev_score, dp_tables.prefix2max[*it]);
					}
					i2 = dp_tables.prev_shifted[i];
					if (i2 != -1)
					{
						prev_score = std::max(prev_score, dp_tables.prefix2_L[i2]);
					}
					dp_tables.prefix2[i][0] = dp_tables.common_shifted_scores[i] + prev_score;

					// s > 0
					for (int s = 0; s < l; ++s)
					{
						int i2 = dp_tables.left_neighbors[i][s];
						if(dp_tables.common_shifted[i2] == -1)
						{
							continue;
						}
						double penalty = 0.0;
						double y = dp_tables.peaks2[dp_tables.common_shifted[i]]+dp_tables.peaks2[dp_tables.common_shifted[i2]];
						if(fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
						{
							penalty = sv_penalty;
						}
						prev_score = dp_tables.prefix2max[i2];
						dp_tables.prefix2[i][s+1] = dp_tables.common_shifted_scores[i] + prev_score + penalty;
					}
				}

				dp_tables.prefix2max[i] = *std::max_element(dp_tables.prefix2[i].begin(), dp_tables.prefix2[i].end());
				if (i != 0)
				{
					dp_tables.prefix2_L[i] = std::max(dp_tables.prefix2_L[i-1], dp_tables.prefix2max[i]);
				}
				else
				{
					dp_tables.prefix2_L[i] = dp_tables.prefix2max[i];
				}
			}

			// fill suffix2
			for(int j = 0; j < n; ++j)
			{
				if(dp_tables.common_shifted[N-j] != -1)
				{
					double next_score = 0.0;
					for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
					{
						next_score = std::max(next_score, dp_tables.suffix2[N-*it]);
					}
					int j_next = dp_tables.next[N-j];
					if(j_next != -1)
					{
						next_score = std::max(next_score, dp_tables.suffix2_R[N-j_next]);
					}
					dp_tables.suffix2[j] = dp_tables.common_shifted_scores[N-j] + next_score;
				}
				if(j != 0)
				{
					dp_tables.suffix2_R[j] = std::max(dp_tables.suffix2_R[j-1], dp_tables.suffix2[j]);
				}
				else
				{
					dp_tables.suffix2_R[j] = dp_tables.suffix2[j];
				}
			}

			// fill D1/M1_R
			for (int i = 0; i < n; ++i)
			{
				if (dp_tables.common[i] == -1)
				{
					if (i != 0)
					{
						dp_tables.M1_L[i] = dp_tables.M1_L[i-1];
					}
					continue;
				}

				for(int j = 0; j < n; ++j)
				{
					if(dp_tables.common_shifted[N-j] == -1)
					{
						if(i != 0)
						{
							dp_tables.M1_L[i][j] = dp_tables.M1_L[i-1][j];
						}
						if(j != 0)
						{
							dp_tables.M1_R[i][j] = dp_tables.M1_R[i][j-1];
						}
						continue;
					}

					double penalty = 0.0;
					double x = dp_tables.peaks[i]+dp_tables.peaks[N-j];
					double y = dp_tables.peaks2[dp_tables.common[i]]+dp_tables.peaks2[dp_tables.common_shifted[N-j]];
					if(fabs(x-(dp_tables.pm_s1/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance || fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
					{
						penalty = sv_penalty;
					}

					if(i >= j)
					{
						double prev_score = dp_tables.suffix2[j]+dif_penalty;
						// If we choose suffix2[j], then we pay the PTM penalty - the penalty is paid only once!
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
						{
							prev_score = std::max(prev_score, dp_tables.D1[*it][j]);
						}
						int i2 = dp_tables.prev[i];
						if(i2 != -1)
						{
							prev_score = std::max(prev_score, dp_tables.M1_L[i2][j]);
						}
						dp_tables.D1[i][j] = dp_tables.common_scores[i] + prev_score + penalty;
					}
					else
					{
						double next_score = dp_tables.prefix1[i]+dif_penalty;
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
						{
							next_score = std::max(next_score, dp_tables.D1[i][N-*it]);
						}
						int j2 = dp_tables.next[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, dp_tables.M1_R[i][N-j2]);
						}
						dp_tables.D1[i][j] = dp_tables.common_shifted_scores[N-j] + next_score + penalty;
					}

					// Compute M1_L
					if(i != 0)
					{
						dp_tables.M1_L[i][j] = std::max(dp_tables.M1_L[i-1][j], dp_tables.D1[i][j]);
					}
					else
					{
						dp_tables.M1_L[i][j] = dp_tables.D1[i][j];
					}

					// Compute M1_R
					if(j != 0)
					{
						dp_tables.M1_R[i][j] = std::max(dp_tables.M1_R[i][j-1], dp_tables.D1[i][j]);
					}
					else
					{
						dp_tables.M1_R[i][j] = dp_tables.D1[i][j];
					}
				}
			}

			// fill D2/D2max/M2_L
			for (int i = 0; i < n; ++i)
			{
				int l = dp_tables.left_neighbors[i].size();
				for(int j = 0; j < n; ++j)
				{
					dp_tables.D2[i][j].resize(l+1, std::numeric_limits<DoubleReal>::min());
					dp_tables.M2_R[i][j].resize(l+1, std::numeric_limits<DoubleReal>::min());
				}
				if(dp_tables.common_shifted[i] == -1)
				{
					if (i > 0)
					{
						dp_tables.M2_L[i] = dp_tables.M2_L[i-1];
					}
					continue;
				}

				for(int j = 0; j < n; ++j)
				{
					if(dp_tables.common_shifted[N-j] == -1)
					{
						if(i != 0)
						{
							dp_tables.M2_L[i][j] = dp_tables.M2_L[i-1][j];
						}
						if(j != 0)
						{
							dp_tables.M2_R[i][j] = dp_tables.M2_R[i][j-1];
						}
						continue;
					}

					double penalty = 0.0;
					double x = dp_tables.peaks[i]+dp_tables.peaks[N-j];
					double y = dp_tables.peaks2[dp_tables.common_shifted[i]]+dp_tables.peaks2[dp_tables.common_shifted[N-j]];
					if(fabs(x-(dp_tables.pm_s1/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance || fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
					{
						penalty = sv_penalty;
					}

					if (i > j)
					{// s = 0
						double prev_score = dp_tables.suffix2[j];
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
						{
							prev_score = std::max(prev_score, dp_tables.D1[*it][j]);
						}
						int i2 = dp_tables.prev[i];
						if(i2 != -1)
						{
							prev_score = std::max(prev_score, dp_tables.M1_L[i2][j]);
						}
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps_shifted[i].begin(); it != dp_tables.left_jumps_shifted[i].end(); ++it)
						{
							prev_score = std::max(prev_score, dp_tables.D2max[*it][j]);
						}
						i2 = dp_tables.prev_shifted[i];
						if (i2 != -1)
						{
							prev_score = std::max(prev_score, dp_tables.M2_L[i2][j]);
						}
						dp_tables.D2[i][j][0] = dp_tables.common_shifted_scores[i] + prev_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int i2 = dp_tables.left_neighbors[i][s];
							if(dp_tables.common_shifted[i2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = dp_tables.peaks2[dp_tables.common_shifted[i]]+dp_tables.peaks2[dp_tables.common_shifted[i2]];
							if(fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += sv_penalty;
							}
							prev_score = dp_tables.D2max[i2][j];
							dp_tables.D2[i][j][s+1] = dp_tables.common_shifted_scores[i] + prev_score + penalty2;
						}
					}
					else
					{// i <= j
						// s = 0
						double next_score = dp_tables.prefix2[i][0];
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
						{
							next_score = std::max(next_score, dp_tables.D2[i][N-*it][0]);
						}
						int j2 = dp_tables.next[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, dp_tables.M2_R[i][N-j2][0]);
						}
						dp_tables.D2[i][j][0] = dp_tables.common_shifted_scores[N-j] + next_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int i2 = dp_tables.left_neighbors[i][s];
							if(dp_tables.common_shifted[i2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = dp_tables.peaks2[dp_tables.common_shifted[i2]]+dp_tables.peaks2[dp_tables.common_shifted[N-j]];
							if(fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += sv_penalty;
							}
							double next_score = dp_tables.prefix2[i][s+1];
							for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
							{
								next_score = std::max(next_score, dp_tables.D2[i][N-*it][s+1]);
							}
							int j2 = dp_tables.next[N-j];
							if(j2 != -1)
							{
								next_score = std::max(next_score, dp_tables.M2_R[i][N-j2][s+1]);
							}
							dp_tables.D2[i][j][s+1] = dp_tables.common_shifted_scores[N-j] + next_score + penalty2;
						}
					}

					// Compute D2max
					dp_tables.D2max[i][j] = *std::max_element(dp_tables.D2[i][j].begin(), dp_tables.D2[i][j].end());

					// Compute M2_L
					if(i != 0)
					{
						dp_tables.M2_L[i][j] = std::max(dp_tables.M2_L[i-1][j], dp_tables.D2max[i][j]);
					}
					else
					{
						dp_tables.M2_L[i][j] = dp_tables.D2max[i][j];
					}

					// compute M2_R
					if(j != 0)
					{
						for(int s = 0; s < l+1; ++s)
						{
							dp_tables.M2_R[i][j][s] = std::max(dp_tables.M2_R[i][j-1][s], dp_tables.D2[i][j][s]);
						}
					}
					else
					{
						dp_tables.M2_R[i][j] = dp_tables.D2[i][j];
					}
				}
			}

			// fill D3/M3_R
			for(int i = 0; i < n; ++i)
			{
				for(int j = 0; j < n; ++j)
				{
					int l = dp_tables.right_neighbors[N-j].size();
					dp_tables.D3[i][j].resize(l+1, std::numeric_limits<DoubleReal>::min());
					dp_tables.M3_L[i][j].resize(l+1, std::numeric_limits<DoubleReal>::min());
				}
				if(dp_tables.common[i] == -1)
				{
					if(i > 0)
					{
						dp_tables.M3_L[i] = dp_tables.M3_L[i-1];
					}
					continue;
				}

				for(int j = 0; j < n; ++j)
				{
					if(dp_tables.common[N-j] == -1)
					{
						if(i != 0)
						{
							dp_tables.M3_L[i][j] = dp_tables.M3_L[i-1][j];
						}
						if(j != 0)
						{
							dp_tables.M3_R[i][j] = dp_tables.M3_R[i][j-1];
						}
						continue;
					}

					int l = dp_tables.right_neighbors[N-j].size();
					double penalty = 0.0;
					double x = dp_tables.peaks[i]+dp_tables.peaks[N-j];
					double y = dp_tables.peaks2[dp_tables.common[i]]+dp_tables.peaks2[dp_tables.common[N-j]];
					if(fabs(x-(dp_tables.pm_s1/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance || fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
					{
						penalty = sv_penalty;
					}

					if(i >= j)
					{ // s = 0
						double prev_score = dp_tables.suffix1[j][0];
						for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
						{
							prev_score = std::max(prev_score, dp_tables.D3[*it][j][0]);
						}
						int i2 = dp_tables.prev[i];
						if(i2 != -1)
						{
							prev_score = std::max(prev_score, dp_tables.M3_L[i2][j][0]);
						}
						dp_tables.D3[i][j][0] = dp_tables.common_scores[i] + prev_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int j2 = dp_tables.right_neighbors[N-j][s];
							if(dp_tables.common[j2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = dp_tables.peaks2[dp_tables.common[i]]+dp_tables.peaks2[dp_tables.common[j2]];
							if(fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += sv_penalty;
							}
							double prev_score = dp_tables.suffix1[j][s+1];
							for(std::vector<int>::const_iterator it = dp_tables.left_jumps[i].begin(); it != dp_tables.left_jumps[i].end(); ++it)
							{
								prev_score = std::max(prev_score, dp_tables.D3[*it][j][s+1]);
							}
							int i2 = dp_tables.prev[i];
							if(i2 != -1)
							{
								prev_score = std::max(prev_score, dp_tables.M3_L[i2][j][s+1]);
							}
							dp_tables.D3[i][j][s+1] = dp_tables.common_scores[i] + prev_score + penalty2;
						}
					}
					else
					{ // s = 0
						double next_score = dp_tables.prefix1[i];
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps[N-j].begin(); it != dp_tables.right_jumps[N-j].end(); ++it)
						{
							next_score = std::max(next_score, dp_tables.D1[i][N-*it]);
						}
						int j2 = dp_tables.next[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, dp_tables.M1_R[i][N-j2]);
						}
						for(std::vector<int>::const_iterator it = dp_tables.right_jumps_shifted[N-j].begin(); it != dp_tables.right_jumps_shifted[N-j].end(); ++it)
						{
							next_score = std::max(next_score, dp_tables.D3max[i][N-*it]);
						}
						j2 = dp_tables.next_shifted[N-j];
						if (j2 != -1)
						{
							next_score = std::max(next_score, dp_tables.M3_R[i][N-j2]);
						}
						dp_tables.D3[i][j][0] = dp_tables.common_scores[N-j] + next_score + penalty;

						// s > 0
						for(int s = 0; s < l; ++s)
						{
							int j2 = dp_tables.right_neighbors[N-j][s];
							if(dp_tables.common[j2] == -1)
							{
								continue;
							}
							double penalty2 = penalty;
							double y = dp_tables.peaks2[dp_tables.common[N-j]]+dp_tables.peaks2[dp_tables.common[j2]];
							if(fabs(y-(dp_tables.pm_s2/* +18 */+Constants::PROTON_MASS_U)) < peak_tolerance)
							{
								penalty2 += sv_penalty;
							}
							double next_score = dp_tables.D3max[i][N-j2];
							dp_tables.D3[i][j][s+1] = dp_tables.common_scores[N-j] + next_score + penalty2;
						}
					}

					// Compute D3max
					dp_tables.D3max[i][j] = *std::max_element(dp_tables.D3[i][j].begin(), dp_tables.D3[i][j].end());

					// Compute M3_L
					if (i != 0)
					{
						for(int s = 0; s < l+1; ++s)
						{
							dp_tables.M3_L[i][j][s] = std::max(dp_tables.M3_L[i-1][j][s], dp_tables.D3[i][j][s]);
						}
					}
					else
					{
						dp_tables.M3_L[i][j] = dp_tables.D3[i][j];
					}

					// Compute M3_R
					if(j != 0)
					{
						dp_tables.M3_R[i][j] = std::max(dp_tables.M3_R[i][j-1], dp_tables.D3max[i][j]);
					}
					else
					{
						dp_tables.M3_R[i][j] = dp_tables.D3max[i][j];
					}
				}
			}

			// find best score
			dp_tables.best_t = ASA_DP_Tables::cell_Invalid;
			for(int i = 0; i < n; ++i)
			{
				for(std::vector<int>::const_iterator it = dp_tables.right_jumps[i].begin(); it != dp_tables.right_jumps[i].end(); ++it)
				{
					int j = N-*it;
					if(j <= n-1)
					{
						if(dp_tables.best_score < dp_tables.D1[i][j])
						{
							dp_tables.best_score = dp_tables.D1[i][j];
							dp_tables.best_i = i;
							dp_tables.best_j = j;
							dp_tables.best_s = 0;
							dp_tables.best_t = ASA_DP_Tables::cell_D1;
						}
						if(dp_tables.best_score < dp_tables.D2max[i][j])
						{
							dp_tables.best_score = dp_tables.D2max[i][j];
							dp_tables.best_i = i;
							dp_tables.best_j = j;
							dp_tables.best_s = (std::max_element(dp_tables.D2[i][j].begin(),dp_tables.D2[i][j].end())-dp_tables.D2[i][j].begin());
							dp_tables.best_t = ASA_DP_Tables::cell_D2;
						}
						if(dp_tables.best_score < dp_tables.D3max[i][j])
						{
							dp_tables.best_score = dp_tables.D3max[i][j];
							dp_tables.best_i = i;
							dp_tables.best_j = j;
							dp_tables.best_s = (std::max_element(dp_tables.D3[i][j].begin(),dp_tables.D3[i][j].end())-dp_tables.D3[i][j].begin());
							dp_tables.best_t = ASA_DP_Tables::cell_D3;
						}
					}
				}

				int j0 = dp_tables.next[i];
				if(j0 != -1)
				{
					int j = std::min(n-1, N-j0);
					if(dp_tables.best_score < dp_tables.M1_R[i][j])
					{
						dp_tables.best_score = dp_tables.M1_R[i][j];
						dp_tables.best_i = i;
						dp_tables.best_j = j;
						dp_tables.best_s = 0;
						dp_tables.best_t =ASA_DP_Tables::cell_M1_R;
					}

					double tmp = *std::max_element(dp_tables.M2_R[i][j].begin(), dp_tables.M2_R[i][j].end());
					if (dp_tables.best_score < tmp)
					{
						dp_tables.best_score = tmp;
						dp_tables.best_i = i;
						dp_tables.best_j = j;
						dp_tables.best_s = (std::max_element(dp_tables.M2_R[i][j].begin(),dp_tables.M2_R[i][j].end())-dp_tables.M2_R[i][j].begin());
						dp_tables.best_t = ASA_DP_Tables::cell_M2_R;
					}

					if (dp_tables.best_score < dp_tables.M3_R[i][j])
					{
						dp_tables.best_score = dp_tables.M3_R[i][j];
						dp_tables.best_i = i;
						dp_tables.best_j = j;
						dp_tables.best_s = 0;
						dp_tables.best_t = ASA_DP_Tables::cell_M3_R;
					}
				}

				if(dp_tables.best_score < dp_tables.prefix1[i])
				{
					dp_tables.best_score = dp_tables.prefix1[i];
					dp_tables.best_i = i;
					dp_tables.best_j = 0;
					dp_tables.best_s = 0;
					dp_tables.best_t = ASA_DP_Tables::cell_prefix1;
				}

				if(dp_tables.best_score < dp_tables.suffix1max[i])
				{
					dp_tables.best_score = dp_tables.suffix1max[i];
					dp_tables.best_j = i;
					dp_tables.best_i = 0;
					dp_tables.best_s = (std::max_element(dp_tables.suffix1[i].begin(),dp_tables.suffix1[i].end())-dp_tables.suffix1[i].begin());
					dp_tables.best_t = ASA_DP_Tables::cell_suffix1;
				}

				if(dp_tables.best_score < dp_tables.prefix2max[i])
				{
					dp_tables.best_score = dp_tables.prefix2max[i];
					dp_tables.best_i = i;
					dp_tables.best_j = 0;
					dp_tables.best_s = (std::max_element(dp_tables.prefix2[i].begin(),dp_tables.prefix2[i].end())-dp_tables.prefix2[i].begin());
					dp_tables.best_t = ASA_DP_Tables::cell_prefix2;
				}

				if(dp_tables.best_score < dp_tables.suffix2[i])
				{
					dp_tables.best_score = dp_tables.suffix2[i];
					dp_tables.best_j = i;
					dp_tables.best_i = 0;
					dp_tables.best_s = 0;
					dp_tables.best_t = ASA_DP_Tables::cell_suffix2;
				}
			}

			return;
		}
		///


		/**
			@brief Method to calculate a antisymmetric alignment

			@param res_1 Spectrum containing the aligned peaks from s1
			@param res_2 Spectrum containing the aligned peaks from s2
			@param score containing the score (best_score/parent_mass s1)
			@param mod_pos the position in res_1 that has to be 'modified' (even if it is the aligned counterpart in res2 that had a subtraction)
			@param s2 the first input spectrum
			@param s2 the second input spectrum

			This method must have the first spectrum have a lesser or equal parent mass than the second spectrum, else a error is thrown. Also the first spectrum must be symmetric, so synthetic peaks with zero intensity might be inserted. Input does not. Also the resulting spectra are intended for consensus making which takes a list of spectra.
		*/
		void getAntisymmetricAlignment(MSSpectrum<PeakT>& res_1, MSSpectrum<PeakT>& res_2, DoubleReal& score, DoubleReal& mod_pos, const MSSpectrum<PeakT>& s1, const MSSpectrum<PeakT>& s2) const
		{
			//~ float sameVertexPenalty=-5, float ptmPenalty=-5
			//~ sameVertexPenalty = -1000000, ptmPenalty = -200; //makes the test fail!
			//~ sameVertexPenalty = 0, ptmPenalty = 0;

			Real peak_tolerance = param_.getValue("peak_tolerance");
			DoubleReal pm_tolerance = param_.getValue("parentmass_tolerance");
			DoubleReal min_dist = param_.getValue("min_dist");

			/// @attention uncharged mass!
			DoubleReal pm_s1 (s1.getPrecursors().front().getUnchargedMass());
			DoubleReal pm_s2 (s2.getPrecursors().front().getUnchargedMass());
			DoubleReal pm_diff = (pm_s2-pm_s1);

			if(pm_diff <= -pm_tolerance)
			{
				throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "prerequisite is the s1-precursor mz is not greater than s2-precursor mz");
			}

			res_1.clear(true);
			res_2.clear(true);
			//~ res_1.clear(); res_1.getFloatDataArrays().clear();res_1.getIntegerDataArrays().clear();res_1.getStringDataArrays().clear();
			//~ res_2.clear(); res_2.getFloatDataArrays().clear();res_2.getIntegerDataArrays().clear();res_2.getStringDataArrays().clear();
			score = -1; mod_pos = -1;

			if(s1.empty() or s2.empty())
			{
				return;
			}

			res_1 = s1;
			res_2 = s2;

			/// @improvement rise max. number of AA jumped possible (now just one AA)
			//jump_masses[i] holds internal mass of aa i - make unique in resolution tolerance (e.g. cause of Q & K)
			const ResidueDB* res_db = ResidueDB::getInstance();
			std::vector<Real> jump_masses;
			for(ResidueDB::ResidueConstIterator it = res_db->beginResidue(); it != res_db->endResidue(); ++it)
			{
				DoubleReal w = (*it)->getMonoWeight() - (*it)->getInternalToFullMonoWeight();
				/// @important: this makes only singly charged ion peaks to potential jump matches and therefor potential alignment-partners (so a test is only reasonable with singly charged ion peaks in da spec)
				if(w>0)
				{
					jump_masses.push_back(w);
				}
			}
			//~ sort jump_masses
			std::sort(jump_masses.begin(),jump_masses.end());
			std::vector<Real>::iterator end_new = std::unique(jump_masses.begin(),jump_masses.end(), EqualInTolerance<Real>(peak_tolerance));
			jump_masses.resize( end_new - jump_masses.begin() );
			DoubleReal jumps_supremum=0;
			jumps_supremum = jump_masses.back()+2*peak_tolerance+0.00001;

//~ remove unneccessary peaks
			//~ charge deconvolution to 'singly charged only' for identified spectra
			filterHighCharges(res_1,2,3);

			//~ change s1 to align antisymmetrical into res_1
			res_1 = getSymetricSpectrum(res_1);

			//~ prune to matching/ shift-matching peaks
			matchPrunePeaks(res_1, res_2, pm_diff);

//~ fill the DP tables
			//~ begin with information collection
			Size m((Size)(((DoubleReal)res_1.getIntegerDataArrays()[1].size())/2.0));
			ASA_DP_Tables dp_tables(res_1.size(), m+1, res_2.size(), (pm_s1+Constants::PROTON_MASS_U), (pm_s2+Constants::PROTON_MASS_U));

			for(Size i=0;i<res_1.size();++i)
			{
				dp_tables.peaks[i]=res_1[i].getMZ();
			}
			for(Size i=0;i<res_2.size();++i) ///@attention was only running up to res_1.size()
			{
				dp_tables.peaks2[i]=res_2[i].getMZ();
			}

			for(Size i=0;i<res_1.size();++i)
			{
				if(res_1.getIntegerDataArrays()[1][i]>-1)
				{
					dp_tables.common[i] = res_1.getIntegerDataArrays()[1][i];
					dp_tables.common_scores[i] = res_1[i].getIntensity()+s2[dp_tables.common[i]].getIntensity();
				}
				else if(res_1.getIntegerDataArrays()[2][i]>-1)
				{
					dp_tables.common_shifted[i] = res_1.getIntegerDataArrays()[2][i];
					dp_tables.common_shifted_scores[i] = res_1[i].getIntensity()+s2[dp_tables.common_shifted[i]].getIntensity();
				}
			}

			dp_tables.delta = std::min(std::max(min_dist-peak_tolerance,pm_diff+peak_tolerance),min_dist-peak_tolerance*2.0f);
			dp_tables.delta2 = jumps_supremum-peak_tolerance;

			for(Size i = 0; i < res_1.size(); ++i)
			{
				DoubleReal lo1 = res_1[i].getMZ()-dp_tables.delta2;
				if(lo1<=0)
				{
					dp_tables.prev[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_lo1 = res_1.MZEnd(res_1.begin(),lo1,res_1.begin()+i);
					it_lo1!=res_1.begin()?dp_tables.prev[i]=((it_lo1-res_1.begin())-1):dp_tables.prev[i]=-1;
				}
				DoubleReal lo2 = res_1[i].getMZ()-std::max(dp_tables.delta,dp_tables.delta2);
				if(lo2<=0)
				{
					dp_tables.prev_shifted[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_lo2 = res_1.MZEnd(lo2);
					it_lo2!=res_1.begin()?dp_tables.prev_shifted[i]=((it_lo2-res_1.begin())-1):dp_tables.prev_shifted[i]=-1;
				}

				DoubleReal hi1 = res_1[i].getMZ()+ dp_tables.delta2;
				if(hi1>=res_1.rbegin()->getMZ())
				{
					dp_tables.next[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_hi1 = res_1.MZBegin(res_1.begin()+i,hi1,res_1.end());
					it_hi1!=res_1.end()?dp_tables.next[i]=(it_hi1-res_1.begin()):dp_tables.next[i]=-1;
				}
				DoubleReal hi2 = res_1[i].getMZ()+std::max(dp_tables.delta,dp_tables.delta2);
				if(hi2>=res_1.rbegin()->getMZ())
				{
					dp_tables.next_shifted[i] = -1;
				}
				else
				{
					typename SpectrumType::Iterator it_hi2 = res_1.MZBegin(res_1.begin()+i,hi2,res_1.end());
					it_hi2!=res_1.end()?dp_tables.next_shifted[i]=(it_hi2-res_1.begin()):dp_tables.next_shifted[i]=-1;
				}

				getAASteps(i, res_1, jump_masses, dp_tables.left_jumps[i], dp_tables.right_jumps[i]);
				Size neighbor_count=0, jumps_count=0;
				dp_tables.left_neighbors[i].resize(dp_tables.left_jumps[i].size());
				dp_tables.left_jumps_shifted[i].resize(dp_tables.left_jumps[i].size());
				for(Size j=0; j<dp_tables.left_jumps[i].size(); ++j)
				{
					if(dp_tables.peaks[i]-dp_tables.peaks[dp_tables.left_jumps[i][j]]<=dp_tables.delta)
					{
						dp_tables.left_neighbors[i][neighbor_count++]=dp_tables.left_jumps[i][j];
					}
					if(dp_tables.peaks[i]-dp_tables.peaks[dp_tables.left_jumps[i][j]]>=dp_tables.delta)
					{
						dp_tables.left_jumps_shifted[i][jumps_count++]=dp_tables.left_jumps[i][j];
					}
				}
				dp_tables.left_neighbors[i].resize(neighbor_count);
				dp_tables.left_jumps_shifted[i].resize(jumps_count);

				neighbor_count=0, jumps_count=0;
				dp_tables.right_neighbors[i].resize(dp_tables.right_jumps[i].size());
				dp_tables.right_jumps_shifted[i].resize(dp_tables.right_jumps[i].size());
				for(Size j=0; j<dp_tables.right_jumps[i].size(); ++j)
				{
					if(dp_tables.peaks[dp_tables.right_jumps[i][j]]-dp_tables.peaks[i]<=dp_tables.delta)
					{
						dp_tables.right_neighbors[i][neighbor_count++]=dp_tables.right_jumps[i][j];
					}
					if(dp_tables.peaks[dp_tables.right_jumps[i][j]]-dp_tables.peaks[i]>=dp_tables.delta)
					{
						dp_tables.right_jumps_shifted[i][jumps_count++]=dp_tables.right_jumps[i][j];
					}
				}
				dp_tables.right_neighbors[i].resize(neighbor_count);
				dp_tables.right_jumps_shifted[i].resize(jumps_count);
			}

			res_1.clear(false);
			res_2.clear(false);
			//~ implicitly copied (and not erased):
			//~ Precursors,RT,MSLevel, MetaDataArrays

			/// @improvement find a good penalty adjustment
			DoubleReal sv_penalty = (DoubleReal)param_.getValue("sv_penalty");
			DoubleReal dif_penalty = (DoubleReal)param_.getValue("dif_penalty");

			//~ go on with tables and traceback, result is aligned peaks list
			calculateAlignementPath(dp_tables, peak_tolerance, sv_penalty, dif_penalty);
			std::pair< std::vector<int>,std::vector<int> > asym_align;
			asym_align = traceback(dp_tables, dif_penalty, score);

			Size num_aligned_peaks = asym_align.first.size()+asym_align.second.size();

			if(num_aligned_peaks==0)
			{
				score = -1; mod_pos = -1;
				res_1.clear(true);
				res_2.clear(true);
				//~ res_1.clear(); res_1.getFloatDataArrays().clear();res_1.getIntegerDataArrays().clear();res_1.getStringDataArrays().clear();
				//~ res_2.clear(); res_2.getFloatDataArrays().clear();res_2.getIntegerDataArrays().clear();res_2.getStringDataArrays().clear();
				return;
			}

			typename SpectrumType::IntegerDataArray ida_synthetic; ida_synthetic.resize(num_aligned_peaks,0);
			res_1.getIntegerDataArrays().insert(res_1.getIntegerDataArrays().begin(),ida_synthetic);
			(res_1.getIntegerDataArrays().begin())->setName("synthetic peaks");
			res_1.getIntegerDataArrays().insert(res_1.getIntegerDataArrays().begin()+1,ida_synthetic);
			(res_1.getIntegerDataArrays().begin()+1)->setName("modification position");
			res_2.getIntegerDataArrays().insert(res_2.getIntegerDataArrays().begin(), ida_synthetic);
			(res_2.getIntegerDataArrays().begin())->setName("synthetic peaks");
			res_2.getIntegerDataArrays().insert(res_2.getIntegerDataArrays().begin()+1, ida_synthetic);
			(res_2.getIntegerDataArrays().begin()+1)->setName("modification position");

			if(asym_align.first.size()==0)
			{
				mod_pos=0;  // No peaks in idx1 => mod at the start
				(res_1.getIntegerDataArrays().begin()+1)->at(0)=-1;
				(res_2.getIntegerDataArrays().begin()+1)->at(0)=-1;
			}
			else
			{
				if(asym_align.second.size()==0)
				{
					mod_pos = /* pm_s1 */-1;  // No peaks in idx2 => mod at the end
					(res_1.getIntegerDataArrays().begin()+1)->at(num_aligned_peaks-1)=-1;
					(res_2.getIntegerDataArrays().begin()+1)->at(num_aligned_peaks-1)=-1;
				}
				else
				{
					mod_pos=dp_tables.peaks[asym_align.second[0]];  // Otherwise mod was placed at the first mass shifted aligned pair (but unshifted version as value)
					(res_1.getIntegerDataArrays().begin()+1)->at(asym_align.first.size())=1;
					(res_2.getIntegerDataArrays().begin()+1)->at(asym_align.first.size())=1;
				}
			}

			typename SpectrumType::IntegerDataArray& ida_ref = *(res_1.getIntegerDataArrays().begin()+2);
			for(Size i=0; i<asym_align.first.size(); ++i)
			{
				PeakT tmp_1, tmp_2;
				tmp_1.setMZ( dp_tables.peaks[asym_align.first[i]] );
				if(dp_tables.common[asym_align.first[i]]>=0)
				{
					tmp_2.setMZ(dp_tables.peaks2[dp_tables.common[asym_align.first[i]]]);
					tmp_2.setIntensity(s2[dp_tables.common[asym_align.first[i]]].getIntensity());
					res_2.getIntegerDataArrays().front().push_back(0);
				}
				else
				{
					tmp_2.setMZ(tmp_1.getMZ());
					//~ intensity stays default constructed 0
					res_2.getIntegerDataArrays().front().push_back(1);
				}
				tmp_1.setIntensity(dp_tables.common_scores[asym_align.first[i]] - tmp_2.getIntensity());
				(ida_ref[asym_align.first[i]]>0)?res_2.getIntegerDataArrays().front()[i]=1:res_2.getIntegerDataArrays().front()[i]=0;
				res_1.push_back(tmp_1);
				res_2.push_back(tmp_2);
			}

			for(Size i=0; i<asym_align.second.size(); ++i)
			{
				PeakT tmp_1, tmp_2;
				tmp_1.setMZ( dp_tables.peaks[asym_align.second[i]] );
				if(dp_tables.common_shifted[asym_align.second[i]]>=0)
				{
					tmp_2.setMZ( dp_tables.peaks2[dp_tables.common_shifted[asym_align.second[i]]] );
					tmp_2.setIntensity( s2[dp_tables.common_shifted[asym_align.second[i]]].getIntensity() );
					res_2.getIntegerDataArrays().front().push_back(0);
				}
				else
				{
					tmp_2.setMZ(tmp_1.getMZ() + pm_diff);
					//~ intensity stays default constructed 0
					res_2.getIntegerDataArrays().front().push_back(1);
				}
				tmp_1.setIntensity( dp_tables.common_shifted_scores[asym_align.second[i]] - tmp_2.getIntensity());
				(ida_ref[asym_align.second[i]]>0)?res_2.getIntegerDataArrays().front()[asym_align.first.size()+i]=1:res_2.getIntegerDataArrays().front()[asym_align.first.size()+i]=0;
				res_1.push_back(tmp_1);
				res_2.push_back(tmp_2);
			}

		}
		///
	};


}
#endif //OPENMS_COMPARISON_SPECTRA_ANTISYMMETRICALIGNMENT_H
