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
#ifndef OPENMS_COMPARISON_SPECTRA_STARCLUSTERS_H
#define OPENMS_COMPARISON_SPECTRA_STARCLUSTERS_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>
#include <list>
#include <map>
#include <set>
#include <limits>
#include <utility>
#include <algorithm>

#define STAR_DEBUG
#undef  STAR_DEBUG

namespace OpenMS
{
	class Peak1D;

	/**
		@brief compares the correlation of two spectra

		@htmlinclude OpenMS_StarClusters.parameters

		@ingroup SpectraClustering
	*/
	template <typename PeakT = Peak1D>
	class StarClusters
		: public DefaultParamHandler
	{
	public:
		typedef PeakT PeakType;
		typedef MSSpectrum<PeakType> SpectrumType;
		typedef typename MSSpectrum<PeakType>::ConstIterator ConstSpectrumIterator;
		typedef typename MSSpectrum<PeakType>::Iterator SpectrumIterator;

	private:
		std::map< Size, std::set<Size> > indices_of_i_in_pairs_; // map key is index i to spec s in original experiment, map value is all indices to s' in pairs (s' is a subset of peaks from s - the aligned ones )
		std::map< Size, std::set<Size> > stars_;
		MSExperiment<PeakType> consensuses_;
		std::vector<Size> consensuses_indices_;
		MSExperiment<PeakType>& experiment_original_;
		MSExperiment<PeakType>& pairs_;
		std::vector< std::pair<Size,Size> >& aligned_;
		std::list<DoubleReal>& modification_positions_;

	public:

		// @name Constructors and Destructors
		// @{
		/**
			@brief Constructor
			@param original a MSExperiment reference to the original input
			@param pairs a MSExperiment reference to chosen spectra pairs from original aligned (containing only the aligned peaks) (due to compab. with ./specalign)
			@param aligned the list of the pairs determined by an alignment (indices in pairs, aligned and modification_positions correspond)
			@param modification_positions a list containing the modification positions of the paired spectra in [da] (indices in pairs, aligned and modification_positions correspond)

			@param original associated experiment in original (but filtered after reading in and maybe in prm-version)
		*/
		StarClusters(MSExperiment<PeakType>& original, MSExperiment<PeakType>& pairs, std::vector< std::pair<Size,Size> >& aligned, std::list<DoubleReal>& modification_positions)
			: DefaultParamHandler("StarClusters"), experiment_original_(original), pairs_(pairs), aligned_(aligned), modification_positions_(modification_positions), consensuses_(),consensuses_indices_(),stars_(), indices_of_i_in_pairs_()
		{
			defaults_.setValue("peak_tolerance", 0.3, "Defines the absolut (in Da) peak tolerance");
			defaults_.setValue("parentmass_tolerance", 3.0, "Defines the absolut (in Da) parent mass tolerance");
			defaults_.setValue("min_shift", 0.0, "Defines the minimal absolut (in Da) shift between the two spectra");
			defaults_.setValue("max_shift", 150.0, "Defines the maximal absolut (in Da) shift between the two spectra");
			defaults_.setValue("min_dist", 57.0, "Defines the minimal distance (in Da) between the two peaks of one sort that may be connected in a sparse matching");
			defaults_.setValue("correlation_scoring", "intensity", "If intensity, correlation scores on basis of matched intensity values, if matchnumber correlation scores solely on basis of the number of matches");
			defaults_.setValidStrings("correlation_scoring", StringList::create("intensity,matchnumber"));
			defaultsToParam_();
			this->createNetwork_();
		}

		/// copy constructor
		StarClusters(const StarClusters& source)
			: DefaultParamHandler(source)
		{
			experiment_original_ = source.experiment_original_;
			pairs_ = source.pairs_;
			aligned_ = source.aligned_;
			modification_positions_ = source.modification_positions_;
			consensuses_ = source.consensuses_;
			consensuses_indices_ = source.consensuses_indices_;
			stars_ = source.stars_;
			indices_of_i_in_pairs_ = source.indices_of_i_in_pairs_;
		}

		/// destructor
		~StarClusters()
		{
		}

		/// assignment operator
		StarClusters& operator = (const StarClusters& source)
		{
			if (this != &source)
			{
				DefaultParamHandler::operator = (source);
				experiment_original_ = source.experiment_original_;
				pairs_ = source.pairs_;
				aligned_ = source.aligned_;
				modification_positions_ = source.modification_positions_;
				consensuses_ = source.consensuses_;
				consensuses_indices_ = source.consensuses_indices_;
				stars_ = source.stars_;
				indices_of_i_in_pairs_ = source.indices_of_i_in_pairs_;
			}
			return *this;
		}

		// @}


		/**
			@brief get the consensuses to the stars

			@return consensuses

			@throw if not created
		*/
		MSExperiment<PeakType>& getStarConsensuses() const
		{
			return consensuses_;
		}

		/**
			@brief Method to calculate ...

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		void getStarIndices() const
		{
			/// @todo
		}


		/**
			@brief Method to calculate ...

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		MSExperiment<PeakType>& getAssociatedMSExperiment() const
		{
			return experiment_original_;
		}

		/**
			@brief Method to calculate ...

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		void reverseSpectrum(MSSpectrum<PeakType>& spec, bool is_prm=false)
		{
			DoubleReal pm = spec.getPrecursors().first().getMZ(); /// @improvement throw error if no precursor
			if(is_prm)
			{
				pm -= 1; /// @improvement check if peakspectra reversion needs subtraction too (charge ?!)
			}

			for(SpectrumIterator p = spec.begin(); p != spec.end(); ++p)
			{
				p->setMZ(pm - p->getMZ());
			}
			spec.sortByPosition();
		}

		/**
			@brief Method to calculate ...

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		void propagateWithIDs(std::vector<PeptideIdentification>& ids, Size hops)
		{

			std::list<Size> potential_propagation_targets;
			std::map< Size, PeptideIdentification > annotations;
			//~ prime first annotation-step with ids
			for(Size i = 0; i < ids.size(); ++i)
			{
				if(ids[i].getHits().size() > 0)
				{
					annotations[consensuses_indices_[i]] = ids[i];
				}
				for(std::set<Size>::const_iterator it_star = stars_[consensuses_indices_[i]].begin(); it_star != stars_[consensuses_indices_[i]].end(); ++it_star)
				{
					potential_propagation_targets.insert(potential_propagation_targets.end(),*it_star);
				}
			}
			potential_propagation_targets.unique();
			potential_propagation_targets.sort();
			for(Size i = 0; i < ids.size(); ++i)
			{
				potential_propagation_targets.remove(consensuses_indices_[i]);
			}

			//~ propagation in |hops| steps:
			for(Size i = 0; i < hops; ++i)
			{
				std::vector<Size> newly_annotated; /// @improvement reserve potential_propagation_targets.size() space
				for(std::list<Size>::const_iterator it_potential = potential_propagation_targets.begin(); it_potential != potential_propagation_targets.end(); ++it_potential)
				{
					const MSSpectrum<PeakType>& current_spectrum  = (experiment_original_[*it_potential]);
					PeptideIdentification annotating_neighbors;
					for(std::set<Size>::const_iterator it_neighbors = stars_[*it_potential].begin(); it_neighbors != stars_[*it_potential].end(); ++it_neighbors)
					{
						for(std::map< Size, PeptideIdentification >::const_iterator it_annotated = annotations.begin(); it_annotated != annotations.end(); ++it_annotated)
						{
							if(it_annotated->first == *it_neighbors)
							{
								/// @todo annotate and create PeptideHit: score shift mod from which neighbor
									/// @todo find pair with (potential_propagation_targets[j]) and (it_annotated->first == *it_neighbors) or annotate with original_experiment_ spectrum if! you have modification mass
								/// @todo span theor spec from it_annotated->first onto neighbor spec and score and find mod
								AASequence sequence(it_annotated->second.getHits().front().getSequence()),sequence_mod;

								std::pair<Size,Size> current_pair;
								DoubleReal current_mod_pos;
								*it_potential<*it_neighbors ? current_pair=std::pair<Size,Size>(*it_potential,*it_neighbors) : current_pair=std::pair<Size,Size>(*it_neighbors,*it_potential);
								std::vector< std::pair<Size,Size> >::iterator current_pair_it = std::find(aligned_.begin(),aligned_.end(),current_pair);
								Size current_pair_index = current_pair_it - aligned_.begin();

								DoubleReal score;
								//~ getPropagation(sequence, current_spectrum, mod_pos, sequence_mod, score); //@ todo siehe blatt im block

								UInt rank = it_annotated->second.getHits().front().getRank() +1; /// @attention rank is reflecting the hop-distance to the originating db-hit! so do not assignRanks() in PeptideIdentification!!!
								int charge = it_annotated->second.getHits().front().getCharge(); /// @improvement also remember from which neighbor the annotation came
								annotating_neighbors.insertHit(PeptideHit(score,rank,charge,sequence_mod));
							}
						}
					}
					/// @improvement find out whether select best annotating_neighbor instead of set all annotating_neighbor is better
					for(Size k = 0; k < annotating_neighbors.getHits().size(); ++k)
					{
						annotations[*it_potential].insertHit(annotating_neighbors.getHits()[k]);
					}
				}

				/// @improvement for dense networks this! will! suck! ... bigtime
				//~ update potential_propagation_targets
				for(Size k = 0; k < newly_annotated.size(); ++k)
				{
					annotations[newly_annotated[k]].sort();
					for(std::set<Size>::const_iterator it_star = stars_[newly_annotated[k]].begin(); it_star != stars_[newly_annotated[k]].end(); ++it_star)
					{
						potential_propagation_targets.insert(potential_propagation_targets.begin(),*it_star);
					}
				}
				newly_annotated.clear();
				//~ erase all annot. from potential_propagation_targets.
				potential_propagation_targets.unique();
				potential_propagation_targets.sort();
				for(std::map< Size, PeptideIdentification >::const_iterator it_annotated = annotations.begin(); it_annotated != annotations.end(); ++it_annotated)
				{
					potential_propagation_targets.remove(it_annotated->first);
				}
			}

			/// @todo copy the annotations to original_experiment and calculate a end score
		}

	private:
		/// default constructor
		StarClusters()
		{
		}


		/**
			@brief Method to calculate ...

			@param ...
			@param ...
			@return ...

			...
			@see ...
		*/
		DoubleReal matchFromAligns(MSSpectrum<PeakType>& s1, MSSpectrum<PeakType>& s2, DoubleReal shiftwindow, DoubleReal stepsize) const
		{
			DoubleReal peak_tolerance = (double)param_.getValue("peak_tolerance");

			/// @improvement for prms score initialization with 0 is no good use limits::min

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
					for(SpectrumIterator it = start; it != end; ++it)
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


		/**
			@brief will create the network as a set of starsclusters which will most probabply overlay so its a network
			...
			@see ...
		*/
		void createNetwork()
		{
			for(Size i = 0; i < aligned_.size(); ++i)
			{
				//~ if(!(stars[aligned[i].first].size())) //if it contains something it sure contains the first, this doesnt has to be added anymore
				//~ {
					//~ stars[aligned[i].first].insert(aligned[i].first);
				//~ }
				//~ if(!(stars[aligned[i].second].size())) //if it contains something it sure contains the first, this doesnt has to be added anymore
				//~ {
					//~ stars[aligned[i].second].insert(aligned[i].second);
				//~ }
				stars_[aligned_[i].first].insert(/* stars_[aligned[i].first].rbegin(), */aligned_[i].second);
				stars_[aligned_[i].second].insert(/* stars_[aligned[i].second].rbegin(), */aligned_[i].first);
				indices_of_i_in_pairs_[aligned_[i].first].insert(i*2);
				indices_of_i_in_pairs_[aligned_[i].second].insert((i*2)+1);
			}

			//orientate the specs - only neccessary if the spectra in pairs came from antisymmetric path alignment
			for(std::map< Size, std::set<Size> >::const_iterator it = indices_of_i_in_pairs_.begin(); it != indices_of_i_in_pairs_.end(); ++it)
			{

				//binning parameters: sz is resolution, sp is parentmass/res
				Real sz = (Real)param_.getValue("peak_tolerance");
				Real parentmass_tolerance = (Real)param_.getValue("parentmass_tolerance");
				Real resolution(0.1); /// @improvement depend stepsize from matchFromAligns with resolution( stepsize = resolution??)
				UInt sp = (UInt)ceil(sz/resolution); /// @improvement is floor or round more appropriate?
				std::vector< MSSpectrum<PeakType> > unmerged;
				unmerged.push_back(pairs_[*(it->second.begin())]);

				MSSpectrum<PeakType> base(pairs_[*(it->second.begin())]), rev_base(pairs_[*(it->second.begin())]);
				reverseSpectrum(rev_base);

				std::set<Size>::const_iterator others = it->second.begin();
				++others;
				for(; others != it->second.end(); ++others)
				{
					MSSpectrum<PeakType> current(pairs_[*others]);
					//find matches to base in -parentmass_tolerance:peakmass_tolerance:+parentmass_tolerance shifts
					/// @improvement if resolution = 0.1 stepsize 0.5!?
					DoubleReal best_score_base = matchFromAligns(base, current, parentmass_tolerance , 0.5);

					//find matches to rev_base in -parentmass_tolerance:peakmass_tolerance:+parentmass_tolerance shifts
					DoubleReal best_score_rev_base = matchFromAligns(rev_base, current, parentmass_tolerance , 0.5);

					//take best, possibly reverse
					if(best_score_base < best_score_rev_base)
					{
						reverseSpectrum(current);
					}
					unmerged.push_back(current);
				}

					//consensus building
					consensuses_.reserve(stars_.size());
					consensuses_.push_back(BinnedSpectrum(sz,sp,unmerged));  /// @improvement search for all spectra?
					consensuses_indices_.push_back(it->first);
			}
			/// @improvement doublettes optimization?! - only possible if spectra are not partially reversed because participating spectra migtht be equal, but orientation of base spectrum differs
		}

		void getPropagation(AASequence& template_sequence, MSSpectrum<PeakType>& current_spectrum, DoubleReal mod_pos, MSSpectrum<PeakType>& modified_sequence, DoubleReal& score)
		{

		}

	};
}
#endif //OPENMS_COMPARISON_SPECTRA_STARCLUSTERS_H
