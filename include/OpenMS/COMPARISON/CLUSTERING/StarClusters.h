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
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
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

	public:

		// @name Constructors and Destructors
		// @{

		/// default constructor
		StarClusters()
				: DefaultParamHandler("StarClusters")
		{
			defaults_.setValue("peak_tolerance", 0.3, "Defines the absolut (in Da) peak tolerance");
			defaults_.setValue("parentmass_tolerance", 3.0, "Defines the absolut (in Da) parent mass tolerance");
			defaults_.setValue("min_shift", 0.0, "Defines the minimal absolut (in Da) shift between the two spectra");
			defaults_.setValue("max_shift", 150.0, "Defines the maximal absolut (in Da) shift between the two spectra");
			defaults_.setValue("min_dist", 57.0214637230, "Defines the minimal distance (in Da) between the two peaks of one sort that may be connected in a sparse matching");
			defaults_.setValue("correlation_scoring", "intensity", "If intensity, correlation scores on basis of matched intensity values, if matchnumber correlation scores solely on basis of the number of matches");
			defaults_.setValidStrings("correlation_scoring", StringList::create("intensity,matchnumber"));
			defaultsToParam_();
		}

		/// copy constructor
		StarClusters(const StarClusters& source)
			: DefaultParamHandler(source)
		{
		}
		///

		/// destructor
		~StarClusters()
		{
		}
		///

		/// assignment operator
		StarClusters& operator = (const StarClusters& source)
		{
			if (this != &source)
			{
				DefaultParamHandler::operator = (source);
			}
			return *this;
		}
		///

		// @}

		/**
			@brief This reverses a spectrum

			@param spec the spectrum to reverse
			@param is_prm if the spectrum is a prm spectrum

			Each peak p is repositioned at pm - p->getMZ(), where pm is the precursor mass of spec (-1 if prm); spec is sortet by position afterwards
		*/
		void reverseSpectrum(MSSpectrum<PeakType>& spec, bool is_prm=false)
		{
			DoubleReal pm = spec.getPrecursors().front().getMZ(); /// @improvement throw error if no precursor
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
		///

		/**
			@brief will clamp the annotation of a pair partner onto the current spectrum (includes possible introduction of a modification)

			@param ...
			@param ...

			...
		*/
		void getPropagation(AASequence& template_sequence, DoubleReal mod_pos, DoubleReal mod_shift, AASequence& modified_sequence, String& sequence_correspondence)
		{
			DoubleReal area_peak_tolerance = (double)param_.getValue("peak_tolerance")*2 + 0.00001;

			/// @important Peptide Sequence with Modifications: the modification is written into the sequence (n- to c-terminal) after the residue that is to be modified (see below)
			if (!template_sequence.isValid())
			{
				modified_sequence = AASequence();
				mod_pos = 0;
				mod_shift = 0;
				sequence_correspondence.clear();
				return;
			}

			DoubleReal charge = 2;
			DoubleReal best_diff = std::numeric_limits<DoubleReal>::max();
			Size best_pos;


			/// @improvement seeking for match consider prm and isotope-peaks
			// "Y-ions"
			for (Size i = 1; i != template_sequence.size(); ++i)
			{
				AASequence ion = template_sequence.getSuffix(i);
				DoubleReal pos = ion.getMonoWeight(Residue::YIon, (Size)charge) / charge;
				DoubleReal diff = abs(mod_pos - pos);
				if(diff < area_peak_tolerance and diff < best_diff)
				{
					best_diff = diff;
					sequence_correspondence = "y"+String(i) + String((Size)charge, "+ shift");
					best_pos = template_sequence.size()-i;
				}
			}
			// "B-ions"
			for (Size i = 1; i != template_sequence.size(); ++i)
			{
				AASequence ion = template_sequence.getPrefix(i);
				DoubleReal pos = ion.getMonoWeight(Residue::BIon, (Size)charge) / charge;
				DoubleReal diff = abs(mod_pos - pos);
				if(diff < area_peak_tolerance and diff < best_diff)
				{
					best_diff = diff;
					sequence_correspondence = "b"+String(i) + String((Size)charge, "+ shift");
					best_pos = i-1;
				}
			}

			modified_sequence = template_sequence;
			modified_sequence.setIndesignatedModification(best_pos, mod_shift);
			/// @important AASequence::getmonoweight(does not use residue::getmonoweight!) changed so above additions get accounted for as indesignated modification
		}
		///

		/**
			@brief Method to calculate ...

			@param ids  - peptide identifications coming from external and arranged in consensuses corresponding order
			@param consensuses  - must correspond to the consensus_indices used in createNetwork()

			...
			@see ...
		*/
		void spanNetwork(MSExperiment<PeakType>& consensuses, ConsensusMap& id_consensuses, std::vector<Size>& indices_id_consensuses, std::vector< std::pair<Size,Size> >& aligned_pairs, std::vector<DoubleReal>& mod_positions, Size hops)
		{

			if(indices_id_consensuses.empty() or indices_id_consensuses.size()!=id_consensuses.size())
			{
				/// @todo throw error
				return;
			}

			//~ set initial network id s
			for(ConsensusMap::Iterator id_consensuses_it = id_consensuses.begin(); id_consensuses_it != id_consensuses.end(); ++id_consensuses_it)
			{
				id_consensuses_it->setMetaValue("network id", id_consensuses_it - id_consensuses.begin());
			}

			std::map< Size, std::set<Size> > stars;
			createNetwork(aligned_pairs, consensuses, stars);

			std::list<Size> potential_propagation_targets;
			//~ prime annotation-step with ids
			for(Size i = 0; i < indices_id_consensuses.size(); ++i)
			{
				for(std::set<Size>::const_iterator it_star = stars[indices_id_consensuses[i]].begin(); it_star != stars[indices_id_consensuses[i]].end(); ++it_star)
				{
					potential_propagation_targets.insert(potential_propagation_targets.end(),*it_star);
				}
			}
			potential_propagation_targets.unique();
			potential_propagation_targets.sort();
			//~ yet annotated dont need re-annotation!
			for(Size i = 0; i < indices_id_consensuses.size(); ++i)
			{
				potential_propagation_targets.remove(i);
			}

			std::map< Size, std::set<Size> > conflict_map;

			//~ propagation in |hops| steps:
			for(Size i = 0; i < hops; ++i)
			{
				std::vector<Size> newly_annotated; /// @improvement reserve potential_propagation_targets.size() space

				for(std::list<Size>::const_iterator it_potential = potential_propagation_targets.begin(); it_potential != potential_propagation_targets.end(); ++it_potential)
				{
					const MSSpectrum<PeakType>& current_spectrum = (consensuses[*it_potential]);
					ConsensusFeature current_annotations;
					current_annotations.setMZ(consensuses[*it_potential].getMZ());
					current_annotations.setRT(consensuses[*it_potential].getRT());
					current_annotations.getPeptideIdentifications().resize(1); // make sure we have the one PeptideIdentification we need
					for(std::set<Size>::const_iterator it_neighbors = stars[*it_potential].begin(); it_neighbors != stars[*it_potential].end(); ++it_neighbors)
					{
						std::vector<Size>::iterator annotating_neighbor = std::find(indices_id_consensuses.begin(),indices_id_consensuses.end(),*it_neighbors);

						if(annotating_neighbor!=indices_id_consensuses.end())
						{
							//~ at this point we have a pair i,j from aligned_pairs:  i=*it_neighbors , j=*it_potential and if i>j: swap(i,j)
							Size annotating_neighbor_index = annotating_neighbor - indices_id_consensuses.begin(); // corresponding ConsensusFeature in ConsensusMap has the same index
							AASequence template_sequence(id_consensuses[annotating_neighbor_index].getPeptideIdentifications().front().getHits().front().getSequence()), modified_sequence; /// @improvement find out if annotate with _all_ hits from annotating_neighbor is better

							DoubleReal pm_diff, mod_pos;
							/// @todo consider special mod_pos, see 543 in AntisymetricAlignment.h
							std::vector< std::pair<Size,Size> >::iterator which_pair_it = std::find(aligned_pairs.begin(),aligned_pairs.end(),std::make_pair((Size)consensuses[*it_potential].getNativeID().toInt(),(Size)consensuses[*it_neighbors].getNativeID().toInt()));
							if(which_pair_it!=aligned_pairs.end())
							{
								//~ modpos always in the first spectrum of a pair, here (*it_potential) = current_spectrum
								//~ but we annotate from (*it_neighbors), so the corresponding mop_pos in (*it_neighbors) has to be |pm_diff| greater and pm_diff will be 0 or negative
								pm_diff = current_spectrum.getPrecursors().front().getMZ() - consensuses[*it_neighbors].getPrecursors().front().getMZ();
								mod_pos = mod_positions[which_pair_it-aligned_pairs.begin()] + abs(pm_diff); /// @attention mod_positions and aligned_pairs have to be of same size else a segfault might occur (if user of this function does not read the docu)
							}
							else
							{
								which_pair_it = std::find(aligned_pairs.begin(),aligned_pairs.end(),std::make_pair((Size)consensuses[*it_neighbors].getNativeID().toInt(),(Size)consensuses[*it_potential].getNativeID().toInt()));
								if(which_pair_it!=aligned_pairs.end())
								{
									//modpos always in the first spectrum of a pair, here (*it_neighbors)
									pm_diff = current_spectrum.getPrecursors().front().getMZ() - consensuses[*it_neighbors].getPrecursors().front().getMZ(); // will be 0 or positive
									mod_pos = mod_positions[which_pair_it-aligned_pairs.begin()];
								}
								//~ else
									/// @todo throw(Exception "youre screwed!!");
							}
							String mod_found_with_ions; //which ions helped to find mod
							getPropagation(template_sequence, mod_pos, pm_diff, modified_sequence,mod_found_with_ions);

							int charge = id_consensuses[annotating_neighbor_index].getPeptideIdentifications().front().getHits().front().getCharge();
							PeptideHit hit(0,0,charge,modified_sequence); /// @attention score and rank have to be set after intensity and peakmatch coverage evaluation (own step)
							hit.setMetaValue("annotator index", *it_neighbors);
							hit.setMetaValue("annotator MZ", consensuses[*it_neighbors].getMZ());
							hit.setMetaValue("annotator RT", consensuses[*it_neighbors].getRT());
							hit.setMetaValue("annotated in hop no.", i+1);
							hit.setMetaValue("annotated with what ions", mod_found_with_ions);
							FeatureHandle annotator;
							annotator.setElementIndex(*it_neighbors);
							annotator.setMZ(consensuses[*it_neighbors].getMZ());
							annotator.setRT(consensuses[*it_neighbors].getRT());
							current_annotations.insert(annotator);
							current_annotations.getPeptideIdentifications().front().insertHit(hit);

							if(!current_annotations.metaValueExists("network id"))
							{
								current_annotations.setMetaValue("network id", id_consensuses[annotating_neighbor_index].getMetaValue("network id"));
							}
							else
							{
								Size old_id = current_annotations.getMetaValue("network id");
								Size new_id = id_consensuses[annotating_neighbor_index].getMetaValue("network id");
								if(new_id < old_id)
								{
									current_annotations.setMetaValue("network id", new_id);
									conflict_map[old_id].insert(new_id);
								}
								else
								{
									if(new_id != old_id)
									{
										conflict_map[new_id].insert(old_id);
									}
								}
							}

						}
					}

					if(!current_annotations.getPeptideIdentifications().front().getHits().empty())
					{
						current_annotations.getPeptideIdentifications().front().sort();
						id_consensuses.push_back(current_annotations);
						indices_id_consensuses.push_back(*it_potential);
						newly_annotated.push_back(*it_potential);
					}
				}

				//~ update potential_propagation_targets
				/// @improvement for dense networks this! will! suck! ... bigtime - so at least dont do in last step cos not needed
				if(i!=(hops-1))
				{
					for(Size k = 0; k < newly_annotated.size(); ++k)
					{
						for(std::set<Size>::const_iterator it_star = stars[newly_annotated[k]].begin(); it_star != stars[newly_annotated[k]].end(); ++it_star)
						{
							potential_propagation_targets.insert(potential_propagation_targets.begin(),*it_star);
						}
					}
					newly_annotated.clear();
					potential_propagation_targets.unique();
					potential_propagation_targets.sort();
					//~ erase all annotated indices from potential_propagation_targets.
					for(Size j = 0; j < indices_id_consensuses.size(); ++j)
					{
						potential_propagation_targets.remove(j);
					}
				}
			}

			/// @improvements find if resolve conflicts in network id's and colorize correspondingly slows everything down
			std::map< Size, std::set<Size> >::iterator conflict_map_it = conflict_map.begin();
			while(conflict_map_it != conflict_map.end())
			{
				for(std::set<Size>::iterator conflicts_it = conflict_map_it->second.begin(); conflicts_it != conflict_map_it->second.end(); ++conflicts_it)
				{
					std::map< Size, std::set<Size> >::iterator others = conflict_map.find(*conflicts_it);
					if(others != conflict_map.end())
					{
						conflict_map_it->second.insert(others->second.begin(), others->second.end());
						conflict_map.erase(others);
					}
				}
			}
			std::map< Size, Size > conter_conflict_map;
			for(std::map< Size, std::set<Size> >::iterator conflict_map_it = conflict_map.begin(); conflict_map_it != conflict_map.end(); ++conflict_map_it)
			{
				for(std::set<Size>::iterator conflicts_it = conflict_map_it->second.begin(); conflicts_it != conflict_map_it->second.end(); ++conflicts_it)
				{
					conter_conflict_map[*conflicts_it] = conflict_map_it->first;
				}
			}
			std::vector<String> colors;
			colors.push_back("#00FFFF");
			colors.push_back("#000000");
			colors.push_back("#0000FF");
			colors.push_back("#FF00FF");
			colors.push_back("#008000");
			colors.push_back("#808080");
			colors.push_back("#00FF00");
			colors.push_back("#800000");
			colors.push_back("#000080");
			colors.push_back("#808000");
			colors.push_back("#800080");
			colors.push_back("#FF0000");
			colors.push_back("#C0C0C0");
			colors.push_back("#008080");
			colors.push_back("#FFFF00");
			for(ConsensusMap::Iterator id_consensuses_it = id_consensuses.begin(); id_consensuses_it != id_consensuses.end(); ++id_consensuses_it)
			{
				std::map< Size, Size >::iterator cur_it = conter_conflict_map.find(id_consensuses_it->getMetaValue("network id"));
				if(cur_it != conter_conflict_map.end())
				{
					id_consensuses_it->setMetaValue("network id", cur_it->second);
					id_consensuses_it->setMetaValue("color", colors[(Size)id_consensuses_it->getMetaValue("network id")%colors.size()]);
				}
				else
				{
					id_consensuses_it->setMetaValue("color", colors[(Size)id_consensuses_it->getMetaValue("network id")%colors.size()]);
				}
			}
			id_consensuses.sortByPosition(); //invalidates consensus indices but each feature has its spectrum in consensus covered
		}
		///

		/**
			@brief will create the network as a set of starsclusters which will most probabply overlay so its a network

			@param aligned_pairs the list of the pairs determined by an alignment (indices correspond to the original experiment)
			@param consensuses_indices are the indices that participate in the network

			star network will be represented by indices 0 to n from the participating spectra as seen in consensus_indices, not as in aligned_pairs (these indices belong to the original dataset)
		*/
		void createNetwork(std::vector< std::pair<Size,Size> >& aligned_pairs, MSExperiment<PeakType>& consensuses, std::map< Size, std::set<Size> >& stars)
		{
			if(aligned_pairs.empty())
			{
				/// @todo throw error
				stars.clear();
				return;
			}

			stars.clear();

			std::map<Size,Size> consensus_indices_map;

			for(Size i = 0; i < consensuses.size(); ++i)
			{
				consensus_indices_map[(Size)(consensuses[i].getNativeID())] = i;
			}

			for(Size i = 0; i < aligned_pairs.size(); ++i)
			{
				stars[consensus_indices_map[aligned_pairs[i].first]].insert(consensus_indices_map[aligned_pairs[i].second]);
				stars[consensus_indices_map[aligned_pairs[i].second]].insert(consensus_indices_map[aligned_pairs[i].first]);
				/// @improvement write in metadata of consensuses[i] its pair partners indices from consensuses - this would replace stars (also in propagation), if metadata r/w speed is sufficient
				/// @improvement maybe also add mod_pos in for each pair also in metadata?!
			}
		}
		///

	};
}
#endif //OPENMS_COMPARISON_SPECTRA_STARCLUSTERS_H
