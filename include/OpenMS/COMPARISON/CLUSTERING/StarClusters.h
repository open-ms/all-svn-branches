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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

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
			defaults_.setValue("min_hop", 0, "Defines the minimal number hops a identification can be propagated");
			defaults_.setValue("max_hop", 6, "Defines the maximal number hops a identification can be propagated (see small-world-phenomena)");
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
			return;
		}
		///

		/**
			@brief will create the network as a set of starsclusters which will most probabply overlay so its a network

			@param aligned_pairs the list of the pairs determined by an alignment (indices correspond to the original experiment)
			@param consensuses_indices are the indices that participate in the network

			star network will be represented by indices 0 to n from the participating spectra as seen in consensus_indices, not as in aligned_pairs (these indices belong to the original dataset)
		*/
		void createNetwork(std::vector< std::pair<Size,Size> >& aligned_pairs, std::map< Size, std::set<Size> >& stars) const
		{
			if(aligned_pairs.empty())
			{
				stars.clear();
				throw Exception::IllegalArgument (__FILE__, __LINE__, __PRETTY_FUNCTION__, "no edges given");
			}

			stars.clear();

			for(Size i = 0; i < aligned_pairs.size(); ++i)
			{
				stars[aligned_pairs[i].first].insert(aligned_pairs[i].second);
				stars[aligned_pairs[i].second].insert(aligned_pairs[i].first);
				/// @improvement write in metadata of consensuses[i] its pair partners indices from consensuses - this would replace stars (also in propagation), if metadata r/w speed is sufficient
				/// @improvement maybe also add mod_pos in for each pair also in metadata?!
			}
			return;
		}
		///

		/**
			@brief will score the networks annotations by explained peaks ratio

			@param experiment the MSExperiment holding the spectra to the annotations
			@param ids the ConsensusMap containing the annotations as ConsensusFeatures and corresponding identifications in associated PeptideIdentification

				*/
		void scoreNetwork(MSExperiment<PeakT>& experiment, ConsensusMap& ids) const
		{

			return;
		}
		///

		/**
			@brief will score a annotation by explained peaks ratio

			@param aligned_pairs the list of the pairs determined by an alignment (indices correspond to the original experiment)
			@param consensuses_indices are the indices that participate in the network

			star network will be represented by indices 0 to n from the participating spectra as seen in consensus_indices, not as in aligned_pairs (these indices belong to the original dataset)
		*/
		void annotationCoverage(MSSpectrum<PeakT>& spectrum, const AASequence& id, DoubleReal& peak_coverage, DoubleReal& intensity_coverage) const
		{
			DoubleReal peak_tolerance = (DoubleReal)param_.getValue("peak_tolerance");
			TheoreticalSpectrumGenerator tsg;
			RichPeakSpectrum ts;
			tsg.getSpectrum(ts, id, 1);
			ts.sortByPosition();
			DoubleReal covered_intensity(0),total_intensity(0);
			std::set<Size> covered;
			for(Size i = 0; i < spectrum.size(); ++i)
			{
				total_intensity += spectrum[i].getIntensity();
			}
			typename MSSpectrum<PeakT>::Iterator lo = spectrum.begin();
			typename MSSpectrum<PeakT>::Iterator hi = spectrum.end();
			for(Size i = 0; i < ts.size(); ++i)
			{
				DoubleReal target(ts[i].getMZ());
				lo = spectrum.MZBegin(lo,target-peak_tolerance-0.00001,spectrum.end());
				hi = spectrum.MZEnd(lo,target+peak_tolerance+0.00001,spectrum.end());
				for(; lo!=hi; ++lo)
				{
					covered.insert(lo-spectrum.begin());
					covered_intensity += lo->getIntensity();
				}
			}
			peak_coverage = (DoubleReal)covered.size()/(DoubleReal)spectrum.size();
			intensity_coverage = covered_intensity/total_intensity;
			return;
		}
		///

		/**
			@brief will clamp the annotation of a pair partner onto the current spectrum (includes possible introduction of a modification)

			@param template_sequence - the sequence to be clamped onto
			@param mod_pos - the position that has to be found otherwise the clamping would be no good anyways
			@param mod_shift - the mass shift introduced
			@param modified_sequence - the clamped sequence with introduced mass shift
			@param sequence_correspondence - which ion species the clamping is based on

			@return true if clamping of annotation was successful(find a peak corresponding to mod_pos) false otherwise

			initially looking for the peak fitting mod_pos m/z in template sequence generated y and b ions, finding is knowing the sequence position the modification took place thus a modified sequence can be constructed.
		*/
		bool getPropagation(AASequence& template_sequence, DoubleReal mod_pos, DoubleReal mod_shift, AASequence& modified_sequence, String& sequence_correspondence) const
		{
			DoubleReal peak_tolerance = (DoubleReal)param_.getValue("peak_tolerance");
			DoubleReal pm_tolerance = (DoubleReal)param_.getValue("parentmass_tolerance");

			/// @important Peptide Sequence with Modifications: the modification is written into the sequence (n- to c-terminal) after the residue that is to be modified (see below)
			if(!template_sequence.isValid())
			{
				modified_sequence = AASequence();
				mod_pos = 0;
				mod_shift = 0;
				sequence_correspondence.clear();
				return false;
			}

			if(fabs(mod_shift) <= /* 2* */pm_tolerance)
			{
				modified_sequence = template_sequence;
				mod_pos = 0;
				mod_shift = 0;
				sequence_correspondence = String("equal");
				/* debug */std::cout << " equal under threshold " << std::endl;
				return true;
			}

			if(mod_pos==0)
			{
				modified_sequence = template_sequence;
				modified_sequence.setIndesignatedModification(0, mod_shift);
				sequence_correspondence = String("equal prefixed");
				/* debug */std::cout << " equal prefixed" << std::endl;
				return true;
			}
			if(mod_pos==-1)
			{
				modified_sequence = template_sequence;
				modified_sequence.setIndesignatedModification(modified_sequence.size()-1, mod_shift);
				sequence_correspondence = String("equal suffixed");
				/* debug */std::cout << " equal suffixed" << std::endl;
				return true;
			}

			DoubleReal charge = 2;
			Size best_pos(0);

			DoubleReal best_diff = std::numeric_limits<DoubleReal>::max();
			/// @improvement seeking for match consider symmetric peaks not contained! (take symmetric info from where? the specs metadata!)
			/// @improvement seeking for match consider prm and isotope-peaks
			/// @improvement seeking for match consider unmodded if template sequence is modded
			// "Y-ions"
			for(Size i = 1; i != template_sequence.size(); ++i)
			{
				AASequence ion = template_sequence.getSuffix(i);
				DoubleReal pos = ion.getMonoWeight(Residue::YIon, (Size)charge) / charge;
				DoubleReal diff = fabs(mod_pos - pos);
				if(diff < (peak_tolerance*2 + 0.00001) and diff < best_diff)
				{
					best_diff = diff;
					sequence_correspondence = String(String("y") + String(i) + String(std::string((Size)charge, '+')));
					best_pos = template_sequence.size()-i;
				}
			}
			// "B-ions"
			for(Size i = 1; i != template_sequence.size(); ++i)
			{
				AASequence ion = template_sequence.getPrefix(i);
				DoubleReal pos = ion.getMonoWeight(Residue::BIon, (Size)charge) / charge;
				DoubleReal diff = fabs(mod_pos - pos);
				if(diff < (peak_tolerance*2 + 0.00001) and diff < best_diff)
				{
					best_diff = diff;
					sequence_correspondence = String(String("b") + String(i) + String(std::string((Size)charge, '+')));
					best_pos = i-1;
				}
			}

			if(best_diff > (peak_tolerance*2 + 0.00001))
			{
				modified_sequence = AASequence();
				mod_pos = 0;
				mod_shift = 0;
				sequence_correspondence.clear();
				/* debug */std::cout << " exceeding tolerance by " << best_diff << std::endl;
				return false;
			}

			modified_sequence = template_sequence;
			modified_sequence.setIndesignatedModification(best_pos, mod_shift);
			/// @important AASequence::getmonoweight(does not use residue::getmonoweight!) changed so above additions get accounted for as indesignated modification
			return true;
		}
		///

		/**
			@brief Method to calculate ...

			@param ids - peptide identifications coming from external and arranged in consensuses corresponding order
			@param consensuses - must correspond to the consensus_indices used in createNetwork()
			@param aligned_pairs - they represent the network edges
			@param mod_positions - must correspond to the aligned_pairs and represent the edges weight

			at the moment only original ids do not get reannotated (i.e. hits added and re-sorted), possible would be a propagation by edge score?
		*/
		void propagateNetwork(MSExperiment<PeakT>& experiment, std::vector< PeptideIdentification > peptide_ids, std::vector< ProteinIdentification > protein_ids, std::vector< std::pair<Size,Size> >& aligned_pairs, std::vector< DoubleReal > mod_positions, ConsensusMap& ids) const
		{
			Size hops = (Size)param_.getValue("max_hop");
			if(peptide_ids.empty())
			{
				throw Exception::IllegalArgument (__FILE__, __LINE__, __PRETTY_FUNCTION__, "no identifications to propagate");
			}

			Real peak_tolerance = (Real)param_.getValue("peak_tolerance");
			DoubleReal pm_tolerance = (DoubleReal)param_.getValue("parentmass_tolerance");

			//aa_masses[i] holds internal mass of aa i
			//make unique in resolution tolerance (e.g. cause of Q & K)
			const ResidueDB* res_db = ResidueDB::getInstance();
			std::vector<Real> aa_masses;
			for(ResidueDB::ResidueConstIterator it = res_db->beginResidue(); it != res_db->endResidue(); ++it)
			{
				DoubleReal w = (*it)->getMonoWeight() - (*it)->getInternalToFullMonoWeight();
				/// @important: this makes only singly charged ion peaks to potential jump matches and therefor potential alignment-partners (so a test is only reasonable with singly charged ion peaks in da spec)
				if(w>0)
				{
					aa_masses.push_back(w);
				}
			}
			//sort aa_masses
			std::sort(aa_masses.begin(),aa_masses.end());
			std::vector<Real>::iterator end_new = std::unique(aa_masses.begin(),aa_masses.end(), EqualInTolerance<Real>(peak_tolerance));
			aa_masses.resize( end_new - aa_masses.begin() );

			for(Size i = 0; i < experiment.size(); ++i)
			{
				ConsensusFeature cf;
				cf.setMZ(experiment[i].getPrecursors().front().getMZ());
				cf.setRT(experiment[i].getRT());
				ids.push_back(cf);
			}

			IDMapper idm;
			idm.annotate(ids,peptide_ids,protein_ids);

			//~ set initial anntoating indices, network id s and original annotation tag
			std::vector<Size> original_annotations;
			Size net_ids(0);
			for(ConsensusMap::Iterator ids_it = ids.begin(); ids_it != ids.end(); ++ids_it)
			{
				if(ids_it->getPeptideIdentifications().size()<1)
				{
					ids_it->getPeptideIdentifications().resize(1); // make sure we have the one PeptideIdentification we need, n.b. not to be confused with the subsequent hits!
					ids_it->getPeptideIdentifications().front().setIdentifier("from SpectralNetworking");
				}
				else
				{
					/* debug */DoubleReal hits_no(0), pi_no(0);
					std::vector<PeptideHit> new_hits;
					/// @todo iterator varianten sind wohl schneller!
					for(Size i = 0; i < ids_it->getPeptideIdentifications().size(); ++i)
					{
						/* debug */++pi_no;
						/* debug */hits_no += ids_it->getPeptideIdentifications()[i].getHits().size();
						//~ pi.getHits().insert(pi.getHits().end()-1,ids_it->getPeptideIdentifications()[i].getHits().begin(),ids_it->getPeptideIdentifications()[i].getHits().end());
						//~ for(Size j = 0; j < ids_it->getPeptideIdentifications()[i].getHits().size(); ++j)
						for(std::vector<PeptideHit>::const_iterator it_hits = ids_it->getPeptideIdentifications()[i].getHits().begin(); it_hits != ids_it->getPeptideIdentifications()[i].getHits().end(); ++it_hits)
						{
							new_hits.push_back(*it_hits);
							DoubleReal pc(0),ic(0);
							annotationCoverage(experiment[ids_it-ids.begin()], new_hits.back().getSequence(), pc, ic);
							new_hits.back().setMetaValue("IdentificationIdentifier", String(ids_it->getPeptideIdentifications()[i].getIdentifier()));
							new_hits.back().setMetaValue("annotated with peak coverage", pc);
							new_hits.back().setMetaValue("annotated with intensity coverage", ic);
							//~ /// @improvement factorized linear combination of pc and ic for advanced parameter - @see below!
							new_hits.back().setScore(((pc+ic)/2.0));
						}
					}
					/* debug */std::cout << " average "<< (hits_no/pi_no) << " hits per " << pi_no << " PeptideIdentifications." << std::endl;
					ids_it->getPeptideIdentifications().resize(1);
					ids_it->getPeptideIdentifications().front().setHits(new_hits);
					ids_it->getPeptideIdentifications().front().setIdentifier("from SpectralNetworking");
				}

				if(!ids_it->getPeptideIdentifications().front().getHits().empty())
				{
					ids_it->getPeptideIdentifications().front().sort();
					//~ ids_it->getPeptideIdentifications().front().assignRanks();

					ids_it->setMetaValue("network id", net_ids);
					++net_ids;
					ids_it->setMetaValue("DBid", DataValue::EMPTY_VALUE);
					original_annotations.push_back(ids_it - ids.begin());
				}

			}
			/* debug */std::cout << " consensusfeature annotated: " << net_ids-1 << std::endl; Size annots(0);

			std::map< Size, std::set<Size> > stars;
			createNetwork(aligned_pairs, stars);

			//~ set potential targets
			std::list<Size> potential_propagation_targets;
			for(Size i = 0; i < original_annotations.size(); ++i)
			{
				for(std::set<Size>::const_iterator it_star = stars[original_annotations[i]].begin(); it_star != stars[original_annotations[i]].end(); ++it_star)
				{
					potential_propagation_targets.insert(potential_propagation_targets.end(),*it_star);
				}
			}
			potential_propagation_targets.unique();
			potential_propagation_targets.sort();

			//~ original ids dont need re-annotation!
			for(Size i = 0; i < original_annotations.size(); ++i)
			{
				potential_propagation_targets.remove(original_annotations[i]);
			}

			std::map< Size, std::set<Size> > conflict_map;
			std::vector<Size> all_annotations(original_annotations);

			//~ propagation in |hops| steps:
			for(Size i = 0; i < hops; ++i)
			{
				std::vector<Size> latest_annotated;
				/// @improvement reserve potential_propagation_targets.size() space

				for(std::list<Size>::const_iterator it_potential = potential_propagation_targets.begin(); it_potential != potential_propagation_targets.end(); ++it_potential)
				{
					const MSSpectrum<PeakType>& current_spectrum = (experiment[*it_potential]);
					ConsensusFeature& current_annotations = (ids[*it_potential]);
					/* debug */std::cout << " at : " << *it_potential;
					if(current_annotations.getPeptideIdentifications().front().getHits().empty())
					{
						//~ current spectrum not annotated yet
						current_annotations.setMetaValue("SNid", DataValue::EMPTY_VALUE);
					}

					for(std::set<Size>::const_iterator it_neighbors = stars[*it_potential].begin(); it_neighbors != stars[*it_potential].end(); ++it_neighbors)
					{
						std::vector<Size>::iterator annotating_neighbor = std::find(all_annotations.begin(),all_annotations.end(),*it_neighbors);
						if(!ids[*it_neighbors].getPeptideIdentifications().front().getHits().empty())
						{
							//~ at this point we have a edge i,j from aligned_pairs:  i=*it_neighbors , j=*it_potential and if i>j: swap(i,j)

							DoubleReal pm_s1 = current_spectrum.getPrecursors().front().getMZ();
							int c_s1 = current_spectrum.getPrecursors().front().getCharge();
							DoubleReal current_pm((pm_s1*c_s1 + (c_s1-1)*Constants::PROTON_MASS_U));
							DoubleReal pm_s2 = experiment[*it_neighbors].getPrecursors().front().getMZ();
							int c_s2 = experiment[*it_neighbors].getPrecursors().front().getCharge();
							/// @attention if the precursor charge is unknown, i.e. 0 best guess is its doubly charged
							(c_s1==0)?c_s1=2:c_s1=c_s1;
							(c_s2==0)?c_s2=2:c_s2=c_s2;
							/// @attention singly charged mass difference!
							DoubleReal pm_diff = (pm_s2*c_s2 + (c_s2-1)*Constants::PROTON_MASS_U)-current_pm;

							std::pair<Size,Size> current_edge(*it_potential,*it_neighbors);

							if(fabs(pm_diff)<pm_tolerance)
							{
								if(current_edge.first>current_edge.second)
								{
									std::swap(current_edge.first,current_edge.second);
								}
							}
							else if(pm_diff<0)
							{
								std::swap(current_edge.first,current_edge.second);
							}

							DoubleReal mod_pos(0.0);
							std::vector< std::pair<Size,Size> >::iterator which_pair_it = std::find(aligned_pairs.begin(),aligned_pairs.end(),current_edge);
							if(which_pair_it!=aligned_pairs.end())
							{
								mod_pos = mod_positions[which_pair_it-aligned_pairs.begin()];
								/// @attention simplify for getPropagation mod_pos so it has no need to know current_spectrum or current_pm neg. meaning the mod is suffixed
								/// @attention mod_positions and aligned_pairs have to be of same size else a segfault might occur (if user of this function does not read the docu)
								if(pm_diff>0 and mod_pos != -1)
								{
									/// @attention pm_diff<0 so we are propagating from light to heavy (modpos always in the first spectrum of a edge (the spectrum with lower mz) - edge was swapped before!, so adapt)
									mod_pos += fabs(pm_diff);
									//~ like this we get the modpos in current spectrum
								}
							}
							else
							{
								String element(String(current_edge.first) + String("-") + String(current_edge.second));
								throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, element);
							}

							/* debug */std::cout << " attempting annotation from spec" << (*it_neighbors) << " on mod_pos " << mod_pos << " and mod_shift " << -pm_diff;

							for(Size k = 0; k < ids[*it_neighbors].getPeptideIdentifications().front().getHits().size() and k < hops; ++k)
							{
								AASequence template_sequence(ids[*it_neighbors].getPeptideIdentifications().front().getHits()[k].getSequence()), modified_sequence; /// @improvement find out if annotate with _all_ hits from annotating_neighbor is better

								String mod_found_with_ions; //which ion species helped to find mod
								/// @attention sideeffect of function fills empty arguments if evaluated true
								/// @attention pm_diff > 0 so we are propagating from heavy to light, so mod_shift must be negative, pm_diff<0 other way'round and mod_shift must be negative!
								if(getPropagation(template_sequence, mod_pos, -pm_diff, modified_sequence, mod_found_with_ions))
								{
									/* debug */ ++annots;
									int charge = ids[*it_neighbors].getPeptideIdentifications().front().getHits()[k].getCharge();
									PeptideHit hit(0,0,charge,modified_sequence); /// @attention score and rank have to be set after intensity and peakmatch coverage evaluation (own step)

									hit.setMetaValue("IdentificationIdentifier", ids[*it_neighbors].getPeptideIdentifications().front().getHits()[k].getMetaValue("IdentificationIdentifier"));

									FeatureHandle annotator;
									//~ define edge
									if(fabs(pm_diff) <= pm_tolerance)
									{
										hit.setMetaValue("annotation relation", "equal");
										//~ annotator.setMetaValue("annotation relation", "equal");
									}
									else if(fabs(pm_diff) <= aa_masses.back()+pm_tolerance)
									{
										for(Size i = 0; i < aa_masses.size(); ++i)
										{
											if(fabs(fabs(pm_diff)-aa_masses[i]) <= pm_tolerance)
											{
												hit.setMetaValue("annotation relation", "ladder");
												//~ annotator.setMetaValue("annotation relation", "ladder");
												break;
											}
										}
									}
									else
									{
										hit.setMetaValue("annotation relation", "modification");
										//~ annotator.setMetaValue("annotation relation", "modification");
									}

									//~ make annotating neighbor visible
									if(ids[*it_neighbors].getFeatures().empty())
									{
										FeatureHandle neighbor_feat;
										neighbor_feat.setElementIndex(*it_neighbors);
										neighbor_feat.setMZ(experiment[*it_neighbors].getPrecursors().front().getMZ());
										neighbor_feat.setRT(experiment[*it_neighbors].getRT());
										neighbor_feat.setIntensity(0);
										ids[*it_neighbors].insert(neighbor_feat);
									}

									//~ color, network id and metadata management
									hit.setMetaValue("annotator index", *it_neighbors);
									hit.setMetaValue("annotator MZ", experiment[*it_neighbors].getPrecursors().front().getMZ());
									hit.setMetaValue("annotator RT", experiment[*it_neighbors].getRT());
									hit.setMetaValue("annotated in hop no.", i+1);
									hit.setMetaValue("annotated with ionspecies", mod_found_with_ions);
									DoubleReal pc(0),ic(0);
									annotationCoverage(experiment[*it_potential], modified_sequence, pc, ic);
									hit.setMetaValue("annotated with peak coverage", pc);
									hit.setMetaValue("annotated with intensity coverage", ic);
									/// @improvement factorized linear combination of pc and ic for advanced parameter
									hit.setScore(((pc+ic)/2.0));

									annotator.setElementIndex(*it_neighbors);
									annotator.setMZ(experiment[*it_neighbors].getPrecursors().front().getMZ());
									annotator.setRT(experiment[*it_neighbors].getRT());
									//~ annotator.setIntensity(((pc+ic)/2.0));

									try
									{
										current_annotations.insert(annotator);
									}
									catch(Exception::InvalidValue& e){
										/* debug */std::cout << " again ";
									}
									current_annotations.getPeptideIdentifications().front().insertHit(hit);

									if(!current_annotations.metaValueExists("network id"))
									{
										current_annotations.setMetaValue("network id", (Size)ids[*it_neighbors].getMetaValue("network id"));
									}
									else
									{
										Size old_id = (Size)current_annotations.getMetaValue("network id");
										Size new_id = (Size)ids[*it_neighbors].getMetaValue("network id");
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
								/* debug */else
								{

								}
								/* debug */std::cout << std::endl;
							}
						}
					}
					if(!current_annotations.getPeptideIdentifications().front().getHits().empty())
					{
						latest_annotated.push_back(*it_potential);
						all_annotations.push_back(*it_potential);
						current_annotations.getPeptideIdentifications().front().sort();
						/// @improvement make featurehandles  DataFilter-able (order same as hits?!)
					}
				}

				//~ update potential_propagation_targets
				/// @improvement for dense networks this! will! suck! ... bigtime - so at least dont do in last step cos not needed
				if(i!=(hops-1))
				{
					potential_propagation_targets.clear();
					for(Size k = 0; k < latest_annotated.size(); ++k)
					{
						for(std::set<Size>::const_iterator it_star = stars[latest_annotated[k]].begin(); it_star != stars[latest_annotated[k]].end(); ++it_star)
						{
							potential_propagation_targets.insert(potential_propagation_targets.begin(),*it_star);
						}
					}
					latest_annotated.clear();
					potential_propagation_targets.unique();
					potential_propagation_targets.sort();
					std::sort(all_annotations.begin(),all_annotations.end());
					/// @improvement find out whats better, remove all annotated or only the originals
					//~ erase all annotated/originally annotated indices from potential_propagation_targets.
					//~ for(Size j = 0; j < original_annotations.size(); ++j)
					for(Size j = 0; j < all_annotations.size(); ++j)
					{
						//~ potential_propagation_targets.remove(original_annotations[j]);
						potential_propagation_targets.remove(all_annotations[j]);
					}
				}
			}
			/// @improvement keep track of used edges and force annotation of 'good' edges (scores not yet available in here :( )

			//~ concluding conflict resolving in network numbering
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
				++conflict_map_it;
			}

			std::map< Size, Size > conter_conflict_map;
			for(std::map< Size, std::set<Size> >::iterator conflict_map_it = conflict_map.begin(); conflict_map_it != conflict_map.end(); ++conflict_map_it)
			{
				for(std::set<Size>::iterator conflicts_it = conflict_map_it->second.begin(); conflicts_it != conflict_map_it->second.end(); ++conflicts_it)
				{
					conter_conflict_map[*conflicts_it] = conflict_map_it->first;
				}
			}
			for(ConsensusMap::Iterator ids_it = ids.begin(); ids_it != ids.end(); ++ids_it)
			{
				if(ids_it->metaValueExists("network id"))
				{
					ids_it->setQuality( ids_it->getPeptideIdentifications().front().getHits().front().getScore() );

					std::map< Size, Size >::iterator cur_it = conter_conflict_map.find((Size)ids_it->getMetaValue("network id"));
					if(cur_it != conter_conflict_map.end())
					{
						ids_it->setMetaValue("network id", cur_it->second);
					}

					/// @todo this makes the annotations stay and not viewable at the moment in toppview - adapt toppview to indesignated modifications
					ids_it->getPeptideIdentifications().front().setIdentifier(ids_it->getPeptideIdentifications().front().getHits().front().getMetaValue("IdentificationIdentifier"));
				}
			}

		}
		///

	};
}
#endif //OPENMS_COMPARISON_SPECTRA_STARCLUSTERS_H
