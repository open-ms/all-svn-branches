// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>

#include<OpenMS/SYSTEM/StopWatch.h>

using namespace std;

namespace OpenMS
{

BaseSweepSeeder::BaseSweepSeeder()
        : BaseSeeder(),
				is_initialized_(false),
				mass_tolerance_alignment_(0),
				scans_to_sumup_(0)
{
		// number of scans used for alignment
		defaults_.setValue("scans_to_sumup",5,"number of scans that are combined in order to improve the signal-to-noise level.");
		// mass tolerance during scan alignment
		defaults_.setValue("mass_tolerance_alignment", 0.1,"mass tolerance for comination of peaks from different scans.");
		
			// minimum number of scans per isotopic cluster
		defaults_.setValue("min_number_scans",5,"min. number of scans");
			// maximum number of scans per isotopic cluster
		defaults_.setValue("max_number_scans",700,"max. number of scans");
		// minimum number of peaks per cluster
		defaults_.setValue("min_number_peaks",20,"minimum number of peaks");		
		// minimum number of peaks per cluster
		defaults_.setValue("min_number_points_per_scan",8,"minimum number of data points per scan");		
		
		// mass tolerance for peak cluster construction
		defaults_.setValue("mz_tolerance_cluster",1.2,"m/z tolerance for looking up a signal in following scans.");
		// rt tolerance for cluster construction (given in number of scans)
		defaults_.setValue("rt_tolerance_cluster",2,"rt tolerance for cluster construction (given in number of scans)");		
		
		// mass tolerance for peak cluster construction
		defaults_.setValue("mz_tolerance_merging",2.0,"m/z tolerance for the merging of overlapping cluster.");
		
		// minimum false discovery rate for a significant cluster
		defaults_.setValue("fdr_alpha",5.0,"minimum false discovery rate for a significant isotopic pattern");
}

BaseSweepSeeder::BaseSweepSeeder(const BaseSweepSeeder& source) : BaseSeeder(source) {}

BaseSweepSeeder::~BaseSweepSeeder() {}

BaseSweepSeeder& BaseSweepSeeder::operator = (const BaseSweepSeeder& source)
{
    if (&source == this)
        return *this;

    BaseSeeder::operator = (source);

    return *this;
}

FeaFiModule::ChargedIndexSet BaseSweepSeeder::nextSeed() throw (NoSuccessor)
{
		if (!is_initialized_)
		{
			sweep_();		// sweep across map and scan for pattern 
			curr_region_  = iso_map_.begin();
			is_initialized_ = true;
		}
		
		if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 )
		{
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, make_pair(0,0));
		}
		
		cout << "Retrieving next region with charge " << (*curr_region_).second.peaks_.charge_; 
		cout << " and size " << (*curr_region_).second.peaks_.size() << endl;
			
		return (curr_region_++)->second.peaks_;
}

void BaseSweepSeeder::updateMembers_()
{
	// params for scan alignment
	mass_tolerance_alignment_ = param_.getValue("mass_tolerance_alignment");
	scans_to_sumup_           			= param_.getValue("scans_to_sumup");

	// sweepline params
	mz_tolerance_cluster_   	= param_.getValue("mz_tolerance_cluster");
	rt_tolerance_cluster_     	= (UInt) param_.getValue("rt_tolerance_cluster");
	
	// peak cluster merging (m/z for overlapping peak cluster that are merged)
	max_mz_dist_merging_   	= param_.getValue("mz_tolerance_merging");
}

void BaseSweepSeeder::sweep_()
{
		// progress logger
		traits_->startProgress(0, traits_->getData().size() , "FeatureFinder");
		
		// copy current scan. 
		// This is necessary as the peak intensities have to be modified in the sumUp_ method
		for (UInt currscan_index = 0; currscan_index < traits_->getData().size(); ++currscan_index)
		{		
			traits_->setProgress(currscan_index);
			
			SpectrumType current_scan  = traits_->getData()[ currscan_index ];		
			
			cout << "---------------------------------------------------------------------------" << endl;
			cout << "RT: " << current_scan.getRT() << " ";
			cout << " spectrum " << (currscan_index + 1) << " of " << traits_->getData().size() << endl;
						
			#ifdef DEBUG_FEATUREFINDER
			// write debug output
			String fname = String("scan_") + current_scan.getRT();;
			ofstream out( fname.c_str() );
			for(UInt k = 0; k<current_scan.size();++k)
			{
				out << current_scan[k].getMZ() << " " << current_scan[k].getIntensity() << endl;
			}
			out.close();
			#endif
			
  		sumUp_(current_scan,currscan_index);
						
			#ifdef DEBUG_FEATUREFINDER
			// write debug output
			fname = String("scan_aligned_") + current_scan.getRT();;
			out.open( fname.c_str() );
			for(UInt k = 0; k<current_scan.size();++k)
			{
				out << current_scan[k].getMZ() << " " << current_scan[k].getIntensity() << endl;
			}
			out.close();
			#endif
						
			// detect isotopic pattern...	
			ScoredMZVector iso_curr_scan = detectIsotopicPattern_(current_scan );
			
			for (ScoredMZVector::const_iterator citer = iso_curr_scan.begin();
						citer != iso_curr_scan.end();
						++citer)
			{
 				traits_->getPeakFlag( make_pair( currscan_index, citer->first ) ) = FeaFiTraits::USED;			
			}
			
			cout << iso_curr_scan.size() << " feature candidates detected." << endl;
			
			// for each m/z position with score: 
			// => check for cluster at similar m/z in previous scans
			// => if (matching cluster found) extend
			for (ScoredMZVector::const_iterator citer = iso_curr_scan.begin();
						citer != iso_curr_scan.end();
						++citer)
			{
				// check if we have another cluster close by
				//TableIteratorType entry_to_insert = checkInPreviousScans_(*citer,currscan_index);
				
				CoordinateType curr_mz = traits_->getPeakMz( make_pair(currscan_index,citer->first) );
				IsotopeClusterScoredCharge isoclust;
        isoclust.scans_.push_back( currscan_index );
				isoclust.first_scan_ = currscan_index;
        //TableIteratorType entry_to_insert = iso_map_.insert( TableType::value_type(curr_mz, isoclust) );
 				pair<CoordinateType,IsotopeClusterScoredCharge> entry = make_pair(curr_mz, isoclust);
				
				#ifdef DEBUG_FEATUREFINDER
			  cout << "Cluster at (" << traits_->getPeakRt( make_pair( currscan_index, citer->first ) );
				cout << " / " << traits_->getPeakMz( make_pair( currscan_index, citer->first ) ) << ")" << endl;
				cout << "Charge estimate: " << citer->second.first << " score " << citer->second.second << endl;
				cout << "Extending...." << endl;
				#endif
				
				UInt this_peak                =  citer->first;
				//CoordinateType start_mz = traits_->getPeakMz( make_pair(currscan_index,this_peak) );
				CoordinateType mz_dist  = 0;
				
				vector<UInt> points;
				points.push_back(this_peak);
													
				// walk to the left (for at most 1.0 Th)
				while (mz_dist < 1.0 && this_peak >= 1)
				{
					if ( traits_->getPeakFlag( make_pair(currscan_index,this_peak) ) == FeaFiTraits::UNUSED )
					{
						points.push_back(this_peak);
						traits_->getPeakFlag( make_pair(currscan_index,this_peak) ) = FeaFiTraits::USED;
					}
					--this_peak;
					mz_dist = ( curr_mz - traits_->getPeakMz( make_pair(currscan_index,this_peak) ) );
				}
			
				// reset
				this_peak =  (citer->first+1);
				mz_dist   = ( traits_->getPeakMz( make_pair(currscan_index,this_peak) )  - curr_mz );
        
        // and to the right (we walk for at most 4.0 m/z )
        CoordinateType dist_to_right = 4.0; // / (double) citer->second.first;
				while (mz_dist <= dist_to_right && this_peak < current_scan.size() )
				{
					if ( traits_->getPeakFlag( make_pair(currscan_index,this_peak) )  == FeaFiTraits::UNUSED )
					{
						points.push_back(this_peak);
						traits_->getPeakFlag( make_pair(currscan_index,this_peak) ) = FeaFiTraits::USED;
					}					
					++this_peak;
					mz_dist = ( traits_->getPeakMz( make_pair(currscan_index,this_peak) )  - curr_mz );
					
				}
				
				if ( points.size() > (UInt) param_.getValue("min_number_points_per_scan") )
				{
						// store charge estimate and score in this scan
					entry.second.scored_charges_.push_back( citer->second );
					// store points
					for (vector<UInt>::const_iterator citer = points.begin(); citer != points.end(); ++citer)
					{
						entry.second.peaks_.insert( make_pair(currscan_index,*citer) );
					}
				}
							
				iso_map_.push_back(entry);		
			}
					
		}	// end loop for all scans
		sort(iso_map_.begin(),iso_map_.end(),TableEntryComparator());
		#ifdef DEBUG_FEATUREFINDER 
		cout << "-----------------------------------------------------------" << endl;
		cout << "List of seeding regions: (before first merging) " << endl;
		for (TableConstIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{
			cout << "m/z " << iter->first << " int(m/z) " << int(iter->first) << " charge: " << iter->second.peaks_.charge_ ;
			cout << " first scan " << iter->second.first_scan_ << /*" last scan " << iter->second.last_scan_ <<*/ endl;
			
			for (vector<UInt>::const_iterator citer = 	iter->second.scans_.begin(); 
						citer != iter->second.scans_.end();
						++citer)
			{
				cout << "# scan : " << *citer << " (";
				IDX tmp;
				tmp.first = *citer;
				cout << traits_->getPeakRt(tmp) << ")" << endl;				
			}
		
		}
		#endif
		
		ofstream of("seed_counts");
		of << iso_map_.size() << endl;
		of.close();
		
		findNeighbours_();
	
		
		// filter hash entries (by number of scans and number of points in the cluster)
    filterHash_();    
    		
    // determine most likely charge state(s) by majority voting
    voteForCharge_();
		
		// debug output of all seeding regions with charge
		#ifdef DEBUG_FEATUREFINDER 
		cout << "-----------------------------------------------------------" << endl;
		cout << "List of seeding regions: " << endl;
		for (TableConstIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{
			cout << "m/z " << iter->first << " charge: " << iter->second.peaks_.charge_ ;
			cout << " first scan " << iter->second.first_scan_ << /*" last scan " << iter->second.last_scan_ <<*/ endl;
			
			for (vector<UInt>::const_iterator citer = 	iter->second.scans_.begin(); 
						citer != iter->second.scans_.end();
						++citer)
			{
				cout << "# scan : " << *citer << " (";
				IDX tmp;
				tmp.first = *citer;
				cout << traits_->getPeakRt(tmp) << ")" << endl;				
			}
		
		}
		#endif
//    		
}

void BaseSweepSeeder::deleteHashEntries_(std::vector<bool>& entries)
{
	cout << "Deleting entries..." << endl;
	TableType new_table;
	for (UInt i=0; i<entries.size();++i)
	{	
		if (!entries[i]) new_table.push_back(iso_map_[i]);
	}
	iso_map_ = new_table;
	cout << "Done ...." << endl;
}

void BaseSweepSeeder::filterForSize_()
{
		// filter for number of scans / significance
		vector<bool> entries_to_delete(iso_map_.size(),false);	
	
		UInt min_number_scans = param_.getValue("min_number_scans");
		UInt max_number_scans = param_.getValue("max_number_scans");
		
		UInt min_number_peaks = param_.getValue("min_number_peaks");
		cout << "Filtering for size: " << endl;
		// Filter point cluster
		
		for (UInt i=0; i< iso_map_.size(); ++i)
		{				
			std::vector<UInt>::iterator new_end = std::unique(iso_map_[i].second.scans_.begin(),iso_map_[i].second.scans_.end());
			iso_map_[i].second.scans_.erase(new_end,iso_map_[i].second.scans_.end());
		
			if (iso_map_[i].second.scans_.size() < min_number_scans || 
					iso_map_[i].second.scans_.size() > max_number_scans || 
			    iso_map_[i].second.peaks_.size() < min_number_peaks)
			{
				entries_to_delete[i] = true;
			}
		}	
		
		deleteHashEntries_(entries_to_delete);		
}

void BaseSweepSeeder::filterForSignificance_()
{
	if (iso_map_.size() == 0) return;

	// filter for number of scans / significance
	vector<bool> entries_to_delete(iso_map_.size(),false);	
			
	ProbabilityType alpha = param_.getValue("fdr_alpha");
	cout << "Filtering for significance with alpha: " << (alpha/iso_map_.size() ) << endl;	
	UInt c = 0;
	for (TableIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter,++c)
	{				
		if (iter->second.peaks_.max_charge_score_ > (alpha/iso_map_.size() ) )
		{
			entries_to_delete[c]=true;	
		}
	}
	
	deleteHashEntries_(entries_to_delete);	
}

void BaseSweepSeeder::mergeIsotopeCluster_(TableIteratorType& it1, const TableIteratorType& it2)
{
	// copy scans
	for (std::vector<UInt>::const_iterator scan_iter = it2->second.scans_.begin(); 
				scan_iter != it2->second.scans_.end();
				++scan_iter) 
		{
			it1->second.scans_.push_back(*scan_iter);
		}
					
		// copy peaks
		for (std::set<IDX>::const_iterator set_iter = it2->second.peaks_.begin();
					set_iter != it2->second.peaks_.end();
					++set_iter)
		{
			it1->second.peaks_.insert(*set_iter);	
		}
					
		// copy charge estimates for each scan
		for (std::vector< ScoredChargeType >::const_iterator sc_iter =  it2->second.scored_charges_.begin();
					sc_iter !=  it2->second.scored_charges_.end();
					++sc_iter)
		{
			it1->second.scored_charges_.push_back(*sc_iter);					
		}
}

void BaseSweepSeeder::findNeighbours_()
{
		// multimap is already sorted, traverse entries 
		vector<bool> seen(iso_map_.size(),false);		
		
		vector<TableIteratorType> entries_to_delete;
		UInt i=0;
		
		for (TableIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter, ++i)
		{
			if (seen[i])
			{ 
				//cout << "skipping..." << endl;
				continue;
			 }
			 
			TableIteratorType tmp_iter = iter;			
			++tmp_iter;
			UInt j = (i+1);
			
			#ifdef DEBUG_FEATUREFINDER
			cout << "Checking " << iter->first << " " << iter->second.first_scan_ << " vs" << endl; 
			#endif
			CoordinateType mz_dist = fabs(tmp_iter->first - iter->first);
			if (mz_dist > mz_tolerance_cluster_) continue;
			
			Int last_scan = iter->second.first_scan_;
			
			for (; (mz_dist <= mz_tolerance_cluster_) &&  (tmp_iter != iso_map_.end()); ++tmp_iter, mz_dist = fabs(tmp_iter->first - iter->first), ++j )
			{
				if (seen[j]) continue;
													
				Int next_scan = tmp_iter->second.first_scan_;
				#ifdef DEBUG_FEATUREFINDER
				cout << " " << tmp_iter->first;
				cout << " m/z dist " << mz_dist << " tol " << mz_tolerance_cluster_ << endl;
				cout << " scan " << tmp_iter->second.first_scan_ << endl;
				cout << " scan dist " << 	 abs(next_scan - last_scan) << " tol " << rt_tolerance_cluster_ << endl;				
				#endif
				if ( ((UInt) abs(next_scan - last_scan)) <=  rt_tolerance_cluster_) 
				{
					//cout << "MERGING..." << endl;
					// merge cluster 
					mergeIsotopeCluster_(iter,tmp_iter);
					last_scan = tmp_iter->second.first_scan_;
					
					// mark second cluster as deleted
					seen[j] = true;
					entries_to_delete.push_back(tmp_iter);
				}
				
			}	// end of while
			//cout << endl;
			//iter->second.first_scan_ = last_scan;
		} // end of for all (table entries)
		
		deleteHashEntries_(seen);
}

void BaseSweepSeeder::filterHash_()
{
		cout << iso_map_.size() << " isotopic clusters were found." << endl;
				
		filterForSize_();
		//filterForSignificance_();
		cout << iso_map_.size() << " clusters remained after filtering." << endl;
}

void BaseSweepSeeder::voteForCharge_()
{
	// charge states > 10 should rareley be encountered
	vector<UInt> charge_scores(10, 0 );

	for (TableIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
	{
		charge_scores.clear();
		
		for (std::vector< ScoredChargeType >::const_iterator scmz_iter = iter->second.scored_charges_.begin();
		     scmz_iter != iter->second.scored_charges_.end();
				++scmz_iter)
		{
			//cout << "Vote for charge " << scmz_iter->first << " score " << scmz_iter->second << endl;
      //cout << "Size of charge vector: " << charge_scores.size() << endl;
      
      if (scmz_iter->first == 0) continue; // zero <=> no charge estimate      
			
			// resize vector
			if ( (scmz_iter->first-1) >= charge_scores.size() || charge_scores.size() == 0)
			{
				charge_scores.resize( scmz_iter->first, 0);
			}
			
			++charge_scores.at( (scmz_iter->first -1) );
		} // end for ( std::vector< ScoredChargeType > )
		//cout << "Done..." << endl;
		
		// search for winning charge
		UInt max_vote = 0;
		UInt max_charge = 0;
		
		for (UInt i = 0; i < charge_scores.size(); ++i)
		{
			//cout << (i+1) << " " << charge_scores[i] << endl;
			if (charge_scores[i] > max_vote)
			{
				max_charge = (i + 1);
				max_vote     = charge_scores[i];
			}
		
		}
		//cout << "And the winner is " << max_charge << " with score " << max_vote << endl;		
		iter->second.peaks_.charge_           = max_charge;
		iter->second.peaks_.max_charge_score_ = max_vote;
	}
}

void BaseSweepSeeder::sumUp_(SpectrumType& scan, UInt current_scan_index)
{	
		for ( UInt i=current_scan_index + 1; i <= current_scan_index + scans_to_sumup_ && i < traits_->getData().size() ; ++i )
    {
				AlignAndSum_(scan,traits_->getData()[i]);
    }
}

void BaseSweepSeeder::addNextScan_(SpectrumType& scan, UInt current_scan_index)
{
	UInt next_scan_index = current_scan_index + scans_to_sumup_;	

	if (next_scan_index < traits_->getData().size() )
	{
		AlignAndSum_(scan,traits_->getData()[ next_scan_index ]);
	}
}

void BaseSweepSeeder::AlignAndSum_(SpectrumType& scan, const SpectrumType& neighbour)
{
    if (scan.size() == 0 || neighbour.size() == 0)
        return;

    UInt index_newscan = 0;
    for (SpectrumType::const_iterator p = neighbour.begin(); p != neighbour.end(); ++p)
    {
        while (scan[index_newscan].getMZ() < p->getMZ() && index_newscan < scan.size())
            ++index_newscan;

        // This seems to happen more frequently than expected -> quit the loop
        if (index_newscan >= scan.size() )
            break;

        if (index_newscan > 0)
        {
            CoordinateType left_diff   = fabs(scan[index_newscan-1].getMZ() - p->getMZ());
            CoordinateType right_diff = fabs(scan[index_newscan].getMZ() - p->getMZ());

            // check which neighbour is closer
            if (left_diff < right_diff && (left_diff < mass_tolerance_alignment_) )
            {
                scan[ (index_newscan-1) ].setIntensity( scan[ (index_newscan-1) ].getIntensity() + p->getIntensity() );
            }
            else if (right_diff < mass_tolerance_alignment_)
            {
                scan[ (index_newscan) ].setIntensity( scan[ (index_newscan) ].getIntensity() + p->getIntensity() );
            }
        }
        else // no left neighbour available
        {
            CoordinateType right_diff = fabs(scan[index_newscan].getMZ() - p->getMZ());
            if (right_diff < mass_tolerance_alignment_)
            {
                scan[index_newscan].setIntensity( scan[index_newscan].getIntensity() + p->getIntensity() );
            }
        }
    } // end for (all peaks in neighbouring scan)
		
}	// end of AlignAndSum_(....)

}	// end of namespace OpenMS

