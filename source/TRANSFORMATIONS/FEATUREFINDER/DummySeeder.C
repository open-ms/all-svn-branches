// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//				   OpenMS Mass Spectrometry Framework
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/DummySeeder.h>

using namespace std;

namespace OpenMS
{
	DummySeeder::DummySeeder()
			: BaseSweepSeeder(),
				is_initialized_(false)
	{
		setName(getProductName());
		
		// minimum sn ratio for "interesting" peaks
		defaults_.setValue("min_snratio",1.1,"minimum signal-to-noise ratio for interesting peaks");
		// maximum distance to neighbouring peaks in the same scan
		defaults_.setValue("max_peak_distance",1.2,"maximum distance to neighbouring peaks in the same scan");
		
		defaultsToParam_();
	}
	
	DummySeeder::~DummySeeder()
	{
	}

  DummySeeder::DummySeeder(const DummySeeder& rhs)
    : BaseSweepSeeder(rhs),
    	is_initialized_(false)
  {
    updateMembers_();
  }
  
  DummySeeder& DummySeeder::operator= (const DummySeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSweepSeeder::operator=(rhs);
    is_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }
  
  void DummySeeder::updateMembers_()
  {
    // update member of base class first
    BaseSweepSeeder::updateMembers_();

    max_dist_mz_ = param_.getValue("max_peak_distance");
    // minimum s/n ratio for an interesting peak
    min_sn_ = param_.getValue("min_snratio");  
    
  }
	
// 	FeaFiModule::ChargedIndexSet DummySeeder::nextSeed() throw (NoSuccessor)
// 	{
// 		if (!is_initialized_)
// 		{
// 			sweep_();
// 			is_initialized_ = true;
// 		}
// 		
// 		if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 )
// 		{
// 			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, make_pair(0,0));
// 		}
// 		
// 		cout << "Retrieving next region... " << endl; 
// 			
// 		return (curr_region_++)->second.peaks_;
// 	}
	
	DummySeeder::ScoredMZVector DummySeeder::detectIsotopicPattern_(SpectrumType& current_scan )
	{
			// search for maximal positions and extract potential peaks
			vector<int> local_maxima;
			filterAndComputeLocalMax_(current_scan, local_maxima);
	
			UInt nr_maxima = local_maxima.size();
			cout << "# local maxima after filtering : " << nr_maxima << endl;
            
      ScoredMZVector scored_positions;
			
			// test for groups of local maxima resembling isotopic pattern
			for (UInt z = 0; z< nr_maxima; ++z)
			{
				// store the m/z of the current peak
				CoordinateType curr_mz = current_scan[ local_maxima[z] ].getMZ();
				
				#ifdef DEBUG_FEATUREFINDER
				cout << "Peak found ! " << endl;
		  	#endif
				
        UInt nr_peaks = 0;
        
        for (UInt c = (z+1); c < nr_maxima; ++c)
      	{
					// there are more local max, check distance
					CoordinateType dist2nextpeak = (current_scan[ local_maxima[ (c) ]  ].getMZ() - curr_mz);
				
					// collect all local maxima which are not too far away
					// and insert them into the same peak cluster (hash entry)
					if (dist2nextpeak > max_dist_mz_)	 break;
					       	
           ++nr_peaks;	
           curr_mz = current_scan[ local_maxima[c] ].getMZ();
          		
			 } // end
	
       if (nr_peaks >= min_peaks_)
       {

        #ifdef DEBUG_FEATUREFINDER
        cout << "Pattern at " << current_scan[ local_maxima[z] ].getMZ() << endl;
        cout  << "There are " << nr_peaks << " supporting local maxima. " << endl;
        #endif
        
        ScoredChargeType sc_charge;
        sc_charge.first = 0;  // no charge estimate
        sc_charge.second = 0; // no score by dummy seeder
        scored_positions.push_back( make_pair(local_maxima[z],sc_charge) );
        
        z += nr_peaks;      
       }
			
		}

  return scored_positions;
} // end of detectIsotopicPattern_()
	
	
void DummySeeder::filterAndComputeLocalMax_(const SpectrumType & vec, std::vector<int>& localmax)
{		
	SignalToNoiseEstimatorMeanIterative< SpectrumType > estimator;
	estimator.init(vec.begin(), vec.end());
		
	int i = 0;
	for(SpectrumType::const_iterator it = vec.begin(); it != vec.end(); ++it)
	{
		// Check for maximum at position i 
		if(  (i > 0) && 
				(vec[i-1].getIntensity() - vec[i].getIntensity()   < 0) &&
				(vec[i].getIntensity() - vec[i+1].getIntensity() > 0) )
		{
       if (estimator.getSignalToNoise(it) > min_sn_) localmax.push_back(i);				
		}
		++i;
	}
		
} // end filterAndComputeLocalMax_(...)

} // end of namespace OpenMS
