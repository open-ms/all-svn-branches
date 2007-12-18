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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ConvexHullExtender.h>

using namespace std;

namespace OpenMS
{

	ConvexHullExtender::ConvexHullExtender() 
	: BaseExtender(),
		last_pos_extracted_()
	{
    setName(getProductName());
        
		defaults_.setValue("intensity_factor",0.03f,"Influences for intensity (ion count) threshold in the feature extension. We include only raw data"
																								"points into this region if their intensity is larger than [intensity_factor * (intensity of the seed)].");

    defaultsToParam_();
	}

	ConvexHullExtender::~ConvexHullExtender()
	{
	}

  ConvexHullExtender::ConvexHullExtender(const ConvexHullExtender& rhs)
    : BaseExtender(rhs)
  {
    updateMembers_();
  }
  
  ConvexHullExtender& ConvexHullExtender::operator= (const ConvexHullExtender& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseExtender::operator=(rhs);
    
    updateMembers_();
    
    return *this;
  }

  void ConvexHullExtender::updateMembers_()
  {
		dist_mz_up_     = param_.getValue("dist_mz_up");
		dist_mz_down_ = param_.getValue("dist_mz_down");
		dist_rt_up_       = param_.getValue("dist_rt_up");
		dist_rt_down_   = param_.getValue("dist_rt_down");
	}


	const FeaFiModule::ChargedIndexSet& ConvexHullExtender::extend(const ChargedIndexSet& seed_region)
	{
		// find maximum of region (seed)
		CoordinateType sum_intensity = 0.0;
		boundary_.clear();
		region_.clear();
    running_avg_.clear();
	
		// pass on charge information
		region_.charge_ = seed_region.charge_;
				
		cout << "Seeding region size: " << seed_region.size() << endl;
		
		vector<Point_2> cgal_points;
		
		CoordinateType max_intensity = 0.0;
		IDX seed;
				
    for (IndexSet::const_iterator cit = seed_region.begin(); cit != seed_region.end(); ++cit)
    {	
			sum_intensity += traits_->getPeakIntensity( *cit);
			
			region_.insert(*cit);
			
			Point_2 p;
			p.attribute = *cit;
			cgal_points.push_back(p);
			
			if (traits_->getPeakIntensity(*cit) > max_intensity)
      {
        seed = *cit;
        max_intensity = traits_->getPeakIntensity(seed);						
			}
    }
		
		// re-compute intensity threshold 
		intensity_threshold_ = (double)param_.getValue("intensity_factor") * traits_->getPeakIntensity(seed);
		
		// compute convex hull
		vector<Point_2> cgal_result;
	  CGAL::convex_hull_2( cgal_points.begin(), cgal_points.end(), std::inserter(cgal_result, cgal_result.begin() ) );
	
		// points on convex hull determine feature boundary
		for (vector< Point_2>::iterator it = cgal_result.begin(); it != cgal_result.end(); ++it)
		{
			boundary_.push_back(it->attribute);
		}
		
		while (boundary_.size() > 0)
		{
			IDX idx = boundary_.back();
			
			// remember last extracted peak
			last_pos_extracted_[RawDataPoint2D::RT] = traits_->getPeakRt(idx);
			last_pos_extracted_[RawDataPoint2D::MZ] = traits_->getPeakMz(idx);

			// Now we explore the neighbourhood of the current peak. Points in this area are included
			// into the boundary if their intensity is not too low and they are not too
			// far away from the seed.			
			// Add position to the current average of positions weighted by intensity
			running_avg_.add(last_pos_extracted_,traits_->getPeakIntensity(idx));
			
			moveMzUp_(idx);
			moveMzDown_(idx);
			moveRtUp_(idx);
			moveRtDown_(idx);

			cout << "ConvexHullExtender: intensity " << traits_->getPeakIntensity(idx) << " threshold : " << (sqrt(sum_intensity) * 0.01) << endl;
			
			if ( traits_->getPeakIntensity(idx) >= sqrt(sum_intensity) * 0.01)
			{
				// set peak flags and add to boundary
				traits_->getPeakFlag(idx) = FeaFiTraits::USED;
				sum_intensity += traits_->getPeakIntensity(idx);
				region_.insert(idx);
			}
		}
		
    return region_;
	} // end of extend


bool ConvexHullExtender::isTooFarFromCentroid_(const IDX& index)
	{
	
		if ( index.first >= traits_->getData().size()) std::cout << "Scan index outside of map!" << std::endl;
		if ( index.second >= traits_->getData()[index.first].size() ) std::cout << "Peak index outside of scan!" << std::endl;
	
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<traits_->getData().size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<traits_->getData()[index.first].size() , "Peak index outside of scan!");

     const FeaFiTraits::PositionType2D& curr_mean = running_avg_.getPosition();
		 
    if ( traits_->getPeakMz(index) > curr_mean[RawDataPoint2D::MZ] + dist_mz_up_   ||
				 traits_->getPeakMz(index) < curr_mean[RawDataPoint2D::MZ] - dist_mz_down_ ||
				 traits_->getPeakRt(index) > curr_mean[RawDataPoint2D::RT] + dist_rt_up_   ||
				 traits_->getPeakRt(index) < curr_mean[RawDataPoint2D::RT] - dist_rt_down_ )
    {
    	//too far
			return true;
    }
		
		//close enough
		return false;
}

	void ConvexHullExtender::moveMzUp_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getNextMz(tmp);
				if (isTooFarFromCentroid_(tmp)) break;
				checkNeighbour_(tmp);
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void ConvexHullExtender::moveMzDown_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getPrevMz(tmp);
				if (isTooFarFromCentroid_(tmp))	break;
				checkNeighbour_(tmp);
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void ConvexHullExtender::moveRtUp_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;

			while (true)
			{
				traits_->getNextRt(tmp);
				if (isTooFarFromCentroid_(tmp)) break;
				checkNeighbour_(tmp);
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void ConvexHullExtender::moveRtDown_(const IDX& index)
	{
    try
    {
			IDX tmp = index;
			while (true)
			{
				traits_->getPrevRt(tmp);
				if (isTooFarFromCentroid_(tmp)) break;
				checkNeighbour_(tmp);
			}
    }
    catch(NoSuccessor)
    {
    }
	}
	
	void ConvexHullExtender::checkNeighbour_(const IDX& index)
	{
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<traits_->getData().size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<traits_->getData()[index.first].size(), "Peak index outside of scan!");
		
    // skip this point if its intensity is too low
    if (traits_->getPeakIntensity(index) <= intensity_threshold_)
		{
		 return;
		}
    if ( traits_->getPeakFlag(index) == FeaFiTraits::UNUSED)
    {
			traits_->getPeakFlag(index) = FeaFiTraits::USED;
			boundary_.push_back(index);			
		}
	}


} // end of class ConvexHullExtender

