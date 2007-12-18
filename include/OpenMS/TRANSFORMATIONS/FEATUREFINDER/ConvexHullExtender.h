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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_CONVEXHULLEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_CONVEXHULLEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseExtender.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Kernel_with_attributes.h>

#include <OpenMS/MATH/STATISTICS/AveragePosition.h>

#include <CGAL/Cartesian.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Point_2.h>

#include <queue>

namespace OpenMS
{

	/**
	  @brief Implements the extension phase of the FeatureFinder.
		
		@ingroup FeatureFinder
	*/
  class ConvexHullExtender 
    : public BaseExtender
  {

  public:
  
		/// Intensity of a data point
  	typedef FeaFiTraits::IntensityType IntensityType;
		/// Coordinates of a point (m/z and rt)
  	typedef FeaFiTraits::CoordinateType CoordinateType;
		/// Priority of a point (see below)
  	typedef DoubleReal ProbabilityType;
		/// Position of a point
		typedef  FeaFiTraits::PositionType2D PositionType2D;
		/// Modified CGAL traits
		typedef CGAL::Kernel_with_attributes <CGAL::Cartesian<double> , CGAL::Kernel_with_attributes_uniform_attributes<IDX> > KernelWithAttributes; 
		/// A CGAL point
		typedef CGAL::Point_2<KernelWithAttributes> Point_2;
		
  	/// Default constructor
    ConvexHullExtender();

    /// destructor
    virtual ~ConvexHullExtender();

    /// Copy constructor
    ConvexHullExtender(const ConvexHullExtender& rhs);
    
    /// Assignment operator
    ConvexHullExtender& operator= (const ConvexHullExtender& rhs);

    /// return next seed
    const ChargedIndexSet& extend(const ChargedIndexSet& seed_region);

		/// returns an instance of this class 
    static BaseExtender* create()
    {
      return new ConvexHullExtender();
    }

		/// returns the name of this module 
    static const String getProductName()
    {
      return "ConvexHullExtender";
    }
                
  protected:
  	virtual void updateMembers_();
  	
  	/// Checks if the current peak is too far from the centroid
  	bool isTooFarFromCentroid_(const IDX& current_index);
   	
   	/// Extends the seed into positive m/z direction   	
  	void moveMzUp_(const IDX& current_peak);
  	
  	/// Extends the seed into negative m/z direction 
  	void moveMzDown_(const IDX& current_peak);
  	
  	/// Extension into positive rt dimension 
  	void moveRtUp_(const IDX& current_peak);
  	
  	/// Extends the seed into negative retention time direction 
  	void moveRtDown_(const IDX& current_peak);
  	
  	/// Checks the neighbours of the current for insertion into the boundary.
  	void checkNeighbour_(const IDX& index);
  	  	
  	/// Position of last peak extracted from the boundary (used to compute the priority of neighbouring peaks)
  	PositionType2D last_pos_extracted_;
	
		/// Mininum intensity of a boundary point. Calculated from 'intensity_factor' and the seed intensity
		IntensityType intensity_threshold_;
		
		/// Maximum distance to seed in positive m/z
		CoordinateType dist_mz_up_; 
		/// Maximum distance to seed in negative m/z
		CoordinateType dist_mz_down_; 
		/// Maximum distance to seed in positive retention time
		CoordinateType dist_rt_up_; 
		/// Maximum distance to seed in negative retention time
		CoordinateType dist_rt_down_;   			
				
		std::vector<IDX> boundary_;
		
		/// keeps an running average of the peak coordinates weighted by the intensities 
  	Math::AveragePosition<2> running_avg_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_CONVEXHULLEXTENDER_H
// 
