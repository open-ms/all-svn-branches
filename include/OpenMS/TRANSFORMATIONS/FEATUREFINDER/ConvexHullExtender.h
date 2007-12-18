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

// #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Kernel_with_attributes.h>

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
		///
		typedef CGAL::Kernel_with_attributes <CGAL::cartesian<double> , CGAL::Kernel_with_attributes_uniform_attributes<int> > KernelWithAttributes; 
		/// A CGAL point
		typedef CGAL::Point_2<KernelWithAttributes> Point_2;

// 		class PointIndex 
// 			: public Point_2
// 		{
// 			public:
// 				/// Base class type
// 				typedef ConvexHullExtender::Point_2 Base;
// 		
// 			PointIndex(CoordinateType rt, CoordinateType mz,  IDX i)
// 				: Base(rt,mz),index_(i)
// 			{ }
// 			
// 			~PointIndex()
//       {}
// 
//       /// Copy constructor
//       PointIndex(const PointIndex& source)
//           : Base(source)
//       {
//         index_ = source.index_;
//       }
// 
//       ///  Assignment operator
//       PointIndex& operator = (const PointIndex& source)
//       {
//         if (this==&source)
//           return *this;
// 
// 				Base::operator=(source);
//         index_ = source.index_;
//        
//         return *this;
//       }
// 		
// 			private:
// 				IDX index_;
// 			
// 		};
		
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
  	
  	/// Computes the priority of a peak as function of intensity and distance from seed. 
  	ProbabilityType computePeakPriority_(const IDX& index);
  	
  	/// Checks the neighbours of the current for insertion into the boundary.
  	void checkNeighbour_(const IDX& index);
  	
  	/// Keeps track of peaks already included in the boundary (value is priority of peak) 
  	std::map<IDX, ProbabilityType> priorities_; 
  	
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
		
		/// Minium priority for points in the feature region (priority is function of intensity and distance to seed)
		ProbabilityType priority_threshold_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_CONVEXHULLEXTENDER_H
// 
