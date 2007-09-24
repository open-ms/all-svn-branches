// -*- C++: make; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>

namespace OpenMS
{
	/** 
		@brief Seeding module which selects single peaks based on their s/n ratio.
		
		Groups of peaks a clustered within a certain distance and traced over consecutive scans.

    @ref DummySeeder_Parameters are explained on a separate page.

		@ingroup FeatureFinder
	*/ 
  class DummySeeder 
    : public BaseSweepSeeder
  {
  	public:	
    
      /// intensity of a peak
      typedef FeaFiTraits::IntensityType IntensityType;
      /// coordinate ( in rt or m/z )
      typedef FeaFiTraits::CoordinateType CoordinateType;
      /// score
      typedef DoubleReal ProbabilityType;

      /// a single MS spectrum
      typedef BaseSweepSeeder::SpectrumType SpectrumType;

      /// charge state estimate with associated score
      typedef BaseSweepSeeder::ScoredChargeType ScoredChargeType;
      /// m/z position in spectrum with charge estimate and score
      typedef BaseSweepSeeder::ScoredMZType ScoredMZType;
      /// container of scored m/z positions
      typedef BaseSweepSeeder::ScoredMZVector ScoredMZVector;
	
		  /// Default constructor
	    DummySeeder();
	
	    /// destructor 
	    virtual ~DummySeeder();

	    /// Copy constructor
	    DummySeeder(const DummySeeder& rhs);
	    
	    /// Assignment operator
	    DummySeeder& operator= (const DummySeeder& rhs);
		
	    static BaseSeeder* create()
	    {
	      return new DummySeeder();
	    }
	
	    static const String getProductName()
	    {
	      return "DummySeeder";
	    }
	
	  protected:
      
      /// keeps member and param entries in synchrony
      virtual void updateMembers_();
	 			
       /// Finds local maxima in the cwt
			void filterAndComputeLocalMax_(const SpectrumType & vec, std::vector<int>& localmax);
	
		  /// Detects isotopic patterns
		  ScoredMZVector detectIsotopicPattern_(SpectrumType& scan);
      
      /// Max distance in mz between consecutive peaks
      CoordinateType max_dist_mz_;
      
      /// Minimum signal/noise ratio for a peak
      IntensityType min_sn_;
      
      /// Minimum number of local maxima 
      UInt min_peaks_;
      
       /// indicates whether this module has been initialized
      bool is_initialized_;
		 
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_DUMMYSEEDER_H
