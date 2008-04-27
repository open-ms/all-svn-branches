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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>

#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <set>

namespace OpenMS
{
	class BaseSeeder;
	class BaseExtender;
	class BaseModelFitter;
	/**
		@brief Traits class for the feature finding algorithm.
		
		This class is rather an "umbrella" for the different modules / steps of the algorithm
		than a traits class in the traditional sense.

		@ingroup FeatureFinder 	
	*/
	class FeaFiTraits : public ProgressLogger
	{
		public:
			/// Int in a MSExperiment ( first index denotes rt, second m/z )
			typedef IsotopeCluster::IDX IDX;
			/// Int set
			typedef IsotopeCluster::IndexSet IndexSet;
			/// index set with associated charge estimate
			typedef IsotopeCluster::ChargedIndexSet ChargedIndexSet;
			/// Flag for each data point
		 enum Flag { UNUSED, USED };
					
			/// The LC/MS map
	    typedef MSExperimentExtern< Peak1D > MapType;
			/// A (single) MS spectrum
			typedef MapType::SpectrumType SpectrumType;
			/// The intensity of a point
	    typedef MapType::IntensityType IntensityType;
	    /// Coordinate of a point (either retention time or m/z)
	    typedef MapType::CoordinateType CoordinateType;
			/// Quality (goodness of fit) of a feature
			typedef Feature::QualityType QualityType;
				
	    /// 2D position type (needed for averagine model)
	    typedef DPosition<2> PositionType2D;
	
			/// No successor exception
	    typedef FeaFiModule::NoSuccessor NoSuccessor;
	
	    /// Default constructor
	    FeaFiTraits();
	
	    /// Destructor
	    virtual ~FeaFiTraits();
			
			/**
				@brief copy input data to external memory and update range information
				
				@p buffer_size is the size of the ring buffer used in the internal MSExperimentExtern
				
				@p sn_threshold is the minimum signal / noise threshold 
				
			*/
			template <class SpectrumIteratorType>
	    void setData(const SpectrumIteratorType& begin, const SpectrumIteratorType& end, UInt buffer_size, IntensityType sn_threshold = -1.0)
		  {
		  	map_.setBufferSize( buffer_size );
				map_.updateBuffer();
		
				{ // sn_estimator has its very own scope.
				
				SignalToNoiseEstimatorMeanIterative< > sn_estimator;
				
				UInt sc=1;
								
				for (SpectrumIteratorType it = begin; it != end; ++it)
				{
					std::cout << "Reading scan " << sc++ << std::endl;
					
					// remove empty scans and tandem spectra
					if (it->getMSLevel() == 1 && it->size() > 0) 
					{
					
						// filter for low intensity points
						MapType::SpectrumType new_spec;		
						new_spec.setRT( it->getRT() );
						//new_spec.resize( it->size() );
				
						if (sn_threshold > 0.0)
						{
							sn_estimator.init(it->begin(),it->end());
						}
						for (SpectrumType::const_iterator cpit = it->begin();
									cpit != it->end();
									++cpit)
						{
							if (sn_threshold > 0.0)
							{
								if (sn_estimator.getSignalToNoise(cpit) >= sn_threshold)
								{
									new_spec.push_back(*cpit);				
								}
							}
							else
							{
								new_spec.push_back(*cpit);											
							}										
						}
						
						//std::cout << new_spec.size() << " of " << it->size() << " points remained." << std::endl;
						map_.push_back(new_spec);				
					}
				}	
		    
				} // end of noise estimator scope
				
				// update range informations
		    map_.updateRanges();
						
				if (map_.getSize() == 0)
				{
					std::cout << "No data with MS level 1 provided. Aborting. " << std::endl;
					return;
				}
									
				std::cout << "This map contains " << map_.size() << " scans and " << map_.getSize() << " data points. " << std::endl;
		
		    // resize internal data structures
		    flags_.resize(map_.size());
				for (UInt i=0; i<map_.size(); ++i)
				{
					flags_[i].assign(map_[i].size(),FeaFiTraits::UNUSED);
				}
		  }
			
			/// Const access to LC-MS map
			inline const MapType& getData() const 
			{ 
				return map_; 
			}
				
	    /// non-mutable access flag with index @p index .
	    inline const Flag& getPeakFlag(const IDX& index) const
	    {
	    	return flags_[index.first][index.second];
	    }
	    /// mutable access flag with index @p index.
	    inline Flag& getPeakFlag(const IDX& index) 
	    { 
	    	return flags_[index.first][index.second];
	    }
		
	    /// access intensity of peak with index @p index.
	    inline IntensityType getPeakIntensity(const IDX& index) const
	    { 
				//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
		    OPENMS_PRECONDITION(index.second<map_[index.first].size(), "Peak index outside of scan!");
			
	    	return map_[index.first][index.second].getIntensity(); 
	    }
	    /// access m/z of peak with index @p index .
	    inline CoordinateType getPeakMz(const IDX& index) const
	    { 
	    	return map_[index.first][index.second].getMZ(); 
	    }
	    /// access retention time of peak with index @p index.
	    inline CoordinateType getPeakRt(const IDX& index) const
	    { 
	    	return map_[index.first].getRT();
	    }
	
	    /// returns the 2D coordinates of a peak (needed for models)
	    inline PositionType2D getPeakPos(const IDX& index) const
			{ 
				return PositionType2D(map_[index.first].getRT(),map_[index.first][index.second].getMZ());
			}
	
	    /// fills @p index with the index of next peak in m/z dimension
	    inline void getNextMz(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
		  	//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
		    OPENMS_PRECONDITION(index.second<map_[index.first].size(), "Peak index outside of scan!");
    
	    	//At the last peak of this spectrum
	      if (index.second==map_[index.first].size()-1)
	      {
	      	throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getNextMz", index);
	      }
	
	      ++index.second;
	    }
	
	    /// fills @p index with the index of previous peak in m/z dimension
	    inline void getPrevMz(IDX& index) const throw (NoSuccessor, Exception::Precondition)
	    {
		  	//Corrupt index
		  	OPENMS_PRECONDITION(index.first<map_.size(), "Scan index outside of map!");
		    OPENMS_PRECONDITION(index.second<map_[index.first].size(), "Peak index outside of scan!");
    
	      //begin of scan
	      if (index.second==0)
	      {
	      	throw NoSuccessor(__FILE__, __LINE__, "FeaFiTraits::getPrevMz", index);
				}
				
	      --index.second;
	    }
	
	    /// fills @p index with the index of nearest peak in m/z dimension in the next spectrum
	    void getNextRt(IDX& index) throw (NoSuccessor, Exception::Precondition);
	
	    /// fills @p index with the index of nearest peak in m/z dimension in the previous spectrum
			void getPrevRt(IDX& index) throw (NoSuccessor, Exception::Precondition);
			
			/// Calculates the convex hull of a index set and adds it to the feature
			void addConvexHull(const IndexSet& set, Feature& f) const;
	
	    /// run main loop
	    const FeatureMap<>& run(const std::vector<BaseSeeder*>& seeders,
	                              const std::vector<BaseExtender*>& extenders,
	                              const std::vector<BaseModelFitter*>& fitters);
	
		protected:
	  	/// Writes gnuplot output (only for debugging purposes)
	    void writeGnuPlotFile_(IndexSet peaks, bool last,UInt nr_feat);
	
	    /// Container for peak data
	    MapType map_;
	
	    /// Flags indicating whether a peak is unused, a seed or inside a feature region
	    std::vector< std::vector<Flag> > flags_;
	
	    /// The features found in the LC/MS map
	    FeatureMap<> features_;
	};

namespace Internal
{
	/// makes operator*() return intensity of the corresponding peak
	struct IntensityIterator 
		: FeaFiTraits::IndexSet::const_iterator
	{
    IntensityIterator ( FeaFiTraits::IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
      : FeaFiTraits::IndexSet::const_iterator(iter),
      	traits_(traits)
    {
    }
    
    FeaFiTraits::IntensityType operator * () const throw()
    {
    	return traits_->getPeakIntensity( FeaFiTraits::IndexSet::const_iterator::operator *() );
    }
    
		protected:
	    FeaFiTraits const * traits_;
	};
	
	/// Makes operator*() return mz of the corresponding peak
	struct MzIterator 
		: FeaFiTraits::IndexSet::const_iterator
	{
    MzIterator ( FeaFiTraits::IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
			: FeaFiTraits::IndexSet::const_iterator(iter), 
				traits_(traits)
    {
    }
    
    FeaFiTraits::CoordinateType operator * () const throw()
    {
    	return traits_->getPeakMz( FeaFiTraits::IndexSet::const_iterator::operator *() );
    }
    
		protected:
	    FeaFiTraits const * traits_;
	};
	
	/// Makes operator*() return retention time of the corresponding peak
	struct RtIterator 
		: FeaFiTraits::IndexSet::const_iterator
	{
    RtIterator ( FeaFiTraits::IndexSet::const_iterator const & iter, FeaFiTraits const * traits )
    	: FeaFiTraits::IndexSet::const_iterator(iter),
    		traits_(traits)
    {
    }
    
    FeaFiTraits::CoordinateType operator * () const throw()
    {
    	return traits_->getPeakRt( FeaFiTraits::IndexSet::const_iterator::operator *() );
    }
		
		protected:
	    FeaFiTraits const * traits_;
};

} // namespace Internal

}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFITRAITS_H
