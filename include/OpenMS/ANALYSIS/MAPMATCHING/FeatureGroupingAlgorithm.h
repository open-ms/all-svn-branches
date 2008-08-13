// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHM_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHM_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

	/**
		@brief Base class for all feature grouping algorithms
		
		These algorithms group corresponding features in one map or across maps.
	*/
	class FeatureGroupingAlgorithm
	 : public FactoryProduct
	{
		public:
			/// Default constructor
			FeatureGroupingAlgorithm();

			/// Destructor
			virtual ~FeatureGroupingAlgorithm();

			///Applies the algorithm. The feature in the input @p maps are grouped and the output is written to the consensus map @p out
			virtual void group(const std::vector< FeatureMap<> >& maps, ConsensusMap& out)=0;

			/// Register all derived classes in this method
			static void registerChildren();
		
		private:
			///Copy constructor is not implemented -> private
			FeatureGroupingAlgorithm(const FeatureGroupingAlgorithm& );
			///Assignment operator is not implemented -> private
			FeatureGroupingAlgorithm& operator=(const FeatureGroupingAlgorithm& );
			
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHM_H
