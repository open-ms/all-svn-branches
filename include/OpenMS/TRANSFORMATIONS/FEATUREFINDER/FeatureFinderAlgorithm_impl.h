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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_IMPL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_IMPL_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSimplest.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSimple.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmWavelet.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmWatershed.h>
// #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmComposite.h>

namespace OpenMS
{

	template<class PeakType, class FeatureType> void FeatureFinderAlgorithm<PeakType,FeatureType>::registerChildren()
	{
		Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct
			(
			 FeatureFinderAlgorithmSimplest<PeakType,FeatureType>::getProductName(),
			 &FeatureFinderAlgorithmSimplest<PeakType,FeatureType>::create
			);
		Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct
			(
			 FeatureFinderAlgorithmSimple<PeakType,FeatureType>::getProductName(),
			 &FeatureFinderAlgorithmSimple<PeakType,FeatureType>::create
			);
		Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct
			(
			 FeatureFinderAlgorithmPicked<PeakType,FeatureType>::getProductName(),
			 &FeatureFinderAlgorithmPicked<PeakType,FeatureType>::create
			);
		Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct
			(
			 FeatureFinderAlgorithmWavelet<PeakType,FeatureType>::getProductName(),
			 &FeatureFinderAlgorithmWavelet<PeakType,FeatureType>::create
			);
		
    Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct
			(
			 FeatureFinderAlgorithmIsotopeWavelet<PeakType,FeatureType>::getProductName(),
			 &FeatureFinderAlgorithmIsotopeWavelet<PeakType,FeatureType>::create
			);
		
    Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct
			(
			 FeatureFinderAlgorithmWatershed<PeakType,FeatureType>::getProductName(),
			 &FeatureFinderAlgorithmWatershed<PeakType,FeatureType>::create
			);
		
    // Factory<FeatureFinderAlgorithm<PeakType,FeatureType> >::registerProduct
		// 	(
		// 	 FeatureFinderAlgorithmComposite<PeakType,FeatureType>::getProductName(),
		// 	 &FeatureFinderAlgorithmComposite<PeakType,FeatureType>::create
		// 	);
	}

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_IMPL_H
