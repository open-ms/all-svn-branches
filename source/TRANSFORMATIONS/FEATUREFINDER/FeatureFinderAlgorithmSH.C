// Florian Zellers SuperHirn Inclusion Test

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>

#include <vector.h>
#include <map.h>

namespace OpenMS
{ 
  
	// This probably fixes the generic template
	FeatureFinderAlgorithmSH<Peak1D,Feature> default_featurefinderalgorithmsh;
}