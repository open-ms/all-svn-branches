
// Florian Zellers SuperHirn Inclusion Test

/*
 
 
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/IsotopicDist.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_elution_peak.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundControl.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Deisotoper.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMSCData.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/simple_math2.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Process_Data.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms2_info.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2_feature.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/featureLCprofile.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/feature.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LC_MS_XML_reader.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS1_feature_merger.h>
 
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PEAK_DETEC_mzXML_reader.h>
 #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FT_PeakDetectController.h>
*/

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