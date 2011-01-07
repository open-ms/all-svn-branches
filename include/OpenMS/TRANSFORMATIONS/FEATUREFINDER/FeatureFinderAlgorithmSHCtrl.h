#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSHCTRL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSHCTRL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

#include <vector>
#include <map>

namespace OpenMS
{
  
  typedef std::map<double, RawData*> MyMap;
  typedef std::vector<MyMap> Vec;
  
  class FeatureFinderAlgorithmSHCtrl
  {
    
  public:
    
    FeatureFinderAlgorithmSHCtrl() { }
    
    std::vector<Feature> extractPeaks(Vec datavec);
    
    void initParams(Param param);
    
  protected:
    
  };
  
}

#endif
