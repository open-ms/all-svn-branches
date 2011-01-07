#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSH_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSH_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSHCtrl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

namespace OpenMS
{
  
  template<class PeakType, class FeatureType> class FeatureFinderAlgorithmSH :
    public FeatureFinderAlgorithm<PeakType, FeatureType>,
    public FeatureFinderDefs
  {
    
  public:
    typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::MapType MapType; // MSExperiment
    typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::FeatureMapType FeatureMapType;
    typedef typename MapType::SpectrumType SpectrumType;
    
    using FeatureFinderAlgorithm<PeakType, FeatureType>::features_;
    
    FeatureFinderAlgorithmSH() : FeatureFinderAlgorithm<PeakType, FeatureType>()
    {
      this->defaults_.setValue("max_inter_scan_retention_time_distance", 0.1, "MS1 max inter scan distance");
//      this->defaults_.setMin ("sweep_line:rt_votes_cutoff", 0);
      this->check_defaults_ =  false;
    }
    
    virtual void run()
    {
      std::cout << "Superhirn integration\n";
      
      map_ = *(FeatureFinderAlgorithm<PeakType, FeatureType>::map_);
      
      MyMap dummyMap;
      Vec* datavec = new Vec(map_.size(), dummyMap);
      unsigned int scanIndex = 0;
      
      for (unsigned int s = 0; s < map_.size(); s++)
      {
        const SpectrumType& spectrum = map_[s];
        double rt = spectrum.getRT();
        
        // determine native id
        //Int valid_id;
        Size num_pos=0;
        String native_id = spectrum.getNativeID();
        
        
        while(!isdigit(native_id[num_pos]) && num_pos<native_id.length())
        {
          ++num_pos;
        }
        if(num_pos==native_id.length())
        {
          std::cout << "Native id could not be determined: " << native_id;
          return;
        }
        else
        {
          scanIndex = native_id.substr(num_pos).toInt() - 1;
//          std::cout << native_id.substr(num_pos).toInt() << ", ";
        }
        
        
        vector<double>* vmzvals = new vector<double>();
        vector<double>* vintvals = new vector<double>();
        
        for (Size p=0; p<spectrum.size(); ++p)
        {
          vmzvals->push_back(spectrum[p].getMZ());
          vintvals->push_back(spectrum[p].getIntensity());
        }
        
        RawData* data = new RawData(*vmzvals, *vintvals);
        
        MyMap m;
        m[rt/60.0] = data;
        //datavec->push_back(m);
        datavec->at(scanIndex) = m;
        //scanIndex++;
      }
      
      FeatureFinderAlgorithmSHCtrl ctrl;
      ctrl.initParams(this->param_);
      std::vector<Feature> thefeatures = ctrl.extractPeaks(*datavec);
      
      for (unsigned int i=0; i<thefeatures.size(); i++)
        features_->push_back(thefeatures[i]);
    }
    
    static FeatureFinderAlgorithm<Peak1D,Feature>* create()
    {
      return new FeatureFinderAlgorithmSH();
    }
    
    static const String getProductName()
    {
      return "superhirn";
    }
    
    
  protected:
    MapType map_;
    
  };
  
}

#endif
