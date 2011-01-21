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
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:data_is_centroided_already", "true", "MS1 data centroid data");
      this->defaults_.setValidStrings("ms1:data_is_centroided_already", StringList::create("true,false"));
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:precursor_detection_scan_levels", IntList::create(1), "Precursor detection scan levels");
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:max_inter_scan_distance", 0, "MS1 max inter scan distance"); // was 0.1, bug!
      this->defaults_.setMinInt(      "ms1:max_inter_scan_distance", 0); // FIXME will be changed to double
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:tr_resolution", 0.01, "MS1 LC retention time resolution");  // seems to have no effect
      this->defaults_.setMinFloat(    "ms1:tr_resolution", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:intensity_threshold", 1000.0, "FT peak detect MS1 intensity min threshold");
      this->defaults_.setMinFloat(    "ms1:intensity_threshold", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:max_inter_scan_rt_distance", 0.1, "MS1 max inter scan distance"); // seems to have no effect
      this->defaults_.setMinFloat(    "ms1:max_inter_scan_rt_distance", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:min_nb_cluster_members", 4, "FT peak detect MS1 min nb peak members");
      this->defaults_.setMinInt(      "ms1:min_nb_cluster_members", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:detectable_isotope_factor", 0.05, "Detectable isotope factor");
      this->defaults_.setMinFloat(    "ms1:detectable_isotope_factor", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:intensity_cv", 0.9, "IntensityCV");
      this->defaults_.setMinFloat(    "ms1:intensity_cv", 0);
      
      
      
      
      this->check_defaults_ =  false;
    }
    
    static unsigned int getNativeScanId(String native_id)
    {
      Size num_pos=0;
      while(!isdigit(native_id[num_pos]) && num_pos<native_id.length())
      {
        ++num_pos;
      }
      if(num_pos==native_id.length())
      {
        std::cout << "Native id could not be determined: " << native_id;
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Cannot convert native id to unsigned integer");
      }
      else
      {
        return native_id.substr(num_pos).toInt();
      }
    }
    
    virtual void run()
    {
      std::cout << "Superhirn integration\n";
      
      map_ = *(FeatureFinderAlgorithm<PeakType, FeatureType>::map_);
      
      MyMap dummyMap;
      Vec* datavec = new Vec(map_.size(), dummyMap);
      unsigned int scanId = 0;
      
      for (unsigned int s = 0; s < map_.size(); s++)
      {
        const SpectrumType& spectrum = map_[s];
        double rt = spectrum.getRT();
        
        bool orderByNativeIds = true; // TODO make this configurable by parameter
        if (orderByNativeIds)
          scanId = getNativeScanId(spectrum.getNativeID());
        else {
          scanId++;
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
        unsigned int scanIndex = scanId - 1;
        datavec->at(scanIndex) = m;
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
