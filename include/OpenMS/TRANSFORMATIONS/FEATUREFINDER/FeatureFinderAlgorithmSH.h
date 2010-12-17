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
    
    //using FeatureFinderAlgorithm<PeakType, FeatureType>::map_;
    
    typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::MapType MapType; // MSExperiment
    typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::FeatureMapType FeatureMapType;
    typedef typename MapType::SpectrumType SpectrumType;
    
    //using FeatureFinderAlgorithm<PeakType, FeatureType>::param_;
    using FeatureFinderAlgorithm<PeakType, FeatureType>::features_;
    //using FeatureFinderAlgorithm<PeakType, FeatureType>::ff_;
    //using FeatureFinderAlgorithm<PeakType, FeatureType>::defaults_;
    //using FeatureFinderAlgorithm<PeakType, FeatureType>::map_;
    
    FeatureFinderAlgorithmSH() : FeatureFinderAlgorithm<PeakType, FeatureType>()
    {
      this->defaults_ = getDefaultParameters();
      this->check_defaults_ =  false;
    }
    
    virtual Param getDefaultParameters() const
    {
      Param tmp;
      return tmp;
    }
    
    virtual void run()
    {
      std::cout << "the old new run\n";
      
      MyMap* datamap = new MyMap();  
      Vec* datavec = new Vec();
      
      map_ = *(FeatureFinderAlgorithm<PeakType, FeatureType>::map_);
      
      for (unsigned int s = 0; s < map_.size(); s++)
      {
        const SpectrumType& spectrum = map_[s];
        double rt = spectrum.getRT();
        
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
        datavec->push_back(m);
      }
      
      FeatureFinderAlgorithmSHCtrl ctrl;
      std::vector<Feature> thefeatures = ctrl.extractPeaks(*datavec);
      
      for (unsigned int i=0; i<thefeatures.size(); i++)
        features_->push_back(thefeatures[i]);
      
      //Feature f;
      //set label
      //f.setMetaValue(3,plot_nr);
      //f.setCharge(c);
      //f.setOverallQuality(final_score);
      //f.setRT(fitter->getCenter());
      //f.setMZ(mono_mz);
      //f.setIntensity(0);
      //add convex hulls of mass traces
      //for (Size j=0; j<traces.size(); ++j)
      //{
      //  f.getConvexHulls().push_back(traces[j].getConvexhull());
      //}
      //features_->push_back(f);
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
