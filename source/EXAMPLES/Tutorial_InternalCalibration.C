#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  InternalCalibration ic;
  PeakMap exp_raw;
  MzMLFile mzml_file;
  mzml_file.load("data/Tutorial_InternalCalibration.mzML",exp_raw);

  std::vector<double> ref_masses;
  ref_masses.push_back(1296.68476942);
  ref_masses.push_back(2465.19833942);

  Param param;
  param.setValue("PeakPicker:thresholds:peak_bound",800.0);
  param.setValue("PeakPicker:thresholds:fwhm_bound",0.1);
  param.setValue("PeakPicker:wavelet_transform:scale",0.12);
  ic.setParameters(param);
  
  ic.calibrate(exp_raw,ref_masses);
  
  return 0;
} //end of main
