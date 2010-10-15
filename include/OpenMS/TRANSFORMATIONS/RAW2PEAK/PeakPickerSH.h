#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERSH_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERSH_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

// TODO: Check if I need this
#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING

namespace OpenMS
{

  class OPENMS_DLLAPI PeakPickerSH
		: public DefaultParamHandler,
			public ProgressLogger
  {
	public:
    PeakPickerSH();
		
    virtual ~PeakPickerSH();
		
    /** 
		 @brief Picks one peak.
		 */
    template <typename PeakType>
    void pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output, float fWindowWidth);
		
    /** 
				@brief Applies the peak-picking algorithm to a map (MSExperiment). This method picks peaks for each scan in the map consecutively. The resulting picked peaks are written to the output map.
    */
//    template <typename PeakType>
//    void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output);
    void pickExperiment(const MSExperiment<>& input, MSExperiment<>& output);
    
	protected:
		QString xmlPath;
	
	}; // end PeakPickerSH
	
}// namespace OpenMS

#endif
