### the directory name
set(directory source/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN)

### list all filenames of the directory here
set(sources_list
BackgroundControl.cpp
BackgroundIntensityBin.cc
CentroidData.cpp
CentroidPeak.cpp
ClusteredMS2ConsensusSpectrum.cpp
Deisotoper.cpp
ExternalIsotopicDistribution.cpp
FT_PEAK_DETEC_mzXML_reader.cpp
FT_PeakDetectController.cc
IsotopicDist.cpp
LCMSCData.cpp
LC_MS.cc
LC_MS_XML_reader.cpp
LC_elution_peak.cpp
MS1_feature_merger.cc
MS2ConsensusSpectrum.cpp
MS2Fragment.cpp
MS2_feature.cc
PeptideIsotopeDistribution.cc
Process_Data.cpp
RawData.cpp
consensIsotopePattern.cc
feature.cc
featureLCprofile.cc
main.cpp
ms2_info.cc
ms_peak.cpp
simple_math2.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\TRANSFORMATIONS\\FEATUREFINDER\\SUPERHIRN" FILES ${sources})
