### the directory name
set(directory source/KERNEL)

### list all filenames of the directory here
set(sources_list
AreaIterator.C
ConsensusFeature.C
ConsensusMap.C
DPeak.C
DRichPeak.C
DSpectrum.C
Feature.C
FeatureHandle.C
FeatureMap.C
MSExperiment.C
MSSpectrum.C
Peak1D.C
Peak2D.C
PeakIndex.C
RangeManager.C
RichPeak1D.C
RichPeak2D.C
StandardTypes.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

