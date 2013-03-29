### the directory name
set(directory source/ANALYSIS/OPENSWATH)

### list all header files of the directory here
set(sources_list
MRMDecoy.C
MRMRTNormalizer.C
TransitionTSVReader.C
OpenSwathHelper.C
ChromatogramExtractor.C
SpectrumAddition.C
MRMTransitionGroupPicker.C
DIAHelper.C
DIAScoring.C
DIAPrescoring.C
MRMFeatureFinderScoring.C
ConfidenceScoring.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\OPENSWATH" FILES ${sources})
set_source_files_properties(${directory}/sources.cmake PROPERTIES HEADER_FILE_ONLY TRUE)
