### the directory name
set(directory source/FILTERING/TRANSFORMERS)

### list all filenames of the directory here
set(sources_list
BernNorm.C
ComplementFilter.C
ComplementMarker.C
FilterFunctor.C
GoodDiffFilter.C
IntensityBalanceFilter.C
IsotopeDiffFilter.C
IsotopeMarker.C
LinearResampler.C
MarkerMower.C
NLargest.C
NeutralLossDiffFilter.C
NeutralLossMarker.C
Normalizer.C
ParentFilter.C
ParentPeakMower.C
PeakDensityFilter.C
PeakMarker.C
PreprocessingFunctor.C
Scaler.C
SqrtMower.C
TICFilter.C
ThresholdMower.C
WindowMower.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

#### source group definition
set(source_group_name source\\ANALYSIS\\DECHARGING)
source_group(${source_group_name} ${sources})

