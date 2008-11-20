### the directory name
set(directory source/COMPARISON/SPECTRA)

### list all filenames of the directory here
set(sources_list
BinnedSharedPeakCount.C
BinnedSpectralContrastAngle.C
BinnedSpectrum.C
BinnedSpectrumCompareFunctor.C
BinnedSumAgreeingIntensities.C
CompareFouriertransform.C
PeakAlignment.C
PeakSpectrumCompareFunctor.C
SpectrumAlignment.C
SpectrumAlignmentScore.C
SpectrumCheapDPCorr.C
SpectrumPrecursorComparator.C
SteinScottImproveScore.C
ZhangSimilarityScore.C
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

