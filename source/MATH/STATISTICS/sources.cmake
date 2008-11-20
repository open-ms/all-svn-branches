### the directory name
set(directory source/MATH/STATISTICS)

### list all filenames of the directory here
set(sources_list
AsymmetricStatistics.C
AveragePosition.C
BasicStatistics.C
GammaDistributionFitter.C
GaussFitter.C
Histogram.C
LinearRegression.C
ROCCurve.C
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

