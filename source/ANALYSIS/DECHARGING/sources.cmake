### the directory name
set(directory source/ANALYSIS/DECHARGING)

### list all filenames of the directory here
set(sources_list
FeatureDecharger.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group(source\\ANALYSIS\\DECHARGING FILES ${sources})

