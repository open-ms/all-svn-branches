### the directory name
set(directory source/FILTERING/BASELINE)

### list all filenames of the directory here
set(sources_list
MorphFilter.C
TopHatFilter.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

