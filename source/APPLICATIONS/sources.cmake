### the directory name
set(directory source/APPLICATIONS)

### list all filenames of the directory here
set(sources_list
INIFileEditorWindow.C
TOPPBase.C
TOPPViewBase.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

