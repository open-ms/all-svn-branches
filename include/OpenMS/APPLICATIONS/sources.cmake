### the directory name
set(directory include/OpenMS/APPLICATIONS)

### list all filenames of the directory here
set(sources_list
INIFileEditorWindow.h
TOPPViewBase.h
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### Apply MOC compiler
QT4_WRAP_CPP(mocced_sources ${sources})

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${mocced_sources})

