### the directory name
set(directory source/ANALYSIS)

### include subdirectories
#add_subdirectory(DIALOGS)
#add_subdirectory(VISUALIZER)
#include(

#set(moc_sources
    #MascotRemoteQuery.C)


		#QT4_WRAP_CPP(mocced_sources ${moc_sources})


### list all filenames of the directory here
set(sources_list)


### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
  message(${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

include(source/ANALYSIS/DECHARGING/sources.cmake)
include(source/ANALYSIS/ID/sources.cmake)
include(source/ANALYSIS/QUANTITATION/sources.cmake)
include(source/ANALYSIS/SVM/sources.cmake)
include(source/ANALYSIS/MAPMATCHING/sources.cmake)

