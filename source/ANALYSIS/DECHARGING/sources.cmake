### the directory name
set(directory source/ANALYSIS/DECHARGING)

### include subdirectories
#add_subdirectory(DIALOGS)
#add_subdirectory(VISUALIZER)
#include(

#set(moc_sources
    #MascotRemoteQuery.C)


		#QT4_WRAP_CPP(mocced_sources ${moc_sources})


### list all filenames of the directory here
set(sources_list
FeatureDecharger.C)


### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
  message(${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

