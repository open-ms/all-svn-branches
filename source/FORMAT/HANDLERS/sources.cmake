### the directory name
set(directory source/FORMAT/HANDLERS)

### list all filenames of the directory here
set(sources_list
ANDIHandler.C
ConsensusXMLHandler.C
FeatureXMLHandler.C
MascotXMLHandler.C
MzDataHandler.C
MzMLHandler.C
MzXMLHandler.C
PTMXMLHandler.C
ParamXMLHandler.C
UnimodXMLHandler.C
XMLHandler.C
XTandemInfileXMLHandler.C
)

### include ANDIHandler only if is set as option
if (USE_ANDIMS)
  list(APPEND sources ANDIHandler.C)
endif()
	

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

