### the directory name
set(directory source/FORMAT/HANDLERS)

### list all filenames of the directory here
set(sources_list
  AcqusHandler.C
  FidHandler.C
	MascotXMLHandler.C
	MzDataHandler.C
	MzIdentMLHandler.C
	MzMLHandler.C
	MzXMLHandler.C
	ParamXMLHandler.C
	PTMXMLHandler.C
	TraMLHandler.C
	UnimodXMLHandler.C
	XMLHandler.C
	XTandemInfileXMLHandler.C
)

if (USE_ANDIMS) 	 
	list(APPEND sources_list ANDIHandler.C) 	 
endif()

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\FORMAT\\HANDLERS" FILES ${sources})

