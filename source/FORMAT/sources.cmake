### the directory name
set(directory source/FORMAT)


### list all filenames of the directory here
set(sources_list
Base64.C
ConsensusXMLFile.C
ControlledVocabulary.C
CVMappingFile.C
CVMappings.C
DTA2DFile.C
DTAFile.C
FASTAFile.C
FastaIterator.C
FastaIteratorIntern.C
FeatureXMLFile.C
FileHandler.C
IdXMLFile.C
InspectInfile.C
InspectOutfile.C
LibSVMEncoder.C
MascotInfile.C
MascotInfile2.C
MascotOutfile.C
MascotRemoteQuery.C
MascotXMLFile.C
MS2File.C
MSPFile.C
MzDataFile.C
MzMLFile.C
MzXMLFile.C
OMSSACSVFile.C
OMSSAXMLFile.C
PeakFileOptions.C
PeakTypeEstimator.C
PepNovoInfile.C
PepNovoOutfile.C
PepXMLFile.C
PersistentObject.C
PTMXMLFile.C
SequestInfile.C
SequestOutfile.C
TextFile.C
TransformationXMLFile.C
UnimodXMLFile.C
UniqueIdGenerator.C
XMLFile.C
XMLValidator.C
XTandemInfile.C
XTandemXMLFile.C
)


if (USE_ANDIMS)
	list(APPEND sources ANDIFile.C)
endif()

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
  message(${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

include(source/FORMAT/HANDLERS/sources.cmake)