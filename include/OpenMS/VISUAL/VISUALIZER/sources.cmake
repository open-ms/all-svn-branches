### the directory name
set(directory include/OpenMS/VISUAL/VISUALIZER)

### list all filenames of the directory here
set(sources_list
AcquisitionInfoVisualizer.h
AcquisitionVisualizer.h
BaseVisualizerGUI.h
ContactPersonVisualizer.h
DataProcessingVisualizer.h
DigestionVisualizer.h
DocumentIdentifierVisualizer.h
ExperimentalSettingsVisualizer.h
GradientVisualizer.h
HPLCVisualizer.h
InstrumentSettingsVisualizer.h
InstrumentVisualizer.h
IonDetectorVisualizer.h
IonSourceVisualizer.h
MassAnalyzerVisualizer.h
MetaInfoDescriptionVisualizer.h
MetaInfoVisualizer.h
ModificationVisualizer.h
PeptideHitVisualizer.h
PeptideIdentificationVisualizer.h
PrecursorVisualizer.h
ProteinHitVisualizer.h
ProteinIdentificationVisualizer.h
SampleVisualizer.h
ScanWindowVisualizer.h
SoftwareVisualizer.h
SourceFileVisualizer.h
SpectrumSettingsVisualizer.h
TaggingVisualizer.h
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

