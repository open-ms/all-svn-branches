### the directory name
set(directory include/OpenMS/VISUAL/DIALOGS)

### list all filenames of the directory here
set(sources_list
DBOpenDialog.h
DataFilterDialog.h
FeatureEditDialog.h
HistogramDialog.h
LayerStatisticsDialog.h
SaveImageDialog.h
Spectrum1DGoToDialog.h
Spectrum1DPrefDialog.h
Spectrum2DGoToDialog.h
Spectrum2DPrefDialog.h
Spectrum3DPrefDialog.h
TOPPViewOpenDialog.h
TOPPViewPrefDialog.h
ToolsDialog.h
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

