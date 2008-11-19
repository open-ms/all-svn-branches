### the directory name
set(directory include/OpenMS/VISUAL)

### list all filenames of the directory here
set(sources_list
AxisWidget.h
ColorSelector.h
EnhancedTabBar.h
HistogramWidget.h
MetaDataBrowser.h
MultiGradientSelector.h
ParamEditor.h
Spectrum1DCanvas.h
Spectrum1DWidget.h
Spectrum2DCanvas.h
Spectrum2DWidget.h
Spectrum3DCanvas.h
Spectrum3DOpenGLCanvas.h
Spectrum3DWidget.h
SpectrumCanvas.h
SpectrumWidget.h
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

