### the directory name
set(directory source/VISUAL)

### include subdirectories
add_subdirectory(DIALOGS)
add_subdirectory(VISUALIZER)

#set(moc_sources
    #MascotRemoteQuery.C)


		#QT4_WRAP_CPP(mocced_sources ${moc_sources})


### list all filenames of the directory here
set(sources_list
		AxisTickCalculator.C
		AxisWidget.C
		ColorSelector.C
		EnhancedTabBar.C
		HistogramWidget.C
		LayerData.C
		MetaDataBrowser.C
		MultiGradient.C
		MultiGradientSelector.C
		ParamEditor.C
		PeakIcon.C
		Spectrum1DCanvas.C
		Spectrum1DWidget.C
		Spectrum2DCanvas.C
		Spectrum2DWidget.C
		Spectrum3DCanvas.C
		Spectrum3DOpenGLCanvas.C
		Spectrum3DWidget.C
		SpectrumCanvas.C
		SpectrumWidget.C)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
  message(${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources} PARENT_SCOPE)

set(uic_sources
		ParamEditor.ui)

set(sources)
foreach(i ${uic_sources})
	list(APPEND sources ${directory}/${i})
	message(${i})
endforeach(i)

set(UI_sources ${UI_sources} ${sources} PARENT_SCOPE)



