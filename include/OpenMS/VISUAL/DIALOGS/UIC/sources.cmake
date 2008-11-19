### the directory name
set(directory include/OpenMS/VISUAL/DIALOGS/UIC)

set(uic_sources
		DataFilterDialog.ui
		FeatureEditDialog.ui
		LayerStatisticsDialog.ui
		Spectrum1DGoToDialog.ui
		Spectrum1DPrefDialog.ui
		Spectrum2DGoToDialog.ui
		Spectrum2DPrefDialog.ui
		Spectrum3DPrefDialog.ui
		TOPPViewOpenDialog.ui
		TOPPViewPrefDialog.ui
		)

set(sources)
foreach(i ${uic_sources})
	list(APPEND sources ${directory}/${i})
	message(${i})
endforeach(i)

QT4_WRAP_UI_OWN(uiced_sources ${sources})

set(OpenMS_sources ${OpenMS_sources} ${uiced_sources})

