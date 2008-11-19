### the directory name
set(directory include/OpenMS/VISUAL/UIC)

set(uic_sources
		#ParamEditor.ui
		)

set(sources)
foreach(i ${uic_sources})
	list(APPEND sources ${directory}/${i})
	message(${i})
endforeach(i)

QT4_WRAP_UI_OWN(uiced_sources ${sources})

set(OpenMS_sources ${OpenMS_sources} ${uiced_sources})

