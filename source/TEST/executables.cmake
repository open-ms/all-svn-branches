### the directory name
set(directory source/TEST)

### list all filenames of the directory here
set(executables_list
ClassTest_test
Exception_Base_test
FactoryBase_test
Factory_test
FactoryProduct_test
SingletonRegistry_test
VersionInfo_test
FuzzyStringComparator_test
)

### add path to the filenames
set(executables)
#foreach(i ${executables_list})
#	list(APPEND executables ${directory}/${i})
#endforeach(i)

### pass source file list to the upper instance
set(TEST_executables ${TEST_executables} ${executables_list})

