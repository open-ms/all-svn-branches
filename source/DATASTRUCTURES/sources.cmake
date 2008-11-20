### the directory name
set(directory source/DATASTRUCTURES)

### list all filenames of the directory here
set(sources_list
BigString.C
ConstRefVector.C
ConvexHull2D.C
DBoundingBox.C
DIntervalBase.C
DPosition.C
DRange.C
DataValue.C
Date.C
DateTime.C
DefaultParamHandler.C
DistanceMatrix.C
DoubleList.C
IntList.C
Map.C
Matrix.C
Param.C
SparseVector.C
String.C
StringList.C
SuffixArray.C
SuffixArrayPeptideFinder.C
SuffixArraySeqan.C
SuffixArrayTrypticCompressed.C
SuffixArrayTrypticSeqan.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

#### source group definition
set(source_group_name source\\ANALYSIS\\DECHARGING)
source_group(${source_group_name} ${sources})

