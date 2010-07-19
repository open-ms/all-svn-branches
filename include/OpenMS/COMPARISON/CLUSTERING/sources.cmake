### the directory name
set(directory include/OpenMS/COMPARISON/CLUSTERING)

### list all header files of the directory here
set(sources_list_h
AverageLinkage.h
CentroidLinkage.h
ClusterAnalyzer.h
ClusterFunctor.h
ClusterHierarchical.h
ClusteringMethod.h
CompleteLinkage.h
EuclideanSimilarity.h
HashClustering.h
QTClustering.h
SingleLinkage.h
)

### add path to the filenames
set(sources_h)
foreach(i ${sources_list_h})
	list(APPEND sources_h ${directory}/${i})
endforeach(i)

### source group definition
source_group("Header Files\\OpenMS\\COMPARISON\\CLUSTERING" FILES ${sources_h})

set(OpenMS_sources_h ${OpenMS_sources_h} ${sources_h})

