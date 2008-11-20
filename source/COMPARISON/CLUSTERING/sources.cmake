### the directory name
set(directory source/COMPARISON/CLUSTERING)

### list all filenames of the directory here
set(sources_list
AverageLinkage.C
ClusterAnalyzer.C
ClusterFunctor.C
ClusterHierarchical.C
CompleteLinkage.C
EuclideanSimilarity.C
HierarchicalClustering.C
SingleLinkage.C
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

