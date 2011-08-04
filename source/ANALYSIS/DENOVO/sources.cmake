### the directory name
set(directory source/ANALYSIS/DENOVO)

### list all filenames of the directory here
set(sources_list
DeNovoIonScoring.C
DeNovoAlgorithm.C
DeNovoPostScoring.C
DeNovoIdentification.C
MassDecomposition.C
MassDecompositionAlgorithm.C
CompNovoIdentificationBase.C
CompNovoIdentificationCID.C
CompNovoIonScoring.C
CompNovoIdentification.C
CompNovoIonScoringBase.C
CompNovoIonScoringCID.C
AntilopeIdentification.C
AntilopeAlgorithm.C
AntilopeSpectrumGraph.C
AntilopeIonScoringBayes.C
AntilopeSubgradientSolver.C
AntilopeBranchBound.C
AntilopePreprocessing.C
AntilopePMCorrect.C
AntilopeLagrangeProblem.C
AntilopeIonScoringRank.C
AntilopePostScoring.C
AntilopeBayesNetwork.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\DENOVO" FILES ${sources})

