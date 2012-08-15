### the directory name
set(directory source/APPLICATIONS/UTILS)

### list all filenames of the directory here
set(UTILS_executables
CVInspector
DeMeanderize
DecoyDatabase
Digestor
DigestorMotif
ERPairFinder
FFEval
FuzzyDiff
IDExtractor
IDMassAccuracy
IDSplitter
IDDecoyProbability
IdXMLEvaluation
ImageCreator
INIUpdater
LabeledEval
MassCalculator
MRMPairFinder
MSSimulator
MapAlignmentEvaluation
OpenMSInfo
SemanticValidator
SequenceCoverageCalculator
SpecLibCreator
SvmTheoreticalSpectrumGeneratorTrainer
TransformationEvaluation
XMLValidator
QCCalculator
)

## all targets with need linkage against OpenMS_GUI.lib - they also need to appear in the list above)
set(UTILS_executables_with_GUIlib
ImageCreator
INIUpdater
)

### add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${UTILS_executables})
	list(APPEND sources_VS "${i}.C")
endforeach(i)
source_group("" FILES ${sources_VS})
