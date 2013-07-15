### the directory name
set(directory source/ANALYSIS/DENOVO/MSNOVOGEN)

### list all filenames of the directory here
set(sources_list
Chromosome.C
GenPool.C
Mater.C
Mutater.C
MSNovoGen.C
Scorer.C
Killer.C
DefaultScorer.C
DefaultMutater.C
DefaultMater.C
DefaultKiller.C
SwappingMutater.C
SubstitutingMutater.C
RandomMutater.C
Seeder.C
RandomSeeder.C
SequenceTagSeeder.C
InvertingMutater.C
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
  list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### source group definition
source_group("Source Files\\ANALYSIS\\DENOVO\\MSNOVOGEN" FILES ${sources})

