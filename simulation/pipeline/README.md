###Pipeline for simulating Hi-C experiments.
HPC simulation pipeline for Hi-C experiments.

####Intro
A user supplies an organised set(s) of reference genomes and community definition(s) as a table of sequence to cell mappings with relative abundance measures.

These are taken by the pipeline and a sweep of experimental conditions are explored. This happens in three stages.

####Stages
1. Generation of Communities
2. Generation of WGS data  
3. Generation of HiC data  
4. WGS and HiC are brought together, subsequently clustered and validated.  

The pipeline is built on top of Nestly and utilises the SCons wrapper. This make-like tool permits the definition of dependencies and manages the need for intermediate recalculation internally.

Generation of WGS and HiC data has been separated out to avoid its duplication within the workflow.

The stages can be run as follows, where the -j option controls concurrency.

####Running pipeline
1. scons -j N -f SConstruct_wgs.py  
2. scons -j N -f SConstruct_hic.py  
3. scons -j N -f SConstruct_map.py  

As a hierarchy of dependences exists, the total pool of tasks are not independent and limits the degree of parallelism obtainable at any point in the simulation.

####Parameter defintion
Parameter sweeps are currently modified within the SConstruct files. Modifying nests high in the tree can result in a large cascading recalculation.

####Reference genomes
The filesystem is used to organise each community, much like how the output results are produced by Nestly.

Each individual community is placed in a folder beneath "references", where both a mulit-fasta and table are required.

Eg. For 2 communities both using a simple uniform distribution.

```
|references
|
|---community1
|   |---genomes.fasta
|   |---uniform.table
|
|---community2
    |---genomes.fasta
    |---uniform.table
```

####Generated data
For each stage generated data and results are stored in folder hierarchies.

- `evo_data` contains communities from sgevolver
- `wgs_data` contains WGS reads and assemblies  
- `hic_data` contains HiC reads  
- `map_data` contains all further results.  

Matt.
