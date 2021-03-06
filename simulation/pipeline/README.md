###Pipeline for simulating Hi-C experiments.
HPC simulation pipeline for Hi-C experiments.

###Prerequisites and Dependencies

####Due to external binaries

- Operating System: Linux x86_64
- Perl 5
- Python 2.7

Perl5, Python 2.7 are required and we reccommend pip to makes installing Python modules very easy. At this point, scons might be simpler to install using your distribution package management rather than pip. Perl5 and Python are frequently installed with any Linux distribution flavour by default, however the pipeline requires Python 2.7 and not the more recent Python 3+.

For with apt-get based package management (Ubuntu, Debian):

```bash
sudo apt-get install perl5 python2.7 scons pip
```

For yum based package management (Redhat, Centos, Scientific):

```bash
sudo yum install perl5 python27 scons pip
```

####Python Modules

- biopython
- nestly
- networkx
- numpy
- pandas
- python-louvain
- pysam
- PyYAML
- scipy
- SCons

On all systems, we recommend using [Pip](https://pip.pypa.io/en/stable/installing/) for installation and dependency management.

For example, on a system requiring root privileges and configured for ```sudo```, try the following. 

```bash
sudo pip install -U biopython nestly networkx numpy pandas python-louvain pysam PyYAML scipy SCons
```

If installing SCons through Pip and you encounter an error, it may be necessary to install it as an Egg.

```bash
sudo pip install --egg SCons
```


It might prove easier to install each module separately if you encounter errors due to other system requirements. Consider updating pip itself if you receive a warning that it is out of date. E.g. ```pip install -U pip```

###Introduction
The basis for the pipeline is a reference sequence in raw format (ASCII nucleotides only, no header), a phylogenetic tree in newick format and an abundance profile table. This information is organised in a reference folder Eg `ref_data` and specified in the configuraiton file `config.yaml`. sgEvolver is used to generate simulated communities of a given diversity as specified by the supplied phylogenetic relationship.

Eg. A basic set up with 1 reference sequence, 2 topologies and 1 profile.
```
|ref_data
|
|---seq.raw
|
|---trees
|   |---ladder.nwk
|   |---star.nwk
|
|---tables
    |---uniform.table
```

The product of this stage are then fed to through 3 further stages in a user specified sweep of experimental conditions.

###Pipeline Stages
1. Generation of communities
2. Generation of WGS data  
3. Generation of HiC data  
4. WGS and HiC are brought together, subsequently clustered and validated.  

The pipeline is built on top of Nestly and utilises the SCons wrapper. Scons is a make-like tool permitting the definition of dependencies and therefore the management of intermediate recalculations.

Generation of each stage has been made separate to avoid excessively complicated make files in the work flow.

The stages can be run as follows, where the -j option controls concurrency.

###Running the pipeline
The following successive commands (1-4) will run the complete pipeline on a single processor.

1. **Community generation**
  - ```scons -f SConstruct_evo.py```
1. **Whole-genome (metagenome) shotgun read generation and assembly**
  - ```scons -f SConstruct_wgs.py```
1. **HiC read-pair generation**
  - ```scons -f SConstruct_hic.py```
1. **WGS and HiC read-mapping**
  - ```scons -f SConstruct_map.py```

####Concurrency

Like ```make```, ```scons``` is capable of concurrent execution by way of the ```-j (integer)``` argument. Therefore, to speed-up any stage of the workflow, specify more processors to scons. 

E.g. ```scons -j 4 -f SConstruct_evo.py``` would use four processors for the community generation stage.

####Hierarchical Dependence

As a hierarchy of dependences exists, the total pool of tasks are not independent and this limits the degree of parallelism obtainable at any point in the simulation. More simply, the four stages must be run successively.

####Execution Target
The execution target can be altered at runtime by specifying local, pbs or sge.
```
scons -j N -f SConstruct_evo.py exec_type=sge
```

###Parameter defintion
Parameters used in sweeps are modified within the configuration file `config.yaml`. Modifying nests high in the tree can result in a large cascading calculation. Any output files generated will not be deleted if their parameter sets are deleted.

###Reference genomes
The filesystem is used to organise each community, much like how the output results are produced by Nestly. The organisation is hierarchical and reflects the phylogenetic tree used in generating the community, followed by the branch length scale factor.

Each individual community is placed in a folder beneath "comm_data", where both a multi-fasta and table are required. 

Eg. For 2 communities both using a simple uniform distribution.

```
|comm_data
|
|---ladder.nwk
|   |---1.000e-1
|       |---genomes.fasta
|
|---star.nwk
    |---1.000e-1
        |---genomes.fasta
```

###Generated data
For each stage generated data and results are stored in folder hierarchies.

- `evo_data` contains communities from sgEvolver
- `wgs_data` contains WGS reads and assemblies  
- `hic_data` contains HiC reads  
- `map_data` contains all further results.  
