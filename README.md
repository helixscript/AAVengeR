# Quick Start

Clone this repository and install one or more precompiled genomes. 

```
%> git clone https://github.com/helixscript/AAVengeR
%> cd AAVengeR

%> ./aavenger.R --list.available.genomes
   Available genomes:

   refGenome   dataSize 
   ----------- ----------
     canFam3     13 GB   
      hg38       17 GB   
       hs1       17 GB   
      mm10       15 GB   
       mm9       15 GB   
    rheMac10     16 GB   
     sacCer3     68 MB   

%> ./aavenger.R --install.genome sacCer3
%> ./aavenger.R --install.genome hg38
```

AAVengeR requires two configuration files. The [first configuration file](config.yml) contains the list of modules to run, module specific processing parameters, paths to resources, and the path to the second configuration file which defines sample specific parameters. The [sample configuration file](sampleData.tsv) contains sample specific information such as barcode sequences for demultiplexing, linker sequences, and reference genome against which reads should be aligned. Typically, only parameters near the top of the configuration file need to be changed such as those pointing to sequencing data, the sample configuration file, module chain, and the name of the output directory.  
  
AAVengeR pushes data through a series of modules defined in the 'modules' section of the configuration. There core modules are: demultiplex, prepReads, alignReads, buildFragments, buildStdFragments and buildSites. Rather than chaining these modules together in the 'modules' section, it is often more advantageous to simply call the 'core' module. This module calls the six core modules and automatically allocates CPUs to each sample depending on the number of demultiplexed reads. When using the core module, progress can be tracked through three key log files. The first log file, assuming your output directory is named 'output', would be 'output/core/demultiplex/log'. This file is updated in real time detailing the progress of the demultiplex module. Once the demultiplexing module finishes, 'output/core/replicateJobTable' will be created and show each sample replicate, the number of CPUs assigned to it, and its position and status in the queue. Finally, once replicate level analyses are complete, the 'output/cores/subjectJobTable' log will be created and report the progress of subject level analyses after which integration sites will be written to 'output/core'.
  
When both configuration files are ready, the pipeline is launched with this command:
  
```
%> ./aavenger.R config.yml
```

The AAVengeR pipeline is written in both R and Python and requires several software libraries and third party tools to run. In order to simplify its installation and standardize its behavior, a [precompiled Docker image](http://bushmanlab.org/data/AAVengeR/docker/aavenger_docker_v3.tar) is available and recommended. 

```
%> docker load < aavenger_docker_v3.tar
```
  
The Docker container expects you to 'bind' a file directory containing all the files needed for the run (AAVengeR, FASTQs, and configuration files) to the container at run time. Within the container, the directory will be bound to */data*. 
For example, if AAVengeR and your data files are all located in your home directory */home/aavenger_user*, these parameters would bind your home directory to the container when it starts:
  
```
--mount type=bind,source=/home/aavenger_user,target=/data
```
  
The container needs to know the location of AAVengeR and its configuration file. These paths are provided in the command to start the container. **Importantly, these paths, and the paths included in your AAVengeR configuration file, need to be written from the container's perspective**. For example, if you installed AAVengeR at */home/aavenger_user/AAVengeR* and bound */home/aavenger_user* to the container's */data* directory, the path to AAVengeR would be */data/AAVengeR*. Likewise, in your configuration files, if the path to your R1 FASTQ file is normally */home/aavenger_user/myRun/R1.fastq.gz*, this should be written from the container's perspective, */data/myRun/R1.fastq.gz*. Putting it all together: 

```
%> docker run --rm --mount type=bind,source=/home/aavenger_user,target=/data -e AAVENGER_DIR=/data/AAVengeR -e AAVENGER_CONFIG_PATH=/data/myRun/config.yml aavenger_docker_v3
```

The pipeline can take a couple of hours to run depending on how many CPUs are allocated to each module. It is reccomended that Docker is called from within a screen session.

