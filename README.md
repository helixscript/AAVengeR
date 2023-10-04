# Overview  

Gene therapy introduces new genetic material to patient cells designed to augment the expression of genes or encode machinery to make specific changes to cellular genomes. Several therapies, specifically those that make use of retroviral vectors, can result in millions of genomic integrations throughout the genome which raises the concern of genotoxity.  Integrated vectors, by means of disrupting regulatory elements, promotor insertion, interrupting transcript splicing, can disrupt the normal transcription patterns of genes. Therapies that make use of non-intgegrating vectors, such as Adenosine-associated Virus (AAV), are not immune to integration events where episomal constructs are captured by nonhomologous end joining pathways during double strand break repair. The identification of integrated genomic positions and estimation of clonal population with specific integrations is critical to the field of molecular medicine.  

AAVengeR is a collection of software modules that provide a general solution for identifying and characterizing DNA integration events where the only requirement is that integrated elements include PCR priming sequences not commonly found in reference genomes. The software makes use of paired-end sequencing data of sheered genomic DNA libraries where one read sequences outward from integrated elements through genomic DNA junctures.  Integrated elements may include lentiviral proviruses, transposons, integrated AAV vectors and oligonucleotide captured via non-homologous end joining double strand break repair. The AAVengeR integration pipeline expects sequenced fragments to include sample specific linker sequences for sample demultiplexing and unique molecular identifier (UMI) sequences for clonal abundance estimations. (Figure 1). These elements can be added post facto to accommodate data generated with different protocols. Due to limitations of library preparations, sequencing depth, and large repetitive regions of reference genomes, AAVengeR typically will not return every integration site but rather return integrations from more abundant clones and a sampling from less abundant clones.  

Figure 1. AAVengeR paired-end sequencing library structure. 
<p align="center"><img src="figures/read_structure1.png" width="70%"></p>

# Approach  

AAVengeR provides six core modules to call integration sites from raw sequencing data (Figure 2). The demultiplex module quality trims reads and assigns reads to samples using both index barcode as well as unique linker sequences ligated onto genomic fragments. Next reads are prepared for alignment to a reference by removing linker sequences, duplicate reads, and reads aligning to the vector.  For integrated elements whose edges have an expected structure, such as retrovirus and transposons, Hidden Markov Models (HMM) are used to model pre-genomic juncture sequences. For elements with irregular or rearranged edges, such as AAV integrations, local sequence alignments against vector sequences are used to model pre-genomic juncture sequences. The next module removes sequences preceding junctures since they can influence alignments to reference genomes and then aligns reads to the reference using the BLAT aligner. BLAT was chosen over faster aligners such as BWA and Bowtie because it is more tolerant of mismatches near the ends of alignment and returns all alignments rather than those deemed likely by aligner algorithms. Since BLAT is not a paired-end read aligner, the following module uses the forward and reverse read alignments to create rationale genomic fragments where read mates align within 20KB of one another and their alignments face one another. The boundaries of genomic fragments are refined in the next module where small boundary variations arising from PCR, sequencing, and alignment errors are standardized. The final core module groups genomic fragments into integration events and estimates clonal abundances by tallying unique UMIs sequences and genomic fragment lengths. AAVengeR includes additional modules for characterizing integration events such as distances to nearest genes, predicting PCR artifacts, and characterizing the rearrangement of pre-genomic juncture sequences. 

Figure 2. AAVengeR core pipeline. 
<p align="center"><img src="figures/pipeline_overview1.png"></p>


# Implementation  

AAVengeR is written in the R programming language and is designed to run on a single server with modest resources while its modular design can be easily adapted to more distributive solutions such as Nextflow and cloud computing platforms. The software is driven by two configuration files, one that defines processing parameters for each module and a second that describes experimental samples. The sample configuration file contains sample specific details such as barcode and linker sequences, vector details, and reference genomes.  AAVengeR modules can be chained together to create custom pipelines and custom modules can be used by simply adding them to module chain lists and adding their parameters to the software configuration file.


# Installing additional AAVengeR genomes and genome annotations   

Due to GitHub size restrictions, only the sacCer3 genome and annotation files are provided 
with the software. Additional genomes and genome annotations are available: hg38, mm9, canFam3, macFas5, and GCA_009914755.4.

Use the commands below to install these genomes after updating the first two lines to reflect 
your installation path and genome of interest.  

```
export AAVengeR_GENOME='hg38'
export AAVengeR_HOME='/home/ubuntu/AAVengeR'
wget -qO- http://bushmanlab.org/data/AAVengeR/genomeData/$AAVengeR_GENOME.tar | tar x
rsync -a $AAVengeR_GENOME/ $AAVengeR_HOME/data/
```

