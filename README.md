# Installation

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
# Installation testing

```
# Build one of the provided test data set configurations.
%> ./buildRunSynSeqData.R data/configFiles/buildSynSeqData_hg38_integrase_100_sites_set1.yml

%> cat data/tests/hg38_integrase_100_sites_set1/result.tsv
      exp	     nSites	 set	percentUniqueRecovery	percentTotalRecovery	unexpectedSites	leaderSeqDist.meanleaderSeqDist.sd
   integrase	       100	  1	         90.0%	              99.0%	               0             	0	             0
```

# Implementation  

AAVengeR is written in the R programming language and is designed to run on a single server with modest resources while its modular design can be easily adapted to more distributive solutions such as Nextflow and cloud computing platforms. The software is driven by two configuration files, one that defines processing parameters for each module and a second that describes experimental samples. The sample configuration file contains sample specific details such as barcode and linker sequences, vector details, and reference genomes.  AAVengeR modules can be chained together to create custom pipelines and custom modules can be used by simply adding them to module chain lists and adding their parameters to the software configuration file.
  
# Usage

AAVengeR requires two configuration files. The [first configuration file](config.yml) contains the list of modules to run, module specific processing paramers, paths to resources, and the path to the second configuration file which defines sample specific parameters. The [sample configuration file](sampleData.tsv) contains sample specific information such as barcode sequences for demultiplexing, linker sequences, and reference genome against which reads should be aligned.  
  
```
aavenger.R config.yml
```


# Overview  

Gene therapy introduces new genetic material to patient cells designed to augment the expression of genes or encode machinery to make specific changes to cellular genomes. Several therapies, specifically those that make use of retroviral vectors, can result in millions of genomic integrations throughout the genome which raises the concern of genotoxity.  Integrated vectors, by means of disrupting regulatory elements, promotor insertion, interrupting transcript splicing, can disrupt the normal transcription patterns of genes. Therapies that make use of non-intgegrating vectors, such as Adenosine-associated Virus (AAV), are not immune to integration events where episomal constructs are captured by nonhomologous end joining pathways during double strand break repair. The identification of integrated genomic positions and estimation of clonal population with specific integrations is critical to the field of molecular medicine.  

AAVengeR is a collection of software modules that provide a general solution for identifying and characterizing DNA integration events where the only requirement is that integrated elements include PCR priming sequences not commonly found in reference genomes. The software makes use of paired-end sequencing data of sheered genomic DNA libraries where one read sequences outward from integrated elements through genomic DNA junctures.  Integrated elements may include lentiviral proviruses, transposons, integrated AAV vectors and oligonucleotide captured via non-homologous end joining double strand break repair. The AAVengeR integration pipeline expects sequenced fragments to include sample specific linker sequences for sample demultiplexing and unique molecular identifier (UMI) sequences for clonal abundance estimations. (Figure 1). These elements can be added post facto to accommodate data generated with different protocols. Due to limitations of library preparations, sequencing depth, and large repetitive regions of reference genomes, AAVengeR typically will not return every integration site but rather return integrations from more abundant clones and a sampling from less abundant clones.  

Figure 1. AAVengeR paired-end sequencing library structure. 
<p align="center"><img src="figures/read_structure1.png" width="70%"></p>

# Approach  

AAVengeR provides six core modules to call integration sites from raw sequencing data (Figure 2). The demultiplex module quality trims reads and assigns reads to samples using both index barcode as well as unique linker sequences ligated onto genomic fragments. Next reads are prepared for alignment to a reference by removing linker sequences, duplicate reads, and reads aligning to the vector.  For integrated elements whose edges have an expected structure, such as retrovirus and transposons, Hidden Markov Models (HMM) are used to model pre-genomic juncture sequences. For elements with irregular or rearranged edges, such as AAV integrations, local sequence alignments against vector sequences are used to model pre-genomic juncture sequences. The next module removes sequences preceding junctures since they can influence alignments to reference genomes and then aligns reads to the reference using the BLAT aligner. BLAT was chosen over faster aligners such as BWA and Bowtie because it is more tolerant of mismatches near the ends of alignment and returns all alignments rather than those deemed likely by aligner algorithms. Since BLAT is not a paired-end read aligner, the following module uses the forward and reverse read alignments to create rationale genomic fragments where read mates align within 20KB of one another and their alignments face one another. The boundaries of genomic fragments are refined in the next module where small boundary variations arising from PCR, sequencing, and alignment errors are standardized. The final core module groups genomic fragments into integration events and estimates clonal abundances by tallying unique UMIs sequences and genomic fragment lengths. AAVengeR includes additional modules for characterizing integration events such as distances to nearest genes, predicting PCR artifacts, and characterizing the rearrangement of pre-genomic juncture sequences. 

Figure 2. AAVengeR core pipeline. 
<p align="center"><img src="figures/pipeline_overview1.png"></p>

# Structure  

The AAvengeR data folder contains four subfolders. 

```
AAVengeR
└── data
    ├── genomeAnnotations
    │   ├── sacCer3.TUs.rds
    │   ├── sacCer3.exons.rds
    │   └── sacCer3.repeatTable.gz
    ├── hmms
    │   ├── HXB2_U5.hmm
    │   └── HXB2_U5.settings
    ├── referenceGenomes
    │   └── sacCer3.2bit
    └── vectors 
        └── HXB2.fasta
```

The hmms folder contains hmm profiles created with the [HMMER](http://hmmer.org) software package using either multiple sequence alignments or single DNA sequences as inputs and are intended to 
recognize the ends of integrated DNA elements. The file names in this folder are referenced in the [sample configuration file](sampleData.tsv). Each profile file has a coresponding settings file 
that contains the default parameters for evaulating and scoring the HMM. These settings are applied if the the *prepReads_useDefaultHMMsetting* parameter in the main configuration
file is set to *TRUE* otherwise the HMM parameters in the main configuration file are used.

Example of an HMM setting file:
```
prepReads_HMMsearchReadStartPos: 1
prepReads_HMMsearchReadEndPos:  16
prepReads_HMMmaxStartPos: 3
prepReads_HMMminFullBitScore: 5
prepReads_HMMmatchEnd: TRUE
prepReads_HMMmatchTerminalSeq: CA
```

The referenceGenomes folder contains [2bit](https://genome.ucsc.edu/goldenPath/help/twoBit.html) formated reference genomes that are referenced in the [sample configuration file](sampleData.tsv). 
These data files are created from FASTA formatted genomes using the [faToTwoBit](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64) conversion tool. 

The genomeAnnotations folder contains annotations for transription unit boundaries (*.TUs.rds) and exon boundaries (*.exons.rds). These boundaries are extracted from [UCSC genome annotations](https://hgdownload.soe.ucsc.edu)
and are stored as GenomicRange objects saved as as R rds files. 

```
GRanges object with 6125 ranges and 12 metadata columns:
         seqnames        ranges strand |       bin           name  cdsStart    cdsEnd exonCount     exonStarts       exonEnds     score       name2 cdsStartStat  cdsEndStat  exonFrames
            <Rle>     <IRanges>  <Rle> | <integer>    <character> <integer> <integer> <integer>    <character>    <character> <integer> <character>  <character> <character> <character>
     [1]     chrI     1806-2169      - |       585 NM_001180043.1      1806      2169         1          1806,          2169,         0        PAU8         cmpl        cmpl          0,
     [2]     chrI     2479-2707      + |       585 NM_001184582.1      2479      2707         1          2479,          2707,         0   YAL067W-A         cmpl        cmpl          0,
     [3]     chrI     7234-9016      - |       585 NM_001178208.1      7234      9016         1          7234,          9016,         0        SEO1         cmpl        cmpl          0,
     [4]     chrI   11564-11951      - |       585 NM_001179897.1     11564     11951         1         11564,         11951,         0     YAL065C         cmpl        cmpl          0,
     [5]     chrI   12045-12426      + |       585 NM_001180042.1     12045     12426         1         12045,         12426,         0   YAL064W-B         cmpl        cmpl          0,
     ...      ...           ...    ... .       ...            ...       ...       ...       ...            ...            ...       ...         ...          ...         ...         ...
  [6121]   chrXVI 939278-939671      - |       592 NM_001184297.1    939278    939671         1        939278,        939671,         0        ARR2         cmpl        cmpl          0,
  [6122]   chrXVI 939921-941136      + |       592 NM_001184298.1    939921    941136         1        939921,        941136,         0        ARR3         cmpl        cmpl          0,
  [6123]   chrXVI 943031-943896      + |       592 NM_001184299.1    943031    943896         2 943031,943198, 943050,943896,         0     YPR202W         cmpl        cmpl        0,1,
  [6124]   chrXVI 943879-944188      + |       592 NM_001184300.1    943879    944188         1        943879,        944188,         0     YPR203W         cmpl        cmpl          0,
  [6125]   chrXVI 944602-947701      + |       592 NM_001184301.1    944602    947701         1        944602,        947701,         0     YPR204W         cmpl        cmpl          0,```
```
<br>

Information about repeat regions is determined by the [RepeatMasker](http://www.repeatmasker.org) software package and is stored as compressed tables (*.repeatTable.gz). 

```
SW_score        percent_div     percent_del     percent_ins     query_seq       query_start     query_end       query_after     strand  repeat_name     repeat_class    repeat_start    repeat_end      repeat_after    ID      alt
34      0       0       0       chrIX   11364   11392   (428496)        +       (TA)n   Simple_repeat   1       29      (0)     1       NA
18      8.5     0       0       chrIX   22808   22832   (417056)        +       (A)n    Simple_repeat   1       25      (0)     2       NA
14      15.9    0       0       chrIX   27205   27232   (412656)        +       (TGA)n  Simple_repeat   1       28      (0)     3       NA
17      14.2    3.2     0       chrIX   28527   28557   (411331)        +       (TA)n   Simple_repeat   1       32      (0)     4       NA
19      28.2    4.9     0       chrIX   39531   39612   (400276)        +       A-rich  Low_complexity  1       86      (0)     5       NA
19      23.3    0       0       chrIX   44119   44168   (395720)        +       (ACCTCC)n       Simple_repeat   1       50      (0)     6       NA
14      9.6     0       0       chrIX   44279   44301   (395587)        +       (CCA)n  Simple_repeat   1       23      (0)     7       NA
15      11.2    3.5     0       chrIX   46885   46913   (392975)        +       (TAA)n  Simple_repeat   1       30      (0)     8       NA
```

# Database

AAVengeR is provided with a [database schema](aavenger.sql) designed to capture intermediate results and final integration sites for report generation and retrospective analyses.  The database consists of four tables:<br>
1.  demultiplex<br>
This table stores the demultiplexing data and analysis parameters provided in sample configuration files and associates the data with sequencing run identifiers.
  
2. fragments  
This table stores replicate level, non-standardized, genomic fragments where one end is anchored to integration events. Fragments are stored in a non-standardized form because standardization should be performed with all available subject data. 
  
3. multihits  
Multihit clusters that arise from reads aligning well to more than on position in a reference genome. These clustered are stored in a compressed tabular format in the database for downstream analyses, 

4. sites  
This table stores the final replicate level integration sites, estimated clonal abundances, and representative leader sequences.

Database credintals need to be stored in a local file named ~/.my.cnf and contain these fields which need to match your database configuration:
 
```
[AAVengeR]
user=admin
password=iAmAdmin
host=174.139.218.44
port=3306
database=AAVengeR
```
In order to have AAVengeR populate its database, the database group identifier in the ~/.my.cnf file needs to be included in the configuration file:

```
databaseGroup: AAVengeR
```



