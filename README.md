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
  
# Usage

AAVengeR requires two configuration files. The [first configuration file](config.yml) contains the list of modules to run, module specific processing parameters, paths to resources, and the path to the second configuration file which defines sample specific parameters. The [sample configuration file](sampleData.tsv) contains sample specific information such as barcode sequences for demultiplexing, linker sequences, and reference genome against which reads should be aligned. 
  
```
aavenger.R config.yml
```

# Advanced analyses

**Vector rearrangements**

AAVengeR is typically used with data sets where sequencing data captures portions of LTR / ITR sequences followed by genomic junctures and flanking genomic DNA. AAVengeR, via its anchorReadRearrangement module, can estimate the percentage of sequencing 
reads associated with rearranged vector forms when sequencing reads face inward rather than outward. To run an inward analysis, first, update the default configuration file with the values shown below.
```
mode: manual

demultiplex_anchorReadsFile: Undetermined_S0_R2_001.fastq.gz    # Change to match your data
demultiplex_adriftReadsFile: Undetermined_S0_R1_001.fastq.gz    # Change to match your data
demultiplex_index1ReadsFile: Undetermined_S0_I1_001.fastq.gz    # Change to match your data
demultiplex_sampleDataFile:  sampleData.tsv                     # Change to match your data

modules:
  - demultiplex
  - anchorReadRearrangements

demultiplex_level: all
demultiplexs_CPUs: 15                  # Change to match your system
anchorReadRearrangements_CPUs: 15      # Change to match your system
```
In your sampleData file, it is important not to include the optional *anchorReadStartSeq* column.

Next, for each vector file defined in your sampleData file, the possible expected inward sequences need to be provided. Each vector may have more than one possible inward sequence if PCR primers can sequence inward from both ends of a vector 
or if samples were transfected with more than one vector. The length of these expected sequences should be at least as long as the sum of your sequencing cycles (R1 + R2). Below is an example of how to add your expected sequences to the
anchorReadRearrangement module.
```
anchorReadRearrangements_expectedSeqs:
  Sabatino_PMC7855056_lightChain_plasmid.fasta:
    - TTCCTTGTAGTTAATGATTAACCCGCCATGCTACTTATCTACGTAGCCATGCTCTAGGAAGATCCAGGTTAATTTTTAAAAAGCAGTCAAAAGTCCAAGTGGCCCTTGGCAGCATTTACTCTCTCTGTTTGCTCTGGTTAATAATCTCAGGAGCACAAACATTCCAGATCCAGGTTAATTTTTAAAAAGCAGTCAAAAGTCCAAGTGGCCCTTGGCAGCATTTACTCTCTCTGTTTGCTCTGGTTAATAATCTCAGGAGCACAAACATTCCAGATCCGGCGCGCCAGGGCTGGAAGCTACCTTTGACATCATTTCCTCTGCGAATGCATGTATAATTTCTACAGAACCTATTAGAAAGGATCACCCAGCCTCTGCTTTTGTACAACTTTCCCTTAAAAAACTGCCAATTCCACTGCTGTTTGGCCCAATAGTGAGAACTTTTTCCTGCTGCCTCTTGGTGCTTTTGCCTATGGCCCCTATTCTGCCTGCTGAAGACACTC
    - TTCCTTGTAGTTAATGATTAACCCGCCATGCTACTTATCTACGTAGCCATGCTCTAGGAAGATCCTTATCGATTTTACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACATCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTATCATGTCTGCTCGAAGCGGCCGCTCTAGATCAGGCGGGCTGCTGGGTGTCGCAGCCCAGGACCTCCAGCCTCAGGGCGATGTGGTGCGCCCAGCTCTGCGGGTGCAGGCGCACGTAGCGAGCCACCAGCGGGGGTTCGAGACGGTTCCGCACAGGCGTGGAGGAGTCCCGGTTTCCCTGGAAGACCTTGACTTTGCCATTCTGAAGAAACAG
  Sabatino_PMC7855056_heavyChain_plasmid.fasta:
    - TTCCTTGTAGTTAATGATTAACCCGCCATGCTACTTATCTACGTAGCCATGCTCTAGGAAGATCCAGGTTAATTTTTAAAAAGCAGTCAAAAGTCCAAGTGGCCCTTGGCAGCATTTACTCTCTCTGTTTGCTCTGGTTAATAATCTCAGGAGCACAAACATTCCAGATCCAGGTTAATTTTTAAAAAGCAGTCAAAAGTCCAAGTGGCCCTTGGCAGCATTTACTCTCTCTGTTTGCTCTGGTTAATAATCTCAGGAGCACAAACATTCCAGATCCGGCGCGCCAGGGCTGGAAGCTACCTTTGACATCATTTCCTCTGCGAATGCATGTATAATTTCTACAGAACCTATTAGAAAGGATCACCCAGCCTCTGCTTTTGTACAACTTTCCCTTAAAAAACTGCCAATTCCACTGCTGTTTGGCCCAATAGTGAGAACTTTTTCCTGCTGCCTCTTGGTGCTTTTGCCTATGGCCCCTATTCTGCCTGCTGAAGACACTC
    - TTCCTTGTAGTTAATGATTAACCCGCCATGCTACTTATCTACGTAGCCATGCTCTAGGAAGATCCTTATCGATTTTACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACATCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTATCATGTCTGCTCGAAGCGGCCGCTCAGCCCTTGTTTCTTTCTATTGCTCCACGTGAATGGTCATCGGCTCTATCTGTGGCCTCTCGGAGATCAGATAAGAACAGTCCACGTGGAGTAGGATTCTGTCCCAACAGCATCAACAAATCACTAGAGGAGACACTTTGTGCTTTAATCAGCTGTGTTCTTTCTCCAGATTGAAGGTCAATCTTCTC
```

It is important that these inward sequences match your experimental data as closely as possible. AAVengeR includes a tool for calculating the most common Kmers found in the output of its demultiplex module. The anchorReadStartSeq module can be 
run immediately after the demultiplex module and its output can be useful for checking your expected inward sequences before starting the anchorReadRearrangement module.

```
modules:
  - demultiplex
  - anchorReadStartSeqs
```

The anchorReadStartSeqs module has different settings including the which Kmers or 'windows' should be calculated.
For example, in R:

```
o <- readRDS('output/anchorReadStartSeqs/result.rds')
subset(o, sample == 'GSTP5877' & window == 100)
 trial            subject  sample      window  count    freq        seq
 Sabatino_AAV_NHP GSTP5877 GSTP5877    100     163202   0.313707022 TTCCTACACAAAAAACCAACACACAGATCTCTAGAGCTCTGATCTTTTATTGCGGCCGCTCAGTACAGATCCTGGGCCTCACAGCCCAGCACCTCCATCC
 Sabatino_AAV_NHP GSTP5877 GSTP5877    100     137422   0.264152684 TTCCTACGCGTGTCTGTCTGCACATTTCGTAGAGCGAGTGTTCCGATACTCTAATCTCCCTAGGCAAGGTTCATATTGACTTAGGTTACTTATTCTCCTT
 Sabatino_AAV_NHP GSTP5877 GSTP5877    100       1713   0.003292730 TTCCTACACAAAAACCAACACACAGATCTCTAGAGCTCTGATCTTTTATTGCGGCCGCTCAGTACAGATCCTGGGCCTCACAGCCCAGCACCTCCATCCT
 Sabatino_AAV_NHP GSTP5877 GSTP5877    100        941   0.001808791 TTCCTACACAAAAAACCAACACACAGATCTCTAGAGCTCCGATCTTTTATTGCGGCCGCTCAGTACAGATCCTGGGCCTCACAGCCCAGCACCTCCATCC
 Sabatino_AAV_NHP GSTP5877 GSTP5877    100        929   0.001785725 TTCCTACGCAAAAAACCAACACACAGATCTCTAGAGCTCTGATCTTTTATTGCGGCCGCTCAGTACAGATCCTGGGCCTCACAGCCCAGCACCTCCATCC
```

The anchorReadRearranagment module works as follows:

1. Demultiplexed R1 and R2 reads are read-in and adapter sequences are removed.

2. Overlapping fragments where R1 and R2 overlap by at least 20 NT are built. By default, only overlapping reads are kept. Overlapped reads provide longer, higher quality reads into the interior of vector bodies.

3. For each vector, the possible inward paths are defined in the configuration file. For each sample, overlapped read sequences must start with the 12 NTs from one possible inward path.

4. In order to protect against calling PCR rearrangements with gDNA as vector rearrangements, the last 12 NTs of overlapped reads must align to its corresponding vector-plasmid sequence otherwise they are removed. This typically removes ~3 - 5% of read pairs.

5. Remaining overlapped reads are clustered with CD-HIT at 90% sequencing ID. CD-HIT is a greedy clustering algorithm that starts with the longest sequences. For each vector, the expected, 500 NT inward sequences found in the configuration file are injected into the clustering and essentially remove non-rearranged forms by collating them into the first cluster.

6. The number of CD-HIT clusters, after removing the clusters associated with expected sequences, is reported as the number of rearranged forms along with the percentage of total reads associated with those forms.  
