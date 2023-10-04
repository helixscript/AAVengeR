# AAVengeR

Due to GitHub size restrictions, only the sacCer3 genome and annotation files are provided 
with the software. Additional genomes and genome annotations are available: hg38, mm9, canFam3, macFas5, and GCA_009914755.4.

Use the commands below to install these genomes after updating the first to lines to reflect 
your installation path and genome of interest.  

```
export AAVengeR_genome='hg38'
export AAVengeR_HOME='/home/ubuntu/AAVengeR'
wget -O update.tar http://bushmanlab.org/data/AAVengeR/genomeData/$AAVengeR_genome.tar
tar xf update.tar
rsync -a $AAVengeR_genome/ $AAVengeR_HOME/data/
```
