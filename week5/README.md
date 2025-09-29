I have updated my data from Week4 for this README file to include an executable bash script.

## 1. Downloading the genome

```bash
#!/bin/bash

# set flags for error handling
set -UEXo pipefail

# 1. The accession number for the genome (methycillin-resistant *Staphylococcus aureus* paper) is NC_007793.

# 2. I downloaded and renamed the genomic data for NC_007793:

echo "Download genome NC_007793.1 from NCBI"
# download FASTA file for genome NC_007793.1
bio fetch NC_007793.1 -format fasta > mrsa.fna

#download GFF file for genome NC_007793.1
bio fetch NC_007793.1 -format gff > mrsa.gff
echo "Download complete"

# 2. Visualizing the genome 

# The genome is 2872769 bp long:
echo "Genome length (bp):"
awk '!/^>/ {sum += length($0)} END {print sum}' mrsa.fna

#2. Obtain sequence features:

echo "sequence features in GFF file:"
grep -v '^#' mrsa.gff | awk -F'\t' 'NF==9 {print $3}' | sort | uniq -c | sort -nr

#3. SAUSA300_RS072 is the longest gene with 31266 bps. <br>
#RNA-seq-based transcriptome analysis states SAUSA300_RS072 is also annotated as *ilvA* and # is a threonine dehydratase (Im *et al*., 2022)

#3. Get the longest gene:

echo "Longest gene:"
grep -v '^#' mrsa.gff | awk -F'\t' '$3=="gene" {len=$5-$4+1; if(len>max){max=len; id=$9}} END {print "Longest gene:", id, "Length:", max}'

# 5. Convert the GFF file to BED to merge overlapping CDS to get the coding density of the genome (86.8056%). <br> 

echo "Longest gene:"
grep -v '^#' mrsa.gff | awk -F'\t' '$3=="gene" {len=$5-$4+1; if(len>max){max=len; id=$9}} END {print "Longest gene:", id, "Length:", max}'

# 6. merge overlapping regions 
echo "Merge and sum regions"
bedtools merge -i cds.bed > cds_merged.bed
#sum CDS lengths"
awk '{sum+=$3-$2} END {print "Total coding bp:", sum}' cds_merged.bed

# 7. approximate coding density
awk 'BEGIN{genome=2880000; coding=2500000; print "Coding density:", coding/genome*100 "%"}'


# 8.  Identify the SRR accession numbers

bio search PRJNA887926 -H --csv

#The BioProject number is PRJNA887926 and the 6 associated SRR accession numbers are:
#- SRR21835896
#- SRR21835897
#- SRR21835898
#- SRR21835899
#- SRR21835900
#- SRR21835901

The output tells me there are 15.5 million reads and 3127.7 million sequenced base pairs associated with this genome. I can also see that the library contains paired end reads, and the platform used was Illumina (Illumina NovaSeq 6000). 

# 9. Download a smaller version of the run for a mini analysis and rename the split files (paired end reads) 

mkdir -p reads 
fastq-dump -X 1000 -F --outdir reads --split-files SRR21835901
mv reads/SRR21835901_1.fastq reads/myreads_1.fastq
mv reads/SRR21835901_2.fastq reads/myreads_2.fastq

# 10.  Run stats on one of the files:

 seqkit stats reads/myreads_1.fastq

 Output:

 file                   format  type  num_seqs  sum_len  min_len  avg_len  max_len
reads/myreads_1.fastq  FASTQ   DNA      1,000  101,000      101      101      101

## 11. 

To obtain 10x coverage, I need to download 300,000 reads (Coverage = total bases sequenced/size of the genome). The genome size is 2.87 Mb bps, and 2.87 Mb x 10 = 28.7 Mb. 

My sequencing reads from the smaller version are 101 bps, so 28.7 Mb / 101 bps = 284,000 reads or ~300,000 reads to round up. 

# 12. Download 300,000 reads for ~10x coverage
fastq-dump -X 300000 --split-files SRR21835901 -O reads/

# 13. Search the SRA database for another genome

# I searched for another genome of *S. aureus* NC_007793 and obtained one sequenced using the Oxford NanoPore Sequencer (BioProject PRJDB20339)

seqkit stats DRR660590_1.fastq
```

The output for the second genome using seqkit stats shows longer average length of each read (4,593.1 bps) vs the average length of the Illumina sequenced run (101 bps). This make sense because Nanopore is a long-read sequencer and Illumina is short-read. 

Another difference is the Nanopore library is genomic from Whole Genome Sequencing, while the Illumina library is transcriptomic from RNA-seq. The Nanopore reads are also single reads (paired for Illumina).

```bash
file               format  type  num_seqs    sum_len  min_len  avg_len  max_len
DRR660590_1.fastq  FASTQ   DNA      1,000  4,593,090      181  4,593.1   46,675
```
