This pipeline reuses the Makefile from Week 6 to generate and compare BAM and Wiggle files of Illumina and Oxford Nanopore sequencing data from *Staphylococcus aureus* strain USA300.

The entire pipeline can be executed using the customizable commands:

```bash
make all SRR=SRR21835901 ALIGNER=bwa

make all SRR=SRR33624387 ALIGNER=minimap2
```
SRR and ALIGNER are customizable to the SRR and platform specific to the sequence input. 

Indivudual targets (can be run separately) for the complete pipeline include:

```bash
make genome
```
*Downloads the reference genome*

```bash
make bioproject
```
*Uses the BioProject number to obtain SRR numbers*

```bash
make fastq
```
*Obtains SRR number*

```bash
make fastqc
```
*Runs FastQC on the downloaded FastQ reads*

```bash
align
```
*Aligns the genomes*

```bash
stats
```
*Runs statistics and calculates coverage data*
```bash
clean
```
*Cleans up intermediate files*

The SRR numbers used for Illumina (paired-end) and Oxford Nanopore (single-end) are below:

```bash
# Illumina 
make all SRR=SRR21835901

# Oxford Nanopore 
make all SRR=SRR33624387
```

Alignment for Oxford Nanopore single-end reads should use minimap2 (better for longer reads) and Illumina paired-end reads should request bwa.

#### Briefly describe the differences between the alignment in both files:

Both have good mapping coverage but different alignment patterns due to read style (long vs. short reads). 

#### Briefly compare the statistics for the two BAM files:

Illumina (SRR21835901) had 99.74% primary mapped from 20,008 total reads, while Oxford Nanopore had 96.88% from 10,791 total reads. There were 0 paired reads for the Nanopore (irrelevant since this is a single-end method). Illumina resulted in 10 singletons from unmapped mates. 

#### How many primary alignments does each of your BAM files contain?

20,000 primary alignments for Illumina and 10,000 for Nanopore (only a subset of 10,000 reads were downloaded).

#### What coordinate has the largest observed coverage (hint samtools depth)

Using:
```bash
samtools depth bam/${SRR} > depth_${platform}.txt
```

Peak coverage:
- SRR21835901 (Illumina): 1408 reads at position 2,500,332
- SRR33624387 (Nanopore): 91 reads at position 36,064

Coverage for Nanopore is lower because only 10,000 reads were downloaded. Coverage for Illumina is higher because they are paired-end and higher in some regions. 

#### Select a gene of interest. How many alignments on a forward strand cover the gene?

```bash
# Count forward strand reads
START=1437069
END=1437701
samtools view -F 16 bam/SRR21835901.bam NC_007793.1:${START}-${END} | wc -l
```
Output: 7 alignments on the forward strand (F = forward, f = reverse)