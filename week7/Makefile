SHELL=/bin/bash

# Execute all commands in a single shell
.ONESHELL:

# Run the shell with strict error checking.
.SHELLFLAGS = -eu -o pipefail -c

# Delete target files if the command fails.
.DELETE_ON_ERROR:

# Warn if a variable is not defined.
MAKEFLAGS += --warn-undefined-variables

# Disable built-in rules.
MAKEFLAGS += --no-builtin-rules

# ---------------------------

# Variables

# user-friendly name for genome
NAME=MRSA-USA300

#NCBI Genome accession number
ACC=NC_007793.1

# SRR accession number 
SRR=SRR21835901

# reference genome file
REF=refs/${NAME}.fa

# reference annotation file
GFF=refs/${NAME}.gff

# reads directory
READS_DIR=reads

# read 1
R1=reads/${SRR}_1.fastq

# read 2
R2=reads/${SRR}_2.fastq

# bam file
BAM=bam/${SRR}.bam

# How many reads to download
N=10000

# ------------ No changes below this line ----------------

# Step 1a: Download the reference genome
genome:
	mkdir -p $(dir ${REF})
	bio fetch ${ACC} -format fasta > ${REF}
	seqkit stats ${REF}
	echo "# Reference genome REF=${REF}"
	bwa index ${REF}
	samtools faidx ${REF}

# Step 2: Search for the BioProject to find SRR numbers
bioproject:
	bio search PRJNA887926 -H --csv

# Step 3: Download reads using SRR number and fastq-dump
fastq:
	mkdir -p $(dir ${READS_DIR})
	# Download the reads
	fastq-dump -X ${N} -F --outdir reads --split-files ${SRR}
	mv reads/SRR21835901_1.fastq ${R1}
	mv reads/SRR21835901_2.fastq ${R2}
	echo "Download complete"
	# run stats
	seqkit stats ${R1} ${R2}

# Step 4: Quality check with FastQC
fastqc: 
	mkdir -p ${READS_DIR}
	fastqc -o ${READS_DIR} ${R1} ${R2}

# Step 5: Index the reference genome with BWA
index:
	bwa index ${REF}

# Step 6: Align reads to the reference genome with BWA-MEM and convert to sorted BAM
align:
	mkdir -p $(dir ${BAM})
	bwa mem -t 4 ${REF} ${R1} ${R2} | samtools sort > ${BAM}
	samtools index ${BAM}

# Step 7: Alignment statistics
stats:
	GENOME_SIZE=$$(seqkit stats -T -a ${REF} | awk 'NR==2 {print $$3}')
	READ_LEN=150
	NUM_READS=$$(samtools view -c ${BAM})
	EXP_COV=$$(echo "scale=2; $$NUM_READS * $$READ_LEN / $$GENOME_SIZE" | bc)
	echo "Genome size: $$GENOME_SIZE bp"
	echo "Number of reads: $$NUM_READS"
	echo "Expected coverage: $$EXP_COV x"

# Step 8: Clean up intermediate files
clean:
	rm -rf ${REF} ${R1} ${R2} ${BAM} ${BAM}.bai

# Step 9: Full pipeline
all: genome bioproject fastq fastqc align alignment_stats