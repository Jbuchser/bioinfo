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

# reference genome
REF=refs/${NAME}.fa

# reference annotation
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

# ----------------

genome:
	# Step 1: Download the reference genome
	mkdir -p $(dir ${REF})
	bio fetch ${ACC} -format fasta > ${REF}
	# stats on the genome
	seqkit stats ${REF}
	echo "# Reference genome REF=${REF}"
	bwa index ${REF}
	samtools faidx ${REF}

analyze: download
	# Calculate expected coverage
	GENOME_SIZE=$$(seqkit stats -T -a ${REF} | awk 'NR==2 {print $$3}')
	READ_LEN=150
	EXP_COV=$$(echo "scale=2; ${N} * $$READ_LEN / $$GENOME_SIZE" | bc)
	echo "Genome size: $$GENOME_SIZE bp"
	echo "Expected coverage: $$EXP_COV x"

bipproject:
	bio search PRJNA887926 -H --csv

fastq:
	# Identify SRR accession numbers
	mkdir -p $(dir ${READS_DIR})
	# Download the reads
	fastq-dump -X ${N} -F --outdir reads --split-files ${SRR}
	mv reads/SRR21835901_1.fastq ${R1}
	mv reads/SRR21835901_2.fastq ${R2}
	echo "Download complete"
	# run stats
	seqkit stats ${R1} ${R2}

fastqc: 
	# run FastQC
	mkdir -p ${READS_DIR}
	fastqc -o ${READS_DIR} ${R1} ${R2}

index:
	bwa index ${REF}

align:
	# make BAM directory
	mkdir -p $(dir ${BAM})
	# align reads to reference genome
	bwa mem -t 4 ${REF} ${R1} ${R2} | samtools sort > ${BAM}
	#index the BAM file
	samtools index ${BAM}

alignment_stats:
	# Expected coverage
	GENOME_SIZE=$$(seqkit stats -T -a ${REF} | awk 'NR==2 {print $$3}')
	READ_LEN=150
	NUM_READS=$$(samtools view -c ${BAM})
	EXP_COV=$$(echo "scale=2; $$NUM_READS * $$READ_LEN / $$GENOME_SIZE" | bc)
	echo "Genome size: $$GENOME_SIZE bp"
	echo "Number of reads: $$NUM_READS"
	echo "Expected coverage: $$EXP_COV x"

	# Percentage of reads aligned
	PCT_ALIGNED=$$(samtools flagstat ${BAM} | awk -F '[()]' '/mapped/ {print $$2; exit}')
	echo "Percent reads aligned: $$PCT_ALIGNED"

	# Observed coverage
	OBS_COV=$$(samtools depth ${BAM} | awk '{sum+=$$3} END {print sum/NR}')
	echo "Observed coverage: $$OBS_COV x"

	# Coverage variation (text summary)
	echo "Coverage distribution (coverage -> # of positions):"
	samtools depth ${BAM} | awk '{print $$3}' | sort | uniq -c | sort -nr
	# Save coverage to a file
	samtools depth ${BAM} | awk '{print $$3}' > coverage_values.txt

clean:
	# clean up generated files
	rm -rf ${REF} ${R1} ${R2} ${BAM} ${BAM}.bai

.PHONY: all refs fastq index align clean stats
	# create necessary directories