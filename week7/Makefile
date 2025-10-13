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

NAME ?=MRSA-USA300
ACC ?=NC_007793.1
SRR ?= SRR21835901
N = 10000
ALIGNER ?= bwa
# choose "bwa" or "minimap2"

REF=refs/${NAME}.fa
GFF=refs/${NAME}.gff
READS_DIR=reads
R1=reads/${SRR}_1.fastq
R2=reads/${SRR}_2.fastq
BAM=bam/${SRR}.bam
BW=bigwig/${SRR}.bw

.PHONY: all genome fastq fastqc index align stats bigwig clean

# ------------ No changes below this line ----------------
# Targets

all: genome fastq fastqc align stats bigwig

# Step 1: Download the reference genome
genome:
	mkdir -p $(dir ${REF})
	bio fetch ${ACC} -format fasta >${REF}
	bio fetch ${ACC} -format gff >${GFF}
	seqkit stats ${REF}
	bwa index ${REF}
	samtools faidx ${REF}
	minimap2 -d ${REF}.mmi ${REF}

# Step 2: Download reads using SRR number and fastq-dump
fastq:
	mkdir -p ${READS_DIR}
	fastq-dump -X ${N} -F --outdir reads --split-files ${SRR}
	mv ${READS_DIR}/${SRR}_1.fastq ${R1} || true
	mv ${READS_DIR}/${SRR}_2.fastq ${R2} || true
	seqkit stats ${R1} ${R2} || true

# Step 3: Quality check with FastQC
fastqc: 
	mkdir -p qc
	fastqc -o qc ${R1} ${R2} || true

# Step 4: Align based on aligner type 
align:
	mkdir -p bam
	if [ "$(ALIGNER)" = "bwa" ]; then \
		bwa mem -t 4 refs/MRSA-USA300.fa reads/$*_1.fastq reads/$*_2.fastq | samtools sort > $@; \
	elif [ "$(ALIGNER)" = "minimap2" ]; then \
	if [ -f reads/${SRR}_2.fastq ]; then \
		 minimap2 --MD -ax sr refs/MRSA-USA300.fa reads/${SRR}_1.fastq reads/${SRR}_2.fastq | samtools sort > bam/${SRR}.bam; \
	else \
		minimap2 --MD -ax map-ont refs/MRSA-USA300.fa reads/${SRR}_1.fastq | samtools sort > bam/${SRR}.bam; \
	fi; \
	else \
		echo "Error: Unknown aligner '$(ALIGNER)'. Use 'bwa' or 'minimap2'."; exit 1; \
	fi
	samtools index bam/${SRR}.bam
	echo "Alignment complete: ${BAM}"

# Step 5: Alignment statistics
stats:
	samtools flagstat ${BAM} > ${BAM}_stats.txt
	cat ${BAM}_stats.txt

# Step 8: Convert BAM to BigWig for visualization
bigwig:
	mkdir -p bigwig
	bedtools genomecov -ibam ${BAM} -bg > ${BAM%.bam}.bedgraph
	sort -k1,1 -k2,2n ${BAM%.bam}.bedgraph > ${BAM%.bam}_sorted.bedgraph
	# extract chromosome sizes from the reference index
	samtools faidx ${REF}
	cut -f1,2 ${REF}.fai > ${NAME}.chrom.sizes
	bedGraphToBigWig ${BAM%.bam}_sorted.bedgraph ${NAME}.chrom.sizes ${BW}
	rm ${BAM%.bam}.bedgraph ${BAM%.bam}_sorted.bedgraph
	echo "BigWig file created: ${BW}"

# Step 8: Clean up intermediate files
clean:
	rm -rf ${REF} ${R1} ${R2} ${BAM} ${BAM}.bai
