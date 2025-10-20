SHELL=/bin/bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# ---------------------------
# Variables
FULL_DESIGN ?= full_design.csv
DESIGN ?= design.csv
READS_DIR = reads
REF_DIR = refs
QC_DIR = qc
BAM_DIR = bam
BW_DIR = bigwig

NAME ?= MRSA-USA300
ACC ?= NC_007793.1
PRJ ?= PRJNA646837
N ?= 100000
ALIGNER ?= bwa

REF = ${REF_DIR}/${NAME}.fa
GFF = ${REF_DIR}/${NAME}.gff

# --------------------------------------------------
.PHONY: all genome reads fastqc index align stats bigwig clean design

# Run everything in sequence
all:
	make design
	make genome
	make reads
	make fastqc
	make index
	make align
	make stats
	make bigwig
	@echo "Pipeline complete"

# ---------------- Design file ----------------
design:
	@echo "Creating ${DESIGN} from ${FULL_DESIGN}..."
	mkdir -p $(dir ${DESIGN})
	awk -F, 'NR==1 {print "sample,SRR"} NR>1 {print $$22"_"$$1","$$1}' ${FULL_DESIGN} > ${DESIGN}
	@head ${DESIGN}

# ---------------- Genome ----------------
genome:
	@echo "Downloading reference genome ${ACC}..."
	mkdir -p ${REF_DIR}
	bio fetch ${ACC} -format fasta > ${REF}
	bio fetch ${ACC} -format gff > ${GFF}
	seqkit stats ${REF}

# ---------------- Reads ----------------
reads:
	mkdir -p ${READS_DIR}
	tail -n +2 ${DESIGN} | parallel --colsep , -j 4 '
		echo "Fetching {2} ({1})";
		fastq-dump -X ${N} -F --split-files --outdir ${READS_DIR} {2};
		mv ${READS_DIR}/{2}_1.fastq ${READS_DIR}/{1}_1.fastq;
		mv ${READS_DIR}/{2}_2.fastq ${READS_DIR}/{1}_2.fastq;
		seqkit stats ${READS_DIR}/{1}_1.fastq ${READS_DIR}/{1}_2.fastq > ${READS_DIR}/{1}_stats.txt
	'

# ---------------- FastQC ----------------
fastqc:
	mkdir -p ${QC_DIR}
	tail -n +2 ${DESIGN} | parallel --colsep , -j 4 '
		fastqc -o ${QC_DIR} ${READS_DIR}/{1}_1.fastq ${READS_DIR}/{1}_2.fastq
	'

# ---------------- Index ----------------
index:
	bwa index ${REF}
	samtools faidx ${REF}

# ---------------- Align ----------------
align:
	mkdir -p ${BAM_DIR}
	tail -n +2 ${DESIGN} | parallel --colsep , -j 4 '
		set -euo pipefail
		echo "Aligning {1} ({2})..."
		bwa mem -t 4 ${REF} ${READS_DIR}/{1}_1.fastq ${READS_DIR}/{1}_2.fastq | \
		samtools sort -o ${BAM_DIR}/{1}.bam
		samtools index ${BAM_DIR}/{1}.bam
	'

# ---------------- Stats ----------------
stats:
	tail -n +2 ${DESIGN} | parallel --colsep , -j 4 '
		samtools flagstat ${BAM_DIR}/{1}.bam > ${BAM_DIR}/{1}_stats.txt
	'

# ---------------- BigWig ----------------
bigwig:
	mkdir -p ${BW_DIR}
	cut -f1,2 ${REF}.fai > ${REF}.chrom.sizes
	tail -n +2 ${DESIGN} | parallel --colsep , -j 4 '
		set -euo pipefail
		sample={1}
		bam="${BAM_DIR}/$$sample.bam"
		bed="${BAM_DIR}/$$sample.bedgraph"
		bw="${BW_DIR}/$$sample.bw"

		if [ -s "$$bam" ]; then
			echo "→ Converting $$sample to BigWig..."
			bedtools genomecov -ibam "$$bam" -split -bg > "$$bed"
			sort -k1,1 -k2,2n "$$bed" -o "$$bed"
			bedGraphToBigWig "$$bed" ${REF}.chrom.sizes "$$bw" || echo "⚠️ bedGraphToBigWig failed for $$sample"
			rm -f "$$bed"
		else
			echo "⚠️ Skipping $$sample: missing BAM"
		fi
	'

# ---------------- Clean ----------------
clean:
	rm -rf ${READS_DIR} ${QC_DIR} ${BAM_DIR} ${BW_DIR}
	rm -f ${REF_DIR}/*.fai ${REF_DIR}/*.bwt ${REF_DIR}/*.pac ${REF_DIR}/*.sa ${REF_DIR}/*.ann ${REF_DIR}/*.amb
	rm -f ${DESIGN}
