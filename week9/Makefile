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
REF_SIZES = ${REF}.chrom.sizes

NAME ?= MRSA-USA300
ACC ?= NC_007793.1
PRJ ?= PRJNA887926
N ?= 100000
ALIGNER ?= bwa

REF = ${REF_DIR}/${NAME}.fa
GFF = ${REF_DIR}/${NAME}.gff

# ----------------------.PHONY----------------------------

.PHONY: all genome reads fastqc index align stats bigwig clean design

# ---------- Design ----------

design:
	mkdir -p $(dir ${DESIGN})
	# Step 1: Retrieve the full run metadata from NCBI SRA
	esearch -db sra -query "${PRJ}" | efetch -format runinfo > "${FULL_DESIGN}"
	@echo "Saved full metadata to ${FULL_DESIGN}"
	awk -F, 'NR==1 {
		for (i=1; i<=NF; i++) {
			if ($$i == "Run") run=i;
			if ($$i == "SampleName" || $$i == "Sample") sample=i;
		}
		print "sample,SRR";
		next;
	}
	{
		if (run && sample) print $$sample"_"$$run","$$run;
	}' "${FULL_DESIGN}" > "${DESIGN}"
	@echo "Created simplified design file: ${DESIGN}"
	@head -n 5 "${DESIGN}"

# ---------------- Genome ----------------
genome:
	@echo "📥 Downloading reference genome ${ACC}..."
	mkdir -p ${REF_DIR}
	bio fetch ${ACC} -format fasta > ${REF}
	bio fetch ${ACC} -format gff > ${GFF}
	seqkit stats ${REF}
	@echo "Reference genome downloaded and indexed in ${REF_DIR}"

# ---------------- Reads ----------------
reads:
	mkdir -p ${READS_DIR}
	echo "Fetching ${SRR} (${SAMPLE})..."
	fastq-dump -X ${N} -F --split-files --outdir ${READS_DIR} ${SRR}
	mv ${READS_DIR}/${SRR}_1.fastq ${READS_DIR}/${SAMPLE}_1.fastq
	mv ${READS_DIR}/${SRR}_2.fastq ${READS_DIR}/${SAMPLE}_2.fastq
	seqkit stats ${READS_DIR}/${SAMPLE}_1.fastq ${READS_DIR}/${SAMPLE}_2.fastq > ${READS_DIR}/${SAMPLE}_stats.txt

# ---------------- FastQC ----------------
fastqc:
	mkdir -p ${QC_DIR}
	fastqc -o ${QC_DIR} ${READS_DIR}/${SAMPLE}_1.fastq ${READS_DIR}/${SAMPLE}_2.fastq

# ---------------- Index ----------------
index:
	bwa index ${REF}
	samtools faidx ${REF}

# ---------------- Align ----------------
align:
	mkdir -p ${BAM_DIR}
	echo "Aligning ${SAMPLE}..."
	bwa mem -t 4 ${REF} ${READS_DIR}/${SAMPLE}_1.fastq ${READS_DIR}/${SAMPLE}_2.fastq | \
		samtools sort -o ${BAM_DIR}/${SAMPLE}.bam
	samtools index ${BAM_DIR}/${SAMPLE}.bam

# ---------------- Stats ----------------
stats:
	samtools flagstat ${BAM_DIR}/${SAMPLE}.bam > ${BAM_DIR}/${SAMPLE}_stats.txt

# ---------------- BigWig ----------------
bigwig:
	@if [ -z "${SAMPLE}" ]; then echo "Must specify SAMPLE"; exit 1; fi
	mkdir -p ${BW_DIR}
	# Create chrom.sizes file if it doesn't exist
	if [ ! -f ${REF}.chrom.sizes ]; then \
		echo "Creating ${REF}.chrom.sizes..."; \
		if [ ! -f ${REF}.fai ]; then \
			echo "Index file ${REF}.fai not found — running samtools faidx..."; \
			samtools faidx ${REF}; \
		fi; \
		cut -f1,2 ${REF}.fai > ${REF}.chrom.sizes; \
	fi

	echo "Converting ${SAMPLE} to BigWig..."
	bedtools genomecov -ibam ${BAM_DIR}/${SAMPLE}.bam -split -bg > ${BAM_DIR}/${SAMPLE}.bedgraph
	sort -k1,1 -k2,2n ${BAM_DIR}/${SAMPLE}.bedgraph -o ${BAM_DIR}/${SAMPLE}.bedgraph
	bedGraphToBigWig ${BAM_DIR}/${SAMPLE}.bedgraph ${REF}.chrom.sizes ${BW_DIR}/${SAMPLE}.bw
	rm -f ${BAM_DIR}/${SAMPLE}.bedgraph

# ---------------- Clean ----------------
clean:
	rm -rf ${READS_DIR} ${QC_DIR} ${BAM_DIR} ${BW_DIR}
	rm -f ${REF_DIR}/*.fai ${REF_DIR}/*.bwt ${REF_DIR}/*.pac ${REF_DIR}/*.sa ${REF_DIR}/*.ann ${REF_DIR}/*.amb
	rm -f ${DESIGN}

# ---------- Run in sequence ------------
all: reads fastqc align stats bigwig
	@echo "Pipeline complete"
