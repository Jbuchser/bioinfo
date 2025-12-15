# =========================
# Makefile — Assignment 134 RNA-Seq Differential Expression Pipeline
# =========================

SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# -------------------------
# CONFIG
# -------------------------
PRJ        ?= PRJNA887926
NAME       ?= MRSA-USA300
THREADS    ?= 4
NREADS     ?= 100000

FULL_DESIGN := full_design.csv
DESIGN      := design.csv

REF_DIR   := refs
READ_DIR  := reads
BAM_DIR   := bam
QC_DIR    := qc
RES_DIR   := results
BW_DIR    := bigwig
COV_DIR   := coverage

# Reference paths
REF        := $(REF_DIR)/$(NAME).fa
GFF        := $(REF_DIR)/$(NAME).gff
FAI        := $(REF).fai
CHROMSIZES := $(REF_DIR)/$(NAME).chrom.sizes

# HISAT2 index paths
HISAT_DIR := $(REF_DIR)/hisat2
HISAT_IDX := $(HISAT_DIR)/$(NAME)

# Sample list from design.csv
SRRS := $(shell tail -n +2 $(DESIGN) | cut -d',' -f2)

# Sample-specific vars (set by `make all SRR=...`)
SRR ?=
R1 := $(READ_DIR)/$(SRR)_1.fastq
R2 := $(READ_DIR)/$(SRR)_2.fastq
BAM := $(BAM_DIR)/$(SRR).bam

# -------------------------
# PHONY TARGETS
# -------------------------
.PHONY: all all_samples design refs index fastq fastqc align stats bigwig counts de clean

# -------------------------
# MASTER (per-sample)
# -------------------------
all: fastq fastqc align stats

# -------------------------
# ALL SAMPLES + BIGWIG + COUNTS + DE
# -------------------------
all_samples: design refs index
	@echo "===== Running per-sample pipeline ====="
	@for srr in $(SRRS); do \
		echo "========== $$srr =========="; \
		$(MAKE) all SRR=$$srr; \
	done
	@echo "===== Generating BigWig files ====="
	$(MAKE) bigwig
	@echo "===== Running featureCounts ====="
	$(MAKE) counts
	@echo "===== Running DESeq2 ====="
	$(MAKE) de

# -------------------------
# DESIGN
# -------------------------
design:
	esearch -db sra -query "$(PRJ)" | efetch -format runinfo > $(FULL_DESIGN)
	awk -F, 'NR==1{ \
		for(i=1;i<=NF;i++){ \
			if($$i=="Run") r=i; \
			if($$i=="SampleName"||$$i=="Sample") s=i \
		} \
		print "sample,condition"; next \
	} \
	{ \
		srr=$$r; \
		sample=$$s; \
		if(sample ~ /^NT/) cond="control"; \
		else if(sample ~ /^NaP/) cond="treatment"; \
		else cond="unknown"; \
		print srr","cond \
	}' $(FULL_DESIGN) > $(DESIGN)

# -------------------------
# REFERENCE
# -------------------------
refs: $(REF) $(GFF) $(FAI) $(CHROMSIZES)

$(REF):
	mkdir -p $(REF_DIR)
	curl -L \
	 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz \
	 -o $(REF).gz
	gunzip -f $(REF).gz

$(GFF):
	mkdir -p $(REF_DIR)
	curl -L \
	 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz \
	 -o $(GFF).gz
	gunzip -f $(GFF).gz

$(FAI): $(REF)
	samtools faidx $(REF)

$(CHROMSIZES): $(FAI)
	cut -f1,2 $(FAI) > $(CHROMSIZES)

# -------------------------
# HISAT2 INDEX
# -------------------------
index: $(HISAT_IDX).1.ht2

$(HISAT_IDX).1.ht2: $(REF)
	mkdir -p $(HISAT_DIR)
	hisat2-build $(REF) $(HISAT_IDX)

# -------------------------
# FASTQ
# -------------------------
fastq: $(R1) $(R2)

$(R1) $(R2):
	@if [ -z "$(SRR)" ]; then echo "ERROR: Set SRR=xxxx"; exit 1; fi
	mkdir -p $(READ_DIR)
	fastq-dump -X $(NREADS) --split-files -O $(READ_DIR) $(SRR)

# -------------------------
# FASTQC
# -------------------------
fastqc:
	@if [ -z "$(SRR)" ]; then echo "ERROR: Set SRR=xxxx"; exit 1; fi
	mkdir -p $(QC_DIR)
	fastqc -o $(QC_DIR) $(R1) $(R2)

# -------------------------
# ALIGNMENT
# -------------------------
align: $(BAM)

$(BAM): $(R1) $(R2) index
	@if [ -z "$(SRR)" ]; then echo "ERROR: Set SRR=xxxx"; exit 1; fi
	mkdir -p $(BAM_DIR)
	set -euo pipefail; \
	hisat2 -p $(THREADS) -x $(HISAT_IDX) \
		-1 $(R1) -2 $(R2) \
		2> $(BAM_DIR)/$(SRR).hisat2.log \
	| samtools sort -@ $(THREADS) -o $(BAM).tmp
	mv $(BAM).tmp $(BAM)
	samtools index $(BAM)

# -------------------------
# ALIGNMENT STATS (per sample)
# -------------------------
stats: $(BAM)
	@echo "===== Alignment stats for $(SRR) ====="
	@samtools flagstat $(BAM) | tee $(BAM_DIR)/$(SRR).flagstat.txt
	@echo "Average depth:"
	@samtools depth $(BAM) | awk '{sum+=$$3; n++} END{if(n>0) print sum/n; else print 0}' \
		| tee $(BAM_DIR)/$(SRR).mean_depth.txt

# -------------------------
# BIGWIG FILES (all samples)
# -------------------------
bigwig: $(CHROMSIZES)
	mkdir -p $(BW_DIR) $(COV_DIR)
	@for srr in $(SRRS); do \
		echo "Creating BigWig for $$srr"; \
		bedtools genomecov -ibam $(BAM_DIR)/$$srr.bam -bg > $(COV_DIR)/$$srr.bedgraph; \
		sort -k1,1 -k2,2n $(COV_DIR)/$$srr.bedgraph > $(COV_DIR)/$$srr.sorted.bedgraph; \
		bedGraphToBigWig $(COV_DIR)/$$srr.sorted.bedgraph $(CHROMSIZES) $(BW_DIR)/$$srr.bw; \
	done

# -------------------------
# COUNTS
# -------------------------
counts:
	mkdir -p $(RES_DIR)
	featureCounts \
		-T $(THREADS) \
		-p \
		-a $(GFF) \
		-t CDS \
		-g ID \
		-o $(RES_DIR)/counts.txt \
		$(BAM_DIR)/*.bam

# -------------------------
# DESEQ2
# -------------------------
de:
	Rscript run_deseq2.R $(RES_DIR)/counts.txt $(DESIGN)

# -------------------------
# CLEAN
# -------------------------
clean:
	rm -rf $(READ_DIR) $(BAM_DIR) $(QC_DIR) $(RES_DIR) $(BW_DIR) $(COV_DIR) $(FULL_DESIGN) $(DESIGN)
