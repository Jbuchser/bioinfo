# =========================
# Makefile — Assignment 13 RNA-Seq Pipeline
# =========================

# ============================ MAKEFILE SETTINGS ============================
SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# ============================= CONFIG =============================

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
BW_DIR    := bigwig
COV_DIR   := coverage
RES_DIR   := results

REF        := $(REF_DIR)/$(NAME).fa
GFF        := $(REF_DIR)/$(NAME).gff
FAI        := $(REF).fai
CHROMSIZES := $(REF_DIR)/$(NAME).chrom.sizes

SRRS := $(shell tail -n +2 $(DESIGN) | cut -d',' -f2)

# Sample-specific vars (set by `make all SRR=xxxx`)
SRR ?=
R1 := $(READ_DIR)/$(SRR)_1.fastq
R2 := $(READ_DIR)/$(SRR)_2.fastq
BAM := $(BAM_DIR)/$(SRR).bam

# ========================== PHONY TARGETS =========================

.PHONY: all all_samples design refs index fastq fastqc align stats bigwig counts clean

# -----------------------------------------------------------------
# MASTER WORKFLOW FOR ONE SAMPLE
# -----------------------------------------------------------------

all: fastq fastqc align stats

# -----------------------------------------------------------------
# RUN PIPELINE FOR ALL SAMPLES
# -----------------------------------------------------------------

all_samples: design refs index
	@echo "====== Running all samples sequentially ======"
	@for srr in $(SRRS); do \
		echo ""; \
		echo "========== Processing $$srr =========="; \
		$(MAKE) all SRR=$$srr; \
	done
	@echo "====== Running featureCounts ======"
	$(MAKE) counts

# -----------------------------------------------------------------
# DESIGN FILE
# -----------------------------------------------------------------

design:
	@echo "Fetching metadata for $(PRJ)..."
	esearch -db sra -query "$(PRJ)" | efetch -format runinfo > $(FULL_DESIGN)
	@echo "Generating $(DESIGN)..."
	awk -F, 'NR==1{ \
		for(i=1;i<=NF;i++){ \
			if($$i=="Run") r=i; \
			if($$i=="SampleName"||$$i=="Sample") s=i \
		} \
		print "sample,SRR,condition"; next \
	} \
	{ \
		sample=$$s"_"$$r; \
		srr=$$r; \
		if(sample ~ /^NT/) cond="control"; \
		else if(sample ~ /^NaP/) cond="treatment"; \
		else cond="unknown"; \
		print sample","srr","cond \
	}' $(FULL_DESIGN) > $(DESIGN)
	@echo "Design preview:"
	@head -n 7 $(DESIGN)

# -----------------------------------------------------------------
# REFERENCE DOWNLOAD + INDEX (RUN ONCE)
# -----------------------------------------------------------------

refs: $(REF) $(GFF) $(FAI) $(CHROMSIZES)

$(REF):
	mkdir -p $(REF_DIR)
	echo "Downloading FASTA..."
	curl -L \
	 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz \
	 -o $(REF).gz
	gunzip -f $(REF).gz

$(GFF):
	echo "Downloading GFF..."
	curl -L \
	 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz \
	 -o $(GFF).gz
	gunzip -f $(GFF).gz

$(FAI): $(REF)
	samtools faidx $(REF)

$(CHROMSIZES): $(FAI)
	cut -f1,2 $(FAI) > $(CHROMSIZES)

index: $(REF).bwt
$(REF).bwt: $(REF)
	bwa index $(REF)

# -----------------------------------------------------------------
# DOWNLOAD FASTQ
# -----------------------------------------------------------------

fastq:
	@if [ -z "$(SRR)" ]; then echo "ERROR: Specify SRR=xxxx"; exit 1; fi
	mkdir -p $(READ_DIR)
	echo "Downloading $(NREADS) reads for $(SRR)..."
	fastq-dump -X $(NREADS) --split-files -O $(READ_DIR) $(SRR)
	seqkit stats $(R1) $(R2)

# -----------------------------------------------------------------
# FASTQC
# -----------------------------------------------------------------

fastqc:
	@if [ -z "$(SRR)" ]; then echo "ERROR: Specify SRR=xxxx"; exit 1; fi
	mkdir -p $(QC_DIR)
	fastqc -o $(QC_DIR) $(R1) $(R2)

# -----------------------------------------------------------------
# ALIGNMENT
# -----------------------------------------------------------------

align: $(BAM)

$(BAM): $(R1) index
	mkdir -p $(BAM_DIR)
	if [ -f "$(R2)" ]; then \
		echo "Paired-end $(SRR)"; \
		bwa mem -t $(THREADS) $(REF) $(R1) $(R2) \
		| samtools sort -@ $(THREADS) -o $(BAM); \
	else \
		echo "Single-end $(SRR)"; \
		bwa mem -t $(THREADS) $(REF) $(R1) \
		| samtools sort -@ $(THREADS) -o $(BAM); \
	fi
	samtools index $(BAM)

# -----------------------------------------------------------------
# ALIGNMENT STATS
# -----------------------------------------------------------------

stats:
	@echo "Alignment statistics for $(SRR):"
	samtools flagstat $(BAM)
	@echo ""
	@echo "Average depth:"
	samtools depth $(BAM) | awk '{sum+=$$3; n++} END{if(n>0) print sum/n; else print 0}'

# -----------------------------------------------------------------
# BIGWIG FILE
# -----------------------------------------------------------------

.PHONY: bigwig
bigwig:
	@mkdir -p $(BW_DIR) $(COV_DIR)
	@for srr in $(SRRS); do \
		echo "Generating BigWig for $$srr"; \
		bedtools genomecov -ibam $(BAM_DIR)/$$srr.bam -bg > $(COV_DIR)/$$srr.bedgraph; \
		sort -k1,1 -k2,2n $(COV_DIR)/$$srr.bedgraph > $(COV_DIR)/$$srr.sorted.bedgraph; \
		bedGraphToBigWig \
			$(COV_DIR)/$$srr.sorted.bedgraph \
			$(CHROMSIZES) \
			$(BW_DIR)/$$srr.bw; \
	done

# -----------------------------------------------------------------
# FEATURE COUNTS
# -----------------------------------------------------------------

counts:

	mkdir -p ${RES_DIR}
	echo "Running featureCounts on all BAM files (paired-end)..."
	featureCounts \
		-T ${THREADS} \
		-p -B -C \
		-a ${GFF} \
		-F GFF \
		-t CDS \
		-g ID \
		-o ${RES_DIR}/counts.txt \
		${BAM_DIR}/*.bam
	echo "Counts saved to ${RES_DIR}/counts.txt"


# -----------------------------------------------------------------
# CLEANUP
# -----------------------------------------------------------------

clean:
	rm -rf $(READ_DIR) $(BAM_DIR) $(QC_DIR) $(BW_DIR) $(COV_DIR) $(RES_DIR)
