# Makefile for cancer genome project (Option 1 - Alignment analysis across different platforms)

# ----------- Makefile settings --------------

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# ------------ Region of interest --------------

# KRAS gene region (chr12:25,205,246-25,250,936)
REGION=chr12:25,205,246-25,250,936

# KRAS region without commas (for tools that don't accept commas)
REGION_NOCOMMAS := $(shell echo $(REGION) | tr -d ',')

# ------------- Reference genome (GRCh38 GIABv3) --------------

REF_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
REF     = refs/GRCh38_GIABv3.fa

# ----------------- Platform BAM URLs (HG008-T, GRCh38-GIABv3) -----------------
# Tumor Samples 

# BCM Illumina-WGS 20240313 (195x)
ILLUMINA_BAM_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam

# PacBio Revio – HiFi 20240125 (116x)
PACBIO_BAM_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Revio_20240125/HG008-T_PacBio-HiFi-Revio_20240125_116x_GRCh38-GIABv3.bam

# Northeastern ONT-UL 20241216 (54x)
ONT_BAM_URL = https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern-ONT-UL-20241216/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.bam

# ----------------------- Directories -----------------------

all: illumina pacbio ont

refs:
	mkdir -p refs

bam:
	mkdir -p bam

# ------------- Reference genome (download & index) --------------

$(REF): | refs
	wget -O $(REF).gz $(REF_URL)
	gunzip -f $(REF).gz
	samtools faidx $(REF)

# ------------- Illumina (short-read) ---------------

illumina: $(REF) | bam
	samtools view -b $(ILLUMINA_BAM_URL) $(REGION) > bam/kras_illumina.bam
	samtools index bam/kras_illumina.bam
	samtools flagstat bam/kras_illumina.bam > bam/kras_illumina.flagstat.txt
	samtools depth -r $(REGION) bam/kras_illumina.bam > bam/kras_illumina.depth.txt

# ---------------- PacBio Revio (HiFi long-read) ----------------

pacbio: $(REF) | bam
	samtools view -b $(PACBIO_BAM_URL) $(REGION) > bam/kras_pacbio.bam
	samtools index bam/kras_pacbio.bam
	samtools flagstat bam/kras_pacbio.bam > bam/kras_pacbio.flagstat.txt
	samtools depth -r $(REGION) bam/kras_pacbio.bam > bam/kras_pacbio.depth.txt

# ----------------- Oxford Nanopore Ultra-long ----------------

ont: $(REF) | bam
	samtools view -b $(ONT_BAM_URL) $(REGION) > bam/kras_ont.bam
	samtools index bam/kras_ont.bam
	samtools flagstat bam/kras_ont.bam > bam/kras_ont.flagstat.txt
	samtools depth -r $(REGION) bam/kras_ont.bam > bam/kras_ont.depth.txt

# -------------- Clean up ----------------

clean:
	rm -rf refs bam

