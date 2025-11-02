This assignment reuses the automated code from Week 9 with the addition of variant calling using BCFtools (Samtools).

The data flow for this pipeline is:

1. Collect the FASTQ reads and reference genome
2. Align FASTQ reads with the reference genome
3. Use alignment software to produce a BAM file
4. Call variants by comparing the aligned BAM file to the reference genome
5. Produce VCF files
6. Merge VCF files 

Note: The reference genome used in this workflow is NC_007793.1 *Staphylococcus aureus* strain USA_300

### Steps to the individual Makefile

#### 1. Create a minimal design.csv (sample name and SRR number only) extracted from the full BioProject dataset

```bash
make design
```
This step downloads the metadata from a BioProject number and then extracts a design.csv file limited to multiple sample IDs and SRR accession numbers. 

#### 3. Download the genome

```bash
make genome
```
This step downloads the reference genome from the BioProject SRR data and saves it with a user-friendly name. 

#### 4. Index the genome 

```bash
make index
```
This step indexes the genome for the next steps (which are done in parallel). 

#### 4. Collect SRR reads

```bash
make reads
```

#### 5. Run FastQC on the reads

```bash
make fastqc
```

#### 6. Align the reads with the reference genome

```bash
make align
```

#### 7. Run basic statistics on the aligned reads

```bash
make stats
```

#### 8. Convert BAM to BigWig format via temporary BedGraph format

```bash
make bigwig
```

#### 9. Call variants using BCFtools

```bash
make call_variants
```

#### **GNU parallel can be used to run each step in parallel (align and process each sample from design.csv):**

```bash
cat design.csv | parallel --colsep , --header : -j 2 "make all SAMPLE={sample} SRR={SRR}"
```

This step uses GNU parallel and the *make all* command to align, process, and call variants for each sample (make all: reads, FastQC, align, stats, bigwig, call variants). 

#### After running multiple samples in parallel, use this command to merge each VCF file. The output of this command provides statistics on variants from all samples. 

```bash
make merge_vcf
```

Output of the 6 merged VCF files:

**Number of samples:** 6 <br>
**Total variant sites:** 1399 <br>
**Total SNPs:** 1341 <br>
**Total indels:** 58 <br>
**Total multiallelic sites:** 7 <br>
**Total transitions/transversions:** 2.19 <br>
**Singletons (variants per 1 sample):** 331 <br>
**Depth distribution:** 1-4x coverage/site <br>

The observed low coverage is likely due to downloading only a subset of reads (N=10000). 