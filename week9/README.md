This assignment revises and improves the automated code from Week 8. 
Inside looping was removed in favor of using GNU parallel from the outside. 

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

#### 4. Use GNU parallel to run each step in parallel (align and process each sample from design.csv):

```bash
cat design.csv | parallel --colsep , --header : -j 2 "make all SAMPLE={sample} SRR={SRR}"
```

This step uses GNU parallel to align and process (run QC, measure stats, and convert to BigWig) each run. 