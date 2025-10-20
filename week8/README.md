This assignment reuses the Makefile from Week 7 (with new changes for parallel automation), and also includes a design.csv file listing minimal SRR data from BioProject PRJNA887926. 

Additionally, the Makefile is reusable with different BioProject and SRR numbers if needed. 

#### 1. Run the full pipeline

```bash
make all 
```
or, run each command individually using the make commands below:

#### 2. Create a minimal design.csv (sample name and SRR number only) extracted from the full BioProject dataset

```bash
make design
```
#### 3. Download the genome and save it with a user-friendly name (MRSA-USA300 in this example)

```bash
make genome
```

#### 4. Download the reads in parallel for each SRR number

```bash
make reads
```

#### 5. Run fastqc on the reads from each SRR number

```bash
make fastqc
```
#### 6. Index the reference genome

```bash
make index
```

#### 7. Align and convert reads to BAM format

```bash
make align
```
#### 8. Run alignment stats on each BAM read

```bash 
make stats
```

#### 9. Convert BAM to BedGraph (temporary), then convert BedGraph to BigWig 

```bash
make bigwig
```
