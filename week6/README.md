
This pipeline allows the user to download, analyze, index, and align the reference genome and FASTq data for *Staphylococcus aureus* USA300 using a Makefile. 

The entire pipeline can be executed using the command:

```bash
make all
``` 

#### Get genome

```bash
make genome
```

#### Download and analyze SRA sequencing reads using FastQ data

```bash
make analyze
```

#### Download fastq data using an SRR accession number 

```bash
make fastq
```
#### Run FastQC

```bash
make fastqc
```

#### Index the reference genome with BWA

```bash
make index
```

#### Align the reads to the reference genome

```bash
make align
```

#### Get the alignment stats including expected and observed coverage data

```bash
make stats
```

#### Clean up intermediate files 

```bash
make clean
```

# Results output:

The genome size is 2,872,769 bps, and the expected coverage is 974.02 x, or 2.8 billion bps (coverage = total bases sequenced / genome size). 

- The expected coverage is 1948.83 x

- The average coverage is 99.74%

- The observed average coverage is 4.36095

- The coverage variation was very high when looking at the full genome in IGV vs. the high coverage stemming from only sampling 10,000 reads. In IGV, there were many spaces where coverage was poor. 

![My Image](igv_snapshot.png)