## 1. Downloading the genome
1. The accession number for the genome (methycillin-resistant *Staphylococcus aureus* paper) is NC_007793.

2. I downloaded and renamed the genomic data for NC_007793:

FASTA:

```bash
bio fetch NC_007793.1 -format fasta > mrsa.fna
```
GFF:
```bash
bio fetch NC_007793.1 -format gff > mrsa.gff
```

## 2. Visualizing the genome 

1. The genome is 2872769 bp long:

```bash
awk '!/^>/ {sum += length($0)} END {print sum}' mrsa.fna
```
2. Sequence features:

```bash
$ grep -v '^#' mrsa.gff | awk -F'\t' 'NF==9 {print $3}' | sort | uniq -c | sort -nr

2713 CDS
2689 gene
  78 pseudogene
  72 exon
  52 tRNA
  15 rRNA
  13 riboswitch
   3 sequence_feature
   3 binding_site
   1 tmRNA
   1 region
   1 ncRNA
   1 SRP_RNA
   1 RNase_P_RNA
```

3. SAUSA300_RS072 is the longest gene with 31266 bps. <br>
RNA-seq-based transcriptome analysis states SAUSA300_RS072 is also annotated as *ilvA* and is a threonine dehydratase (Im *et al*., 2022)

```bash
grep -v '^#' mrsa.gff | awk -F'\t' '$3=="gene" {len=$5-$4+1; if(len>max){max=len; id=$9}} END {print "Longest gene:", id, "Length:", max}'
```

4. I browsed the genome using IGV and picked another gene called *murT* encoding the murT subunit of the murT-gatD complex. This complex is essential for proper peptidoglycan synthesis and cell wall integrity (particularly important for Gram-positive *S. aureus*) (Gonçalves *et al*., 2019).

5. I converted the GFF file to BED to merge overlapping CDS to get the coding density of the genome (86.8056%). <br> 
-   I needed to do this because just counting the CDS gave me 10x the size of the actual genome (29513151 bps)

```bash
awk 'BEGIN{genome=2880000; coding=2500000; print "Coding density:", coding/genome*100 "%"}'

bedtools merge -i cds.bed > cds_merged.bed

awk '{sum+=$3-$2} END {print "Total coding bp:", sum}' cds_merged.bed

awk 'BEGIN{genome=2880000; coding=2500000; print "Coding density:", coding/genome*100 "%"}'
```
86.8% coding density matches what can be observed with IGV, with minimal intergenic space/tightly packed genes as expected for *S. aureus*. 

## 3. Alternative genome of interest

CP003029 is the accession number for a vancomycin resistant strain of *S. aureus*. <br>
This strain could be used to compare resistance-associated genes with MRSA strains (like NC_007793.1). <br>
Vancomycin is a last-resort antibiotic used to treat severe/systemic MRSA infections, and increasing vancomycin resistance is a major concern. 