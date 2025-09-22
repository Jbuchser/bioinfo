# Week 2 Assignment 

1. The organism I am working with is *Felis catus* (the domestic cat).

- I tried to use the bio fetch command but was unsuccessful because of a problem with the GFF file already being unzipped. 
- Instead, I manually downloaded the public GFF3 data locally, unzipped it, and ensured it was present on my system:

```bash
wget ftp://ftp.ensembl.org/pub/release-115/gff3/felis_catus/Felis_catus.F.catus_Fca126_mat1.0.115.gff3.gz

gunzip Felis_catus.F.catus_Fca126_mat1.0.115.gff3.gz

ls -lh Felis_catus.F.catus_Fca126_mat1.0.115.gff3
```
- I analyzed the first 20 lines of the data to get more information about the organism:

```bash 
head -20 Felis_catus.F.catus_Fca126_mat1.0.115.gff3
##gff-version 3
##sequence-region   A1 1 239367248
##sequence-region   A2 1 169388855
##sequence-region   A3 1 140443288
##sequence-region   B1 1 205367284
##sequence-region   B2 1 151959158
##sequence-region   B3 1 148491486
##sequence-region   B4 1 142168536
##sequence-region   C1 1 221611373
##sequence-region   C2 1 158578461
##sequence-region   D1 1 115366950
##sequence-region   D2 1 88083857
##sequence-region   D3 1 94435393
##sequence-region   D4 1 95154158
##sequence-region   E1 1 61876196
##sequence-region   E2 1 61988844
##sequence-region   E3 1 41437797
##sequence-region   F1 1 69239673
##sequence-region   F2 1 83466477
##sequence-region   JAFEKA010000039.1 1 89503
```
- I noticed that this downloaded data was missing full annotations so I again manually downloaded the data (This time I took the genomic.gff file directly from NCBI), and again asked for the head (this time only -5). 
- I also asked to remove the comment annotations and any blank lines: 

```bash
grep -v "^#" genomic.gff | grep -v "^\s*$" | head -n 5
```

2. This code returned that the domestic cat has 70 sequence regions, but this is not expected (the domestic cat has 18 paired autosomes and 1 pair of sex chromosomes, or 38 total chromosomes):

```bash
grep "##sequence-region" Felis_catus.F.catus_Fca126_mat1.0.115.gff3 | wc -l
      70
```
- Next, I asked for all unique sequence regions instead:

```bash
awk '{if($1 !~ /^#/) print $1}' Felis_catus.F.catus_Fca126_mat1.0.115.gff3 | sort | uniq
```
This returned all 19 chromosome pairs and extra scaffolding data, so I filtered out the scaffolding data by saving only the chromosomes in their own file called "main_chroms.txt."

```bash
echo -e "A1\nA2\nA3\nB1\nB2\nB3\nB4\nC1\nC2\nD1\nD2\nD3\nD4\nE1\nE2\nE3\nF1\nF2\nX" > main_chroms.txt
```
3. To determine the total features, I asked for all features (including from the scaffolding data) which returned 2,301,949 features:

```bash
grep -v "^#" Felis_catus.F.catus_Fca126_mat1.0.115.gff3 | wc -l
 2301949
```
vs. asking for the total features in the main chromosomes only, which returned 2298398 features. 

```bash
awk '$1 ~ /^A|^B|^C|^D|^E|^F|^X/' Felis_catus.F.catus_Fca126_mat1.0.115.gff3 | wc -l
 2298398
```

4. I asked for the number of genes using data from the main_chroms.txt file, which returned 19,165 genes:

```bash
awk '$3=="gene"' Felis_catus.F.catus_Fca126_mat1.0.115.gff3 | grep -F -f main_chroms.txt | wc -l
   19165
```

5. I asked for all unique sequence features (listed below):

```bash
awk '{if($1 !~ /^#/) print $3}' Felis_catus.F.catus_Fca126_mat1.0.115.gff3 | sort | uniq
CDS
C_gene_segment
J_gene_segment
V_gene_segment
Y_RNA
biological_region
exon
five_prime_UTR
gene
lnc_RNA
mRNA
miRNA
ncRNA_gene
pseudogene
pseudogenic_transcript
rRNA
region
scRNA
snRNA
snoRNA
three_prime_UTR
transcript
```
> **snoRNA** : snoRNA are small nucleolar RNAs found in the eukaryotic nucleolus and assist in the modification of other RNAs.

- There are 638 snoRNAs annotated for *Felis catus*:

```bash
awk '$3=="snoRNA"' Felis_catus.F.catus_Fca126_mat1.0.115.gff3 | wc -l

     638
```
6. I identified the top ten annotated features:

```bash 
grep -v "^#" genomic.gff | awk '{print $3}' | sort | uniq -c | sort -nr | head -n 10
1194611 exon
964710 CDS
71463 mRNA
35687 gene
20463 lnc_RNA
6354 cDNA_match
5520 transcript
3553 pseudogene
1124 snRNA
 763 tRNA
```

7. This appears to be a complete and fully-annotated organism based on the analyzed GFF data. 

8. Additional interesting data I was able to locate is that the *Felis catus* genome uploaded to NCBI was from a female cat (based on the GFF data sex=female). I was also able to ID that the genomic data came from fibroblast tissue. 
