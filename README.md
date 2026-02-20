# PhyloIntroBiochem

Hands-on introduction to phylogenetics (for USU Biochem students)

This exercise uses NCBI for sequence retrieval plus R for alignment and maximum-likelihood phylogenetics.

Mini-lecture slides

## Data

This assumes you have the DNA sequence from a gene and organism of interest in a fasta file, and that you want to estimate of phylogeny centered on this gene and species. Such a phylogeny could include the same gene in additional species or populations or other members of the same gene family or both. In this example, I am starting with an opsin gene from *Drosophila melanogaster* (for arbitrary reasons). 

[DmelOpsin.fasta](https://github.com/user-attachments/files/25432149/sequence.txt)

We begin by obtaining additional DNA sequences from NCBI.

1. Navigate to the NCBI nucleotide blast page: [blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome). (You can use blastp if you have amino acide data)
2. Enter (paste) the contents of your fasta file under **Enter Query Sequence**.
3. Choose the **Core nucelotide database** and set the Organism field to limit your search. I used Drosophilidae (taxid: 7214).
4. Select **Uncultured/environmental sample sequences** to reduce noise
5. Choose megablast, discontiguous megablast or blastn 

Select a subset of sequences to use. Download the sequences as a fasta file. You can edit the headers in your fasta file to have more managable names downstream. Here is an example: [dseqs.fasta](https://github.com/user-attachments/files/25432277/dseqs.txt)


Keep sequences that:

- match the same gene name in the annotation (or obvious homolog)

- are similar length (avoid very short fragments if possible)

- come from species of interest; you might or might not want multiple representatives of the same species

Avoid:

- "partial" sequences that are extremely short

- records with unclear gene identity

NCBI has simple tools for viewing a sequence alignment (MSA viewer) and making a tree from your results. You can check these out but we will not use them as our final product.

## Analysis

Install the following R packages.

```r
install.packages("ape")
install.packages("phangorn")
install.packages("BiocManager")
BiocManager::install("DECIPHER")
```
