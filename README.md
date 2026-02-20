# PhyloIntroBiochem

Hands-on introduction to phylogenetics (for USU Biochem students)

This exercise uses NCBI for sequence retrieval plus R for alignment and maximum-likelihood phylogenetics.

Mini-lecture slides

## Data
This assumes that you have DNA sequence data for a gene and organism of interest in a FASTA file, and that you want to estimate a phylogeny centered on this gene and species. Such a phylogeny could include the same gene from additional species or populations, other members of the same gene family, or both. In this example, I begin with an opsin gene from *Drosophila melanogaster* (chosen arbitrarily).

[DmelOpsin.fasta](https://github.com/user-attachments/files/25432149/sequence.txt)

Begin by obtaining additional DNA sequences from NCBI.

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

We will use R to align the DNA sequences and estimate the phylogeny. Install the following R packages.

```r
install.packages("ape")
install.packages("phangorn")
install.packages("BiocManager")
BiocManager::install("DECIPHER")
```

Now lets align the sequences. We will use the R package `DECIPHER`. [MUSCLE](https://www.ebi.ac.uk/jdispatcher/msa/muscle) and [ClustalW](https://www.genome.jp/tools-bin/clustalw) are other options.

```r
## load libraries
library("DECIPHER")
library("ape")
library("phangorn")

## read in the sequence data
fa<-readDNAStringSet("Downloads/dseqs.txt")
fa
width(fa)

## align the sequences
aln<-AlignSeqs(fa)
aln_mat <- as.matrix(aln)
phydat <- phyDat(aln_mat, type = "DNA")
phydat
```

We will start with a simple neighbor-joining tree. We will use this a starting point to optimize a maximum likelihood tree.

```r
## calcualte the distance matrix, uses JC69 by default
dm <- dist.ml(phydat)
## make the NJ tree
tree_nj <- NJ(dm)

## ML fit, optimize under GTR
## we could use other models or try model selection
fit <- pml(tree_nj, data = phydat)
fit_gtr <- optim.pml(fit,
  model = "GTR",
  rearrangement = "stochastic",
  control = pml.control(trace = 0)
)
## plot the tree
## this should have approximate SH-like (aLRT) branch support
plot(fit_gtr)
```
You can also generate bootstrap support values

```r
bs <- bootstrap.pml(fit_gtr, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

## pot with bootstrap labels (>= 50%)
tree_bs <- plotBS(fit_gtr$tree, bs, p = 50, type = "phylogram")
plot(tree_bs)
```
We can compare alternative models.
```r
## model testing
modelTest(phydat)
```

And we can reroot the tree
```r
root(tree, outgroup="ACCESSION", resolve.root=TRUE)
```
