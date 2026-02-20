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

bs <- bootstrap.pml(fit_gtr, bs = 100, optNni = TRUE, control = pml.control(trace = 0))

## pot with bootstrap labels (>= 50%)
tree_bs <- plotBS(fit_gtr$tree, bs, p = 50, type = "phylogram")
plot(tree_bs,show.node.label=TRUE)

## model testing
mo<-modelTest(phydat)
## extract the best model
as.pml(mo)

rt<-root(fit_gtr$tree, outgroup="Drosophila affinis", resolve.root=TRUE)
plot(rt)
