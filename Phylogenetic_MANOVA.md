Phylogenetic Manova
================

## packages

``` r
library(ape)
library(phytools)
library(vegan)
```

## Import data and Phylogeny

``` r
cdata <- read.csv("Traitdata.csv", row.names = 1)
species.tree<- read.tree("Species_tree.nexus")
```

## Phylogenetic Manova set up

``` r
# add column called species to cdata for later use
cdata$Species<- rownames(cdata)
# (excluding species name, genus name, and what ever you're wanting to compare eg OS)
data.use <- cdata[, !names(cdata) %in% c("OS","WF","Species")]
head(data.use) # correct. 


# Check that all trait species are in the phylogeny
rownames(data.use) %in% species.tree$tip.label # Yes! Great. 



head(data.use) # All looks good. 
summary(data.use) # Looks good. 


#### Use PCoA to summarise the data. ####
# 1. Create a distance matrix
head(data.use)

# scale and centre the traits (z transformation)
scaled.data <- data.frame(scale(data.use))


# Create the distance matrix
dist.mat <- dist(scaled.data)

# 2. Perform PCoA

pcoa.result <- pcoa(dist.mat)

# 3. Extract PCoA traits
pcoa.traits <- pcoa.result$vectors

# 4. Check variation explained by the PCoA
pcoa.result$values # 1st PCoA axis explains 17% of the original data. 
# 2nd axis explains 16%. Etc. Look for the "Relative eigenvalue" column. 
```

## Perform Multivariate MANOVA WITH phylogenetic correction

``` r
phylo.pca <- phyl.pca(tree = species.tree, Y = data.use)
phylo.axes <- phylo.pca$S
ppca.summ <- summary(phylo.pca)
ppca.imp <- cumsum(ppca.summ$importance[2, ])

PhyloMANOVA <- NULL
for(i in 1:ncol(phylo.axes)){
  trait.dis <- dist(phylo.axes[, 1:i])
  
  FG <- cdata[cdata$Species %in% rownames(phylo.axes), "OS"] ####change here to explanatory variable of interest
  
  manova.result <- adonis2(trait.dis ~ FG) # this is the MANOVA
  
  r2 <- manova.result$R2[1]
  temp <- data.frame(r2 = r2, i = i)
  PhyloMANOVA <- rbind(PhyloMANOVA, temp)
  print(i)
  
}

manova.result
```
