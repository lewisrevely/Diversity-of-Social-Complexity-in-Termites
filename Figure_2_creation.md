Figure 2 Creation
================

## packages

``` r
library(ggtree)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(viridis)
```

## Import data and Phylogeny

``` r
df <- read.csv("Traitdata.csv", row.names = 1)
species.tree<- read.tree("Species_tree.nexus")
```

## Create circular phylogeny without tiplabels

``` r
circ <- ggtree(species.tree, layout = "rectangular")
```

## Modify trait data

We need to scale the discrete data into scale between 0 and 1 and
factorise the binary data

``` r
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df$CS<-range01(df$CS)
df$HM<-range01(df$HM)
df$NC<-range01(df$NC)

df$OS<-as.factor(df$OS)
df$WF<-as.factor(df$WF)
```

## Create the heatmap from Figure 2

``` r
p1 <- gheatmap(circ, df, offset=0, width=.25,
               colnames_angle=90, colnames_offset_y = .01,font.size=2)
```

    ## Warning: attributes are not identical across measure variables;
    ## they will be dropped

``` r
## Change the colours to be better in black and white if need be 

p2<-p1 + scale_fill_viridis(discrete = TRUE, name = "Complexity level", labels =c('Low','','','','','','','','High'), direction = -1, option= "C")
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

``` r
plot(p2)
```

![](Figure_2_creation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

I then put this into Inkscape for adding the photographs and modifying
the legend and labels
