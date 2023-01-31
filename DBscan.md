DBscan
================

## packages

``` r
library(fpc)
library(dbscan)
library(KneeArrower)
```

## Import data

``` r
dat <- read.csv("Morphometric_data.csv")
```

## Clean up the dataset

``` r
str(dat)
dat<-dat[,-1]
dat$HF <- as.numeric(dat$HF)
dat$FT <- as.numeric(dat$FT)
dat$HW <- as.numeric(dat$HW)
dat <-na.omit(dat)
```

## Looking at an individual species and caste and running the clustering analysis

``` r
##subsetting to the species and caste of interest
subset<- subset(dat, id == "Macrotermes gilvus Worker")

##Manually finding the best eps using the position at which the curve changes the most
knndist<- as.matrix(kNNdist(subset[,-4], k = 4, all = FALSE))

## the $y should provide the best cutoff and therefore the best eps value
cuttoff<- findCutoff(row(knndist),knndist, method = "curvature")
cuttoff

##best to visualise the curve sometimes as $y can be a bit off
kNNdistplot(subset[,-4], k =  4)+
abline(h =cuttoff$y, lty = 2)
```

![](DBscan_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Fitting DBScan clustering Model
# to training dataset
set.seed(220) # Setting seed
Dbscan_cl <- dbscan(subset[,-4], eps = 0.25, MinPts = 4)
Dbscan_cl$cluster


# Checking cluster
Dbscan_cl$cluster

# 
table(Dbscan_cl$cluster, subset$id)

# Plotting Cluster
plot(subset[,-4])
```

![](DBscan_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
plot(subset[,-4], col=Dbscan_cl$cluster+1, main="DBSCAN")
```

![](DBscan_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
#get a single number of clusters 
max(Dbscan_cl$cluster)
```

Now run a for loop over all the species and castes to get automated
morph number. This needs to be manually checked for obvious errors which
likely will be due to eps choice issues. In that case do manually. for
extra care it should be done manually so that you can use the above
diagnostic plots to sense check

``` r
DBscan_results<- NA
eps_results<- NA
for (i in unique(dat$id)){ 
  try({
    data <- dat[dat$id == i,]
    knndist<- as.matrix(kNNdist(data[,-4], k = 4, all = FALSE))
    eps<-findCutoff(row(knndist),knndist, method = "curvature")
    Dbscan_cl <- dbscan(data[,-4], eps = eps$y, MinPts = 4)
    DBscan_clust <- as.data.frame(max(Dbscan_cl$cluster))
    eps_clust<- as.data.frame(eps$y)
    DBscan_clust$id<- i
    eps_clust$id<- i
    DBscan_results<-rbind(DBscan_results, DBscan_clust)
    eps_results<-rbind(eps_results, eps_clust)
  })
}

DBscan_results
```

Once this has been produced, we compared with our picutres and existing
literature to sense check our results then compiled morph number of both
soldier and worker into a single trait called ‘helper polymorphism’
within excel and brought it together with other traits into a single
file “Traitdata.csv”
