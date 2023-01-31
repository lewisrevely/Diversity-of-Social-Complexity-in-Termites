MCMCglmm & figure 3 Creation
================

## packages

``` r
library(tidyverse)
library(ape)
library(geiger)
library(treeplyr)
library(MCMCglmm)
library(coda)
```

## Import data and Phylogeny

``` r
cdata <- read.csv("Traitdata.csv", row.names = 1)
species.tree<- read.tree("Species_tree.nexus")
```

### Formatting phylogeny and data for the analysis

``` r
# modify the data for this analysis
cdata$Species<- rownames(cdata)
cdata$tiplabel<-cdata$Species
cdata$genus <- sapply(strsplit(cdata$Species, split = "*_"), "[[", 1)
# remove species which do not appear in the phylogeny
cdata <- cdata[rownames(cdata) %in% species.tree$tip.label, ]

# Combine and match the tree and data
data <- make.treedata(tree = species.tree,  data = cdata, 
                           name_column = "tiplabel")

# Look at the tree
data$phy


# Look at the data
glimpse(data$dat)

# Make a new column called tiplabel with the tip labels in it
data$dat$tiplabel <- data$phy$tip.label
# Force mydata to be a data frame
mydata <- as.data.frame(data$dat)
# Save tree as mytree
mytree <- data$phy

# Remove zero length branches and replace with polytomies
mytree2 <- di2multi(mytree)
# Remove node labels 
mytree2$node.label <- NULL
# Get the inverse vcv matrix for the phylogeny
inv.phylo <- inverseA(mytree2, nodes = "TIPS", scale = FALSE)$Ainv
```

## Setting up and running the MCMCglmm

``` r
# Set up priors for MCMCglmm
gelmanprior<-list(B=list(mu=c(0,0,0,0),V=gelman.prior(~CS+HM+NC, data=mydata,
                                                      scale=1+1+pi^2/3)), R=list(V=1,fix=1),G=list(G1=list(V=diag(1)*0.1, nu=1)))




  
#universal priors  
nitt <- 1000000
thin <- 500
burnin <- 1000


mydata$WF<- as.factor(mydata$WF)
mydata$OS<-as.factor(mydata$OS)

# Fit MCMCglmm model for binary
model_mcmcglmm <- MCMCglmm(WF~CS+NC+HM, 
                           data = mydata, 
                           random = ~ tiplabel,
                           family = "threshold",
                           ginverse = list(tiplabel = inv.phylo), 
                           prior = gelmanprior,
                           nitt = nitt, thin = thin, burnin = burnin,
                           verbose = TRUE,
                           trunc = TRUE)
```

## diagnostics

``` r
# Plot model diagnostics for MCMCglmm
# For fixed effects
plot(model_mcmcglmm$Sol)  

plot(model_mcmcglmm$VCV)  

effectiveSize(model_mcmcglmm$Sol[, 1:model_mcmcglmm$Fixed$nfl, 
                                 drop = FALSE])[[1]]

# Look for autocorrelation
autocorr(model_mcmcglmm$VCV)

# Look at summary of fixed effect results
summary(model_mcmcglmm)
```

## Figure 3 creation

Using the model results from the OS and WF runs we can create figure 3

``` r
##Bring the two runs back and extract the relevant parts 
OSglmm<-readRDS("OS_MCMCglmm_results.rda")
WFglmm<-readRDS("WF_MCMCglmm_results.rda")
mod1<-OSglmm
mod2<-WFglmm

    summaryTableB<-matrix(ncol=4,nrow=4)
rownames(summaryTableB)<-c("Intercept","CS","NC","HM")
colnames(summaryTableB)<-c("Effect","2.5HPD","97.5HPD","P")

    summaryTableB1<-matrix(ncol=4,nrow=4)
rownames(summaryTableB1)<-c("Intercept","CS","NC","HM")
colnames(summaryTableB1)<-c("Effect","2.5HPD","97.5HPD","P")

##For the Obligate sterility run create a summary table of the results 
summaryTableB[,1]<-summary(mod1)$solutions[,1]
summaryTableB[,2]<-summary(mod1)$solutions[,2]
summaryTableB[,3]<-summary(mod1)$solutions[,3]
summaryTableB[,4]<-summary(mod1)$solutions[,5]
summaryTableB<-summaryTableB[-1,]   
summary<-as.data.frame(summaryTableB)
summary$Traits<-rownames(summary)
summary$response<-rep('Obligate Sterility',3)

##For the Foraging run create a summary table of the results
summaryTableB1[,1]<-summary(mod2)$solutions[,1]
summaryTableB1[,2]<-summary(mod2)$solutions[,2]
summaryTableB1[,3]<-summary(mod2)$solutions[,3]
summaryTableB1[,4]<-summary(mod2)$solutions[,5]
summaryTableB1<-summaryTableB1[-1,]
summary1<-as.data.frame(summaryTableB1)
summary1$Traits<-rownames(summary)
summary1$response<-rep('Foraging',3)

##Merge together so that we have one dataset of results with effect size and confidence limits
summary2<-rbind(summary1,summary)

##Create the plot 
p<-summary2%>%
  mutate(Low= `2.5HPD` ,High=`97.5HPD`) %>%
  ggplot(aes(x= Traits , y=Effect, color= Traits, shape = response)) + 
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbar(aes(ymin=Low,ymax=High),width=0.1,position = position_dodge(0.5), size= 1)+
  ggtitle("MCMCglmm")+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_classic()
p2<-  p + scale_color_manual(values=c("#D85859", "#FCB76F", "#E88766"))


ggsave(p2, filename='CL_overlap_plot.pdf')
```
