MCMCglmm & figure 3 Creation
================

## Load Required Packages
Load libraries necessary for data manipulation, phylogenetic analysis, and MCMCglmm.

``` r
library(tidyverse)
library(ape)
library(geiger)
library(treeplyr)
library(MCMCglmm)
library(coda)
```

## Import Data and Phylogeny
Import trait data and the species phylogenetic tree.

``` r
cdata <- read.csv("Traitdata.csv", row.names = 1)
species.tree<- read.tree("Species_tree.nexus")
```

### Formatting Phylogeny and Data for Analysis
Prepare data and phylogeny for compatibility and analysis.

``` r
# Add species and genus information to the data
cdata$Species <- rownames(cdata)
cdata$tiplabel <- cdata$Species
cdata$genus <- sapply(strsplit(cdata$Species, split = "_"), "[[", 1)

# Exclude species not present in the phylogeny
cdata <- cdata[rownames(cdata) %in% species.tree$tip.label, ]

# Combine tree and data ensuring they match
data <- make.treedata(tree = species.tree, data = cdata, name_column = "tiplabel")

# Visual inspection of the tree and data
data$phy
glimpse(data$dat)

# Prepare the data frame for analysis and clean up the phylogenetic tree
mydata <- as.data.frame(data$dat)
mydata$tiplabel <- data$phy$tip.label
mytree <- data$phy
mytree2 <- di2multi(mytree)  # Convert zero length branches to polytomies
mytree2$node.label <- NULL   # Remove node labels
inv.phylo <- inverseA(mytree2, nodes = "TIPS", scale = FALSE)$Ainv  # Inverse vcv matrix for phylogeny
```

## Setup and Run MCMCglmm

Configure and execute the MCMCglmm model for analysis.

``` r
# Define priors for MCMCglmm - Gelman prior for fixed effects and inverse Wishart for random effects (variance set to 1 since it is binary data)
gelmanprior<-list(B=list(mu=c(0,0,0,0),V=gelman.prior(~CS+HM+NC, data=mydata,
                                                      scale=1+1+pi^2/3)), R=list(V=1,fix=1),G=list(G1=list(V=diag(1)*0.1, nu=1)))

# Set parameters for MCMCglmm 
nitt <- 100000
thin <- 50
burnin <- 5000


# Convert factors for the model
mydata$WF<- as.factor(mydata$WF)
mydata$OS<-as.factor(mydata$OS)


# Fit the MCMCglmm model  
model_mcmcglmm <- MCMCglmm(OS~CS+NC+HM, 
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
Assess the model's performance and diagnostics.

``` r
# Plot model diagnostics for fixed effects
plot(model_mcmcglmm$Sol)  

plot(model_mcmcglmm$VCV)  

# Evaluate effective size and check for autocorrelation
effectiveSize(model_mcmcglmm$Sol[, 1:model_mcmcglmm$Fixed$nfl, 
                                 drop = FALSE])[[1]]
autocorr(model_mcmcglmm$VCV)

# Summary of fixed effect results
summary(model_mcmcglmm)

# Convergence diagnostics using Gelman plots
listmc <- mcmc.list(model_mcmcglmm$Sol,model_mcmcglmm1$Sol,model_mcmcglmm2$Sol)
gelman.plot(listmc)
gelman.diag(listmc)
summary(listmc)

# Compile model chains for figure creation
dfmc <- rbind.data.frame(model_mcmcglmm$Sol,model_mcmcglmm1$Sol,model_mcmcglmm2$Sol)

sumtable<-apply(dfmc,2,quantile, probs=c(0.05,0.95))
mean<-summary(listmc)
mean<-mean$statistics[,1]
effect<-effectiveSize(listmc)

WF_chain_summary<-t(rbind(mean,sumtable,effect))
```

## Create Figure 3

Generate a visual representation of the model results.

``` r
# Combine model results for figure creation
OS_chain_summary<-readRDS("OS_MCMCglmm_results.rda")
WF_chain_summary<-readRDS("WF_MCMCglmm_results.rda")
mod1<-OS_chain_summary
mod2<-WF_chain_summary


# Prepare summary data for plotting
summary<-mod1[-1,]
summary<-as.data.frame(summary)
summary$Traits<-rownames(summary)
summary$response<-rep('Obligate Sterility',3)

summary1<-mod2[-1,]
summary1<-as.data.frame(summary1)
summary1$Traits<-rownames(summary1)
summary1$response<-rep('Foraging',3)

summary2<-rbind(summary,summary1)

# Plot the results 
p<-summary2%>%
  mutate(Low= `5%` ,High=`95%`) %>%
  ggplot(aes(x= Traits , y=mean, shape = response)) + 
  geom_point(position = position_dodge(0.5), size = 3) +
  geom_errorbar(aes(ymin=Low,ymax=High),width=0.1,position = position_dodge(0.5), size= 1)+
  ylab("Effect")+
  scale_x_discrete(labels = c("Colony Size", "Helper Polymorphism", "Nest Complexity"))+
  geom_hline(yintercept = 0, linetype="dotted")+
  theme_classic()

# save the plot
ggsave(p, filename='CL_chained_overlap_plot.pdf')
```
