Tree Creation
================

\## packages

``` r
library(RRphylo)
library(ape)
library(ggtree)
library(geiger)
```

## Import and tidy the phylogeny

We use the data file and the phylogeny file which was produced using
pyphawld to refine and tidy the phylogeny for our work

``` r
setwd("~/OneDrive - University College London/PhD/Traits/")
cdata <- read.csv("Traitdata.csv", row.names = 1)
  tree <- read.tree(file = "pyphlawdtree.nexus")

#adding column for genus names and species as rownames
cdata$Species<- rownames(cdata) # set species names as row names too. 
cdata$genus <- sapply(strsplit(cdata$Species, split = "*_"), "[[", 1)

####quick fix for masto problem
tree$tip.label[tree$tip.label == "36984"] <- "Mastotermes_darwiniensis"


# Which genera are on the tree? 
treedata <- data.frame(tree_taxa = tree$tip.label)
treedata$genus <- sapply(strsplit(treedata$tree_taxa, split = "*_"), "[[", 1)


table(cdata$genus %in% unique(treedata$genus)) 

# Which tree labels do we need to keep? 
treedata2 <- treedata[treedata$genus %in% cdata$genus, ]

#this is for omitted genus'
treedata_omit<- cdata[! cdata$genus %in% treedata$genus , ]
treedata_omit[unique(treedata_omit$genus), ]
treedata_omit <- treedata_omit[!duplicated(treedata_omit$genus), ]

# Keep only one representative species per genus
treedata2[unique(treedata2$genus), ]
treedata2 <- treedata2[!duplicated(treedata2$genus), ]
nrow(treedata2)  
tips2keep <- treedata2$tree_taxa
tips2drop <- treedata$tree_taxa[!treedata$tree_taxa %in% tips2keep]

# Trim the tree
tree.trim <- drop.tip(tree, tip = tips2drop, trim.internal = TRUE)

# Alter tree tip labels
tree.trim$tip.label <- sapply(strsplit(tree.trim$tip.label, 
                                       split = "*_"), "[[", 1)

# How many species of each genus do we need to add? 
table(treedata2$genus)
table(table(treedata2$genus) == 1) # Each genus only has a single species!

# This makes our job easier. We can just rename the genera on the phylogeny:
species.tree <- tree.trim
species.tree$tip.label <- 
  treedata2[treedata2$genus == species.tree$tip.label, "tree_taxa"]


#change tip genus names so that they're identical to our species

species.tree$tip.label[species.tree$tip.label == "Stolotermes_brunneicornis"] <- "Stolotermes_victoriensis"
species.tree$tip.label[species.tree$tip.label == "Zootermopsis_laticeps"] <- "Zootermopsis_nevadensis"
species.tree$tip.label[species.tree$tip.label == "Kalotermes_brouni"] <- "Kalotermes_flavicollis"
species.tree$tip.label[species.tree$tip.label == "Neotermes_meruensis"] <- "Neotermes_sanctaecrucis"
species.tree$tip.label[species.tree$tip.label == "Rugitermes_panamae"] <- "Rugitermes_rugosus"
species.tree$tip.label[species.tree$tip.label == "Incisitermes_snyderi"] <- "Incisitermes_milleri"
species.tree$tip.label[species.tree$tip.label == "Bifiditermes_improbus"] <- "Bifiditermes_durbanensis"
species.tree$tip.label[species.tree$tip.label == "Glyptotermes_fuscus"] <- "Glyptotermes_brevicornis"
species.tree$tip.label[species.tree$tip.label == "Schedorhinotermes_actuosus"] <- "Schedorhinotermes_medioobscurus"
species.tree$tip.label[species.tree$tip.label == "Parrhinotermes_queenslandicus"] <- "Parrhinotermes_aequalis"
species.tree$tip.label[species.tree$tip.label == "Psammotermes_voeltzkowi"] <- "Psammotermes_hybostoma"
species.tree$tip.label[species.tree$tip.label == "Prorhinotermes_simplex"] <- "Prorhinotermes_inopinatus"
species.tree$tip.label[species.tree$tip.label == "Microtermes_pakistanicus"] <- "Microtermes_magnocellus"
species.tree$tip.label[species.tree$tip.label == "Allodontermes_tenax"] <- "Allodontermes_morogorensis"
species.tree$tip.label[species.tree$tip.label == "Protermes_minutus"] <- "Protermes_prorepens"
species.tree$tip.label[species.tree$tip.label == "Macrotermes_jinghongensis"] <- "Macrotermes_michaelseni"
species.tree$tip.label[species.tree$tip.label == "Dicuspiditermes_makhamensis"] <- "Dicuspiditermes_nemorosus"
species.tree$tip.label[species.tree$tip.label == "Cristatitermes_tutulatus"] <- "Cristatitermes_Froggatti"
species.tree$tip.label[species.tree$tip.label == "Xylochomitermes_occidualis"] <- "Xylochomitermes_melvillensis"
species.tree$tip.label[species.tree$tip.label == "Neocapritermes_taracua"] <- "Neocapritermes_parvus"
species.tree$tip.label[species.tree$tip.label == "Cornitermes_acignathus"] <- "Cornitermes_bequaerti"
species.tree$tip.label[species.tree$tip.label == "Amitermes_foreli"] <- "Amitermes_meridionalis"
species.tree$tip.label[species.tree$tip.label == "Globitermes_brachycerastes"] <- "Globitermes_sulphureus"
species.tree$tip.label[species.tree$tip.label == "Fulleritermes_coatoni"] <- "Fulleritermes_tenebricus"
species.tree$tip.label[species.tree$tip.label == "Nasutitermes_glabritergus"] <- "Nasutitermes_arborum"
species.tree$tip.label[species.tree$tip.label == "Velocitermes_velox"] <- "Velocitermes_heteropterus"
species.tree$tip.label[species.tree$tip.label == "Trinervitermes_oeconomus"] <- "Trinervitermes_bettonianus"
species.tree$tip.label[species.tree$tip.label == "Reticulitermes_yaeyamanus"] <- "Reticulitermes_lucifugus"
species.tree$tip.label[species.tree$tip.label == "Postelectrotermes_howa"] <- "Postelectrotermes_militaris"
species.tree$tip.label[species.tree$tip.label == "Procryptotermes_leewardensis"] <- "Procryptotermes_falcifer"
species.tree$tip.label[species.tree$tip.label == "Proneotermes_macondianus"] <- "Proneotermes_perezi"
species.tree$tip.label[species.tree$tip.label == "Epicalotermes_mkuzii"] <- "Epicalotermes_kempae"
species.tree$tip.label[species.tree$tip.label == "Calcaritermes_temnocephalus"] <- "Calcaritermes_emarginicollis"
species.tree$tip.label[species.tree$tip.label == "Alyscotermes_trestus"] <- "Alyscotermes_kilimandjaricus"
species.tree$tip.label[species.tree$tip.label == "Adaiphrotermes_choanensis"] <- "Adaiphrotermes_cuniculator"
species.tree$tip.label[species.tree$tip.label == "Ruptitermes_arboreus"] <- "Ruptitermes_reconditus"
species.tree$tip.label[species.tree$tip.label == "Grigiotermes_hageni"] <- "Grigiotermes_metoecus"
species.tree$tip.label[species.tree$tip.label == "Armitermes_peruanus"] <- "Armitermes_cerradoensis"
species.tree$tip.label[species.tree$tip.label == "Syntermes_molestus"] <- "Syntermes_grandis"
species.tree$tip.label[species.tree$tip.label == "Araujotermes_caissara"] <- "Araujotermes_parvellus"
species.tree$tip.label[species.tree$tip.label == "Bulbitermes_laticephalus"] <- "Bulbitermes_singaporiensis"
species.tree$tip.label[species.tree$tip.label == "Hospitalitermes_bicolor"] <- "Hospitalitermes_hospitalis"
species.tree$tip.label[species.tree$tip.label == "Cortaritermes_intermedius"] <- "Cortaritermes_silvestri"
species.tree$tip.label[species.tree$tip.label == "Coptotermes_formosanus"] <- "Coptotermes_acinaciformis"
species.tree$tip.label[species.tree$tip.label == "Anacanthotermes_pacificus"] <- "Anacanthotermes_ochraceus"
```

## transferring current tree to treegraph2 for manual additions

write tree out and put in treegraph2 for changing individual species
positions using cited literature. treegraph2 can be downloaded on this
website with the following instructions:
<http://treegraph.bioinfweb.info/Download/HowToInstall>

``` r
write.tree(species.tree, file = "~/OneDrive - University College London/PhD/Traits/Tree_toedit2.nexus")
```

## bringing manually edited tree back in for final editing

Once this is brought back in we can add the species which already have a
a representative of their genus present in the phylogeny

``` r
species.tree <- read.nexus(file = "Tree_edited.nexus")


## function to add new species

## Function for adding a cherry to a tree where a single tip was before, this allows us to add a none zero value for the tips of new node
add.cherry <- function(tree, tip, new.tips) {
  
  ## Find the edge leading to the tip
  tip_id <- match(tip, tree$tip.label)
  
  ## Create the new cherry
  tree_to_add <- ape::stree(length(c(tip, new.tips)))
  
  ## Naming the tips
  tree_to_add$tip.label <- c(tip, new.tips)
  
  ## Add 0.001 branch length
  tree_to_add$edge.length <- rep(0.001, Nedge(tree_to_add))
  
  ## Binding both trees
  return(bind.tree(tree, tree_to_add, where = tip_id))
}

###this is for adding many species firstly for small phylogeny
## List of four tips to modify
tips_to_modify <- list("Reticulitermes_lucifugus", "Cryptotermes_cynocephalus",
                       "Ancistrotermes_crucifer",   "Bifiditermes_durbanensis", "Glyptotermes_brevicornis",
                       "Neotermes_sanctaecrucis",   "Pericapritermes_nigerianus",       "Trinervitermes_bettonianus",
                       "Amitermes_meridionalis",    "Microtermes_magnocellus",  "Odontotermes_transvaalensis",  
                       "Macrotermes_michaelseni",   "Nasutitermes_arborum", "Cubitermes_tenuiceps",
                       "Microcerotermes_parvus",    "Coptotermes_acinaciformis",    
                       "Bulbitermes_singaporiensis", "Neocapritermes_parvus")

## List of tips to add to these four tips
tips_to_add <- list("Reticulitermes_flavipes",  c("Cryptotermes_domesticus",    "Cryptotermes_brevis"), 
                    c("Ancistrotermes_cavithorax",  "Ancistrotermes_latinotus"),    c("Bifiditermes_mutabae",
                                                                                  "Bifiditermes_sibayensis"),   c("Glyptotermes_brevicaudatus", "Glyptotermes_ueleensis"),  
                    c("Neotermes_aburiensis",   "Neotermes_camerunensis"),  c("Pericapritermes_dolichocephalus",
                                                                           "Pericapritermes_nitobei"),  c("Trinervitermes_trinervius",  "Trinervitermes_dispar"),   c("Amitermes_laurensis",    
                                                                                                                                                                  "Amitermes_beaumonti",    "Amitermes_evuncifer",  "Amitermes_dentatus"),  c("Microtermes_traghardi",  
                                                                                                                                                                                                                                         "Microtermes_subhyalinus", "Microtermes_alluaudanus"), c("Odontotermes_sarawakensis",  "Odontotermes_badius",
                                                                                                                                                                                                                                                                                                  "Odontotermes_zambesiensis"), c("Macrotermes_gilvus", "Macrotermes_carbonarius",  "Macrotermes_natalensis",
                                                                                                                                                                                                                                                                                                                                  "Macrotermes_bellicosus"),    c("Nasutitermes_elegantulus",   "Nasutitermes_matangensis",
                                                                                                                                                                                                                                                                                                                                                                "Nasutitermes_corniger", "Nasutitermes_novarumhebridarum"), c("Cubitermes_glebae",  "Cubitermes_oculatus",  
                                                                                                                                                                                                                                                                                                                                                                                                                              "Cubitermes_serverus",    "Cubitermes_inclitus"), c("Microcerotermes_fuscotibialis",  "Microcerotermes_edentatus",    "Microcerotermes_Serrula",  "Microcerotermes_strunckii"),   
                    c("Coptotermes_sjostedti",  "Coptotermes_testaceus",    "Coptotermes_truncatus",    "Coptotermes_curvignathus"),    "Bulbitermes_perpusillus", c("Neocapritermes_opacus", "Neocapritermes_utiariti"))
  
  
  ####large phylogeny need to make lists of the larger species/add more to the above when possible
  
## Adding the tips all together
new_tree <- species.tree
for(one_tip in seq_along(tips_to_modify)) {
  new_tree <- add.cherry(new_tree, tip = tips_to_modify[[one_tip]], new.tips = tips_to_add[[one_tip]])
}

species.tree<-new_tree
```

## Resolving polytomies and creating a dichotomous ultrametric tree

``` r
resolve_species.tree<-  fix.poly(species.tree,type=c("resolve"))

tree_ultra=chronos(species.tree, lambda=0)  
is.ultrametric(tree_ultra) #now it's ultrametric.
tree.unmatched <- multi2di(tree_ultra, random=TRUE)  #This makes the tree dichotomous 
species.tree<-tree.unmatched

cdata$tiplabel<-cdata$Species

# Check whether the names match in the data and the tree
check <- name.check(phy =species.tree, data = cdata, 
                    data.names = cdata$tiplabel)
# Look at check
check

# save tree
write.tree(species.tree, file = "Species_tree.nexus")
```
