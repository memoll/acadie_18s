###############################################################
# Data preparation for network analysis - SPARCC              #
# Data: Miseq-18S & 16S - L'Acadie (ACA)                      #
# Mona Parizadeh & Mathieu landry - 2020-2021                 #
###############################################################

#Data preparation steps ####
#1. agglomorate rarefied data at family-level (bacteria, nematodes)
#2. merge rarefied bacteria and nematodes 
#3. keep only the mutual samples

#Load libraries ####
library(vegan); packageVersion("vegan") #‘2.5.7’
library(permute); packageVersion("permute") #‘0.9.5’
require(lattice)
require(nlme)
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(picante); packageVersion("picante") #‘1.8.2’
library(SpiecEasi); packageVersion("SpiecEasi") #‘1.1.1’

# Import rarefied data #### 
setwd("../mp/aca_18s/files/")
#Nematodes ####
ps.nema = readRDS("18S_nema_aca_soil1000.rds")
ps.nema
ps.nema.fam = tax_glom(ps.nema, "Family") #agglomerate at family level
ps.nema.fam
#Bacteria ####
ps.bac = readRDS("16S_aca_soil10000.rds")
ps.bac = phyloseq(otu_table(ps.bac),tax_table(ps.bac),sample_data(ps.bac)) #remove tree
ps.bac
ps.bac.fam = tax_glom(ps.bac, "Family") #agglomerate at family level
ps.bac.fam

#Merge nematode and bacteria ####
# keep only the samples in common between both becatria and nematode datasets
ps.bac2 = prune_samples(sample_names(ps.nema.fam), ps.bac.fam)
ps.bac2 = prune_taxa(taxa_sums(ps.bac2)>0, ps.bac2); ps.bac2
taxa_names(ps.bac2) = paste0("16S_", taxa_names(ps.bac2));taxa_names(ps.bac2)
ps.nem2 = prune_samples(sample_names(ps.bac.fam), ps.nema.fam)
ps.nem2 = prune_taxa(taxa_sums(ps.nem2)>0, ps.nem2); ps.nem2
taxa_names(ps.nem2) = paste0("18S_", taxa_names(ps.nem2));taxa_names(ps.nem2)

otu.bac = t(otu_table(ps.bac2))
otu.nem = otu_table(ps.nem2)
#setdiff(colnames(otu.bac),colnames(otu.nem))
identical(colnames(otu.bac),colnames(otu.nem)) #sampleIDs are not in the same order

tax.bac = tax_table(ps.bac2)
tax.nem = tax_table(ps.nem2)
#identical(colnames(tax.bac),colnames(tax.nem))  

sample_data(ps.bac2)$LibrarySize = NULL
sample_data(ps.bac2)$is.neg = NULL
sample_data(ps.bac2)$treatment = NULL
meta.bac = sample_data(ps.bac2)
meta.nem = sample_data(ps.nem2)
identical(rownames(meta.bac),rownames(meta.nem)) #sampleIDs are not in the same order

#match rownames of bacteria matrix w/ nematode matrix
otu.bac2 = otu.bac[,match(colnames(otu.nem), colnames(otu.bac))]; dim(otu.bac2)
identical(colnames(otu.bac2),colnames(otu.nem)) #TRUE

otu = rbind(otu.bac2,otu.nem); dim(otu)
tax = rbind(tax.bac,tax.nem); dim(tax)
identical(colnames(otu),rownames(meta.nem)) #TRUE
#merge tables
ps.all = phyloseq(t(otu_table(otu, taxa_are_rows = TRUE)), sample_data(meta.nem), 
               tax_table(tax))
ps.all

#Subset treatments ####
ps.ctl = subset_samples(ps.all,sample_data(ps.all)$neonic == "N")
ps.ctl = prune_taxa(taxa_sums(ps.ctl)>0,ps.ctl);ps.ctl
ps.neo = subset_samples(ps.all,sample_data(ps.all)$neonic == "Y")
ps.neo = prune_taxa(taxa_sums(ps.neo)>0,ps.neo);ps.neo
#make picante objects ####
#control
comm.ctl = otu_table(ps.ctl)
taxo.ctl = tax_table(ps.ctl)
meta.ctl = sample_data(ps.ctl)
#neonic
comm.neo = otu_table(ps.neo)
taxo.neo = tax_table(ps.neo)
meta.neo = sample_data(ps.neo)

#SparCC####
#%CONTROL ####
#occurence filtration #### 
#transform into presence/absence then look at distribution of columns totals) 
pa.ctl = decostand(comm.ctl,method = "pa")
hist(apply(pa.ctl,2,sum))
#try occurence in 5 samples
comm.ctl5 = comm.ctl[,apply(pa.ctl,2,sum) > 5] #asv presented in more than 5 samples
quantile(apply(comm.ctl5, 1,sum))
meta.ctl5 = meta.ctl[rownames(comm.ctl5),]
taxo.ctl5 = taxo.ctl[colnames(comm.ctl5),]

#make a phyloseq object ####
OTU.ctl = otu_table(comm.ctl5, taxa_are_rows = FALSE)
TAX.ctl = tax_table(as.matrix(taxo.ctl5))
DATA.ctl = sample_data(meta.ctl5)
ps.ctl1 = phyloseq(OTU.ctl, TAX.ctl, DATA.ctl)
ps.ctl1
saveRDS(ps.ctl1,"ps.ctl5.rds")
#calculate bootstraps of SparCC correlation coefficients ####
#Bootstrapping uses the observed data to simulate resampling from the population. 
#This produces a large number of bootstrap resamples. 
#We can calculate a statistic for each bootstrap resample and 
#use the distribution of the simulated statistics to approximate characteristics of the population. 
tp0.ctl = proc.time()
set.seed(106)
sparcc.boot.ctl1000 = sparccboot(otu_table(ps.ctl1), R=1000, ncpus=10)
tp1.ctl1000 <- proc.time()
tp1.ctl1000 - tp0.ctl
saveRDS(sparcc.boot.ctl1000,"sparcc.boot.ctl1000.rds")

#%NEONIC ####
#occurence filtration #### 
#transform into presence/absence then look at distribution of columns totals) 
pa.neo = decostand(comm.neo,method = "pa")
hist(apply(pa.neo,2,sum))
#try occurence in 5 samples
comm.neo5 = comm.neo[,apply(pa.neo,2,sum) > 5]
quantile(apply(comm.neo5, 1,sum))
meta.neo5 = meta.neo[rownames(comm.neo5),]
taxo.neo5 = taxo.neo[colnames(comm.neo5),]

#make a phyloseq object ####
OTU.neo = otu_table(comm.neo5, taxa_are_rows = FALSE)
TAX.neo = tax_table(as.matrix(taxo.neo5))
DATA.neo = sample_data(meta.neo5)
ps.neo1 = phyloseq(OTU.neo, TAX.neo, DATA.neo)
ps.neo1
saveRDS(ps.neo1,"ps.neo5.rds")
#calculate bootstraps of SparCC correlation coefficients ####
tp0.neo = proc.time()
set.seed(107)
sparcc.boot.neo1000 = sparccboot(otu_table(ps.neo1), R=1000, ncpus=10)
tp1.neo1000 <- proc.time()
tp1.neo1000 - tp0.neo
saveRDS(sparcc.boot.neo1000,"sparcc.boot.neo1000.rds")


