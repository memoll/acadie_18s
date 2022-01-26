###############################################################
# Cleaning and denoising 18S data                             #
# Data: Miseq-18S - all Runs - Subset L'Acadie (ACA)          #
# Mona Parizadeh - 2020-2021                                  #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(vegan); packageVersion("vegan") #‘2.5.7’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(tidyverse); packageVersion("tidyverse") #‘1.3.0’

# Import data #### 
setwd("../mp/aca_18s/files/")
ps = readRDS("ps.rds") 

#Explore controls ####
#%Positive controls ####
ct.posID = sample_names(ps)[grep("CTL..00",sample_names(ps))] 
ps.pos.ctl = subset_samples(ps, sample_data(ps)$sampleid %in% ct.posID)
ps.pos.ctl = prune_taxa(taxa_sums(ps.pos.ctl)>0, ps.pos.ctl)
#plot 
ps.pos.ctl %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #remove NAs
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Genus") +
  labs(title = "Nematode genera present in the positive controls") +
  xlab("Positive Control IDs") + ylab("Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "left")
#remove positive controls
ps.neg = subset_samples(ps, !sample_data(ps)$sampleid %in% ct.posID)
ps.neg = prune_taxa(taxa_sums(ps.neg)>0, ps.neg); ps.neg

#%Negative controls ####
ps.neg.ctl = subset_samples(ps.neg, sample_data(ps.neg)$sample_or_control == "control")
ps.neg.ctl = prune_taxa(taxa_sums(ps.neg.ctl)>0, ps.neg.ctl); ps.neg.ctl
#plot
ps.neg.ctl %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #keep NAs, glom based on genus
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Genus") +
  labs(title = "Bacteria families present in the negative controls \nAbsolute abundance") +
  xlab("Negative Control IDs") + ylab("Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1))

#Explore data ####
subset_samples(ps.neg, sample_data(ps.neg)$neonic == "Y") #neonic-treated
subset_samples(ps.neg, sample_data(ps.neg)$neonic == "N") #control (non-treated)

#number of seq per sample
summary(sample_sums(ps.neg))  #or: summary(apply(comm,1,sum))
sd(sample_sums(ps.neg), na.rm=TRUE)/sqrt(length(sample_sums(ps.neg)[!is.na(sample_sums(ps.neg))])) #SE
head(sort(sample_sums(ps.neg),TRUE))
hist(sample_sums(ps.neg))
#distribution of ASVs
hist(log10(taxa_sums(ps.neg))) 
#ASV richness 
summary(estimate_richness(ps.neg, measures = "Observed")) 
#ASV (OTU) table 
taxa_names(ps.neg) = paste0("ASV", seq(ntaxa(ps.neg))) #replace sequence w/ ASV
otu_mat = function(ps) as(otu_table(ps), "matrix")
otu_mat(ps.neg)[1:5,1:5] 
#taxonomic table 
tax_mat = function(ps) as(tax_table(ps), "matrix")
tax_mat(ps.neg)[1:5,]

#Correct names ####
#1. remove D_...
#remove D_0__
tax_table(ps.neg)[,1] = stringr::str_replace(tax_table(ps.neg)[,1], 'D_0__', '')
#remove D_1__
tax_table(ps.neg)[,2] = stringr::str_replace(tax_table(ps.neg)[,2], 'D_1__', '')
#remove D_2__
tax_table(ps.neg)[,3] = stringr::str_replace(tax_table(ps.neg)[,3], 'D_2__', '')
#remove D_3__
tax_table(ps.neg)[,4] = stringr::str_replace(tax_table(ps.neg)[,4], 'D_3__', '')
#remove D_4__
tax_table(ps.neg)[,5] = stringr::str_replace(tax_table(ps.neg)[,5], 'D_4__', '')
#remove D_5__
tax_table(ps.neg)[,6] = stringr::str_replace(tax_table(ps.neg)[,6], 'D_5__', '')
#remove D_6__
tax_table(ps.neg)[,7] = stringr::str_replace(tax_table(ps.neg)[,7], 'D_6__', '')
#2. capitalize the first letter
tax_table(ps.neg)[,1:7] = stringr::str_to_title(tax_table(ps.neg))
#3. replace "Na" taxa w/ "Unclassified"
tax_table(ps.neg)[,2][tax_table(ps.neg)[,2] == "Na"] = "Unclassified" #phylum
tax_table(ps.neg)[,5][tax_table(ps.neg)[,5] == "Na"] = "NA" #family
tax_table(ps.neg)[,6][tax_table(ps.neg)[,6] == "Na"] = "NA" #genus
#4. replace NA taxa w/ "Unclassified"
tax_table(ps.neg)[,2][is.na(tax_table(ps.neg)[,2])] = "Unclassified"
tax_table(ps.neg)[,5][is.na(tax_table(ps.neg)[,5])] = "NA"
tax_table(ps.neg)[,6][is.na(tax_table(ps.neg)[,6])] = "NA"

#Denoising data ####
#1. Keep samples with at least 1,000 reads ####
ps.neg.ctrl_1Krds = prune_samples(sample_sums(ps.neg)>=1000, ps.neg.ctrl)
ps.neg.ctrl_1Krds = prune_taxa(taxa_sums(ps.neg.ctrl_1Krds)>0, ps.neg.ctrl_1Krds); ps.neg.ctrl_1Krds
setdiff(as.vector(sample_data(ps.neg.ctrl)$sampleid),as.vector(sample_data(ps.neg.ctrl_1Krds)$sampleid)) 
which(sample_data(ps.neg.ctrl)$sample_or_control == "control")
any(sample_data(ps.neg.ctrl_1Krds)$sample_or_control == "control")

#2. Filter ASVs w/ less than 10 reads ####
ps.neg.ctrl_seq10 = prune_taxa(taxa_sums(ps.neg.ctrl_1Krds) > 10, ps.neg.ctrl_1Krds)
ps.neg.ctrl_seq10 = prune_samples(sample_sums(ps.neg.ctrl_seq10)>0, ps.neg.ctrl_seq10); ps.neg.ctrl_seq10
100 - (ntaxa(ps.neg.ctrl_seq10)/ntaxa(ps.neg.ctrl))*100 

#3.Remove problematic annotations ####
ps.neg.ctrl_seq10 %>% 
  transform_sample_counts(function(otu) otu/sum(otu)) %>%
  psmelt() %>%
  ggplot(aes(x=Family, y=Abundance/sum(Abundance), fill= Order, color=Order)) + 
  theme_bw() +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        legend.text=element_text(size=8)) +
  guides(fill=guide_legend(ncol=15))

#3.1.Remove ASVs assigned to genera that have mistakenly been assigned to the Nematoda phylum ####
#bilateria is not really a genus belonging to the nematoda phyla, it is a clade that includes nematoda and some other animals
rmGen1 = tax_table(ps.neg.ctrl_seq10)[which(tax_table(ps.neg.ctrl_seq10)[,6] =="Bilateria"),]; length(rownames(rmGen1))
ps.neg.ctrl_rmGen1 = subset_taxa(ps.neg.ctrl_seq10, !rownames(tax_table(ps.neg.ctrl_seq10)) %in% rownames(rmGen1))
#3.2. Remove ASV assigned to human parasitic nematodes ####
rmGen2 = tax_table(ps.neg.ctrl_rmGen1)[which(tax_table(ps.neg.ctrl_rmGen1)[,6] =="Eucoleus"),]; length(rownames(rmGen2))
ps.neg.ctrl_rmGen2 = subset_taxa(ps.neg.ctrl_rmGen1, !rownames(tax_table(ps.neg.ctrl_rmGen1)) %in% rownames(rmGen2)); ps.neg.ctrl_rmGen2

#4. Remove outliers ####
#4.1 NMDS ####
ps.ra = transform_sample_counts(ps.neg.ctrl_rmGen2, function(otu) otu/sum(otu))
nmds = ordinate(ps.neg.ctrl_rmGen2, method = "NMDS", k = 2, try = 100, distance = "bray") 
plot_ordination(ps.ra, nmds, color = "neonic") + 
  theme_bw() +
  geom_point() + ggtitle("nMDS") + geom_point(size = 5) +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3, nudge_y = -0.1)  #nudge_y to seperate id from point
#4.2 Alpha diversity ####
shn = estimate_richness(ps.neg.ctrl_rmGen2, split=TRUE, measures="Shannon") 
plot_richness(ps.neg.ctrl_rmGen2, "month","neonic", measures = "Shannon") +
  geom_text(aes(label = sampleid))
##remove outliers
out = c("SDA7661","SDA8681","SDA9842")
ps2 = prune_samples(!sample_data(ps.neg.ctrl_rmGen2)$sampleid %in% out, ps.neg.ctrl_rmGen2)
ps2 = prune_taxa(taxa_sums(ps2)>0, ps2); ps2

#Ordination - NMDS ####
ps2.ra = transform_sample_counts(ps2, function(otu) otu/sum(otu))
nmds2 = ordinate(ps2, method = "NMDS", k = 2, try = 100, distance = "bray") 
plot_ordination(ps2.ra, nmds2, color = "year", shape = "neonic") + 
  theme_bw() +
  geom_point() + ggtitle("nMDS") + geom_point(size = 4) +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3, nudge_y = -0.05) 
#Richeness - alpha diversity ####
plot_richness(ps2, "year","neonic", measures = "Shannon") +
  geom_text(aes(label = sampleid))

#Rarefaction ####
#rarefaction curve
source("../mp/scripts/my_functions.R")
rare_curve = calculate_rarefaction_curves(ps2,c('Observed', 'Shannon','Simpson'), c(100,200,500,1000,2000))
summary(rare_curve)
rare_curve_summary = ddply(rare_curve, 
                           c('Depth', 'Sample', 'Measure'), summarise, 
                           Alpha_diversity_mean = mean(Alpha_diversity))
obs = rare_curve_summary[rare_curve_summary$Measure=="Observed",]
#FIG.1S-A####
# Plot
p.nema.rare = ggplot(data = obs,mapping = aes(x = Depth,y = Alpha_diversity_mean,colour=Sample,group=Sample)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,2000,100))+ geom_line() + theme_bw() +theme(legend.position = "none") + 
  labs(x = "\nNumber of sequences", y = "Observed number of ASVs\n") + geom_point(size = 1) +
  theme(axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 20, face = "bold"))
p.nema.rare 
library(cowplot); packageVersion("cowplot") #"‘1.0.0’"
save_plot("/data/users/mona/miseq_18S/Mona_18S_all/article2/graphs/ParizadehM_Fig.S1A_rare.pdf",p.nema.rare, ncol = 2, nrow = 2)

ps.rare = rarefy_even_depth(ps2, sample.size = 1000, rngseed = 5006, trimOTUs = TRUE, replace = TRUE)
ps.rare
rarecurve(t(otu_table(ps.rare)), step=100,label=FALSE, col = "darkred",
          main = "Rarefy to 1,000 reads per sample")

nsamples(ps2) - nsamples(ps.rare)
ntaxa(ps2) - ntaxa(ps.rare)
#lost samples and ASVs after data denoising and rarefaction at 4,000
100-(nsamples(ps.rare)/nsamples(ps.neg.ctrl))*100 
100-(ntaxa(ps.rare)/ntaxa(ps.neg.ctrl))*100 
100-(sum(taxa_sums(ps.rare))/sum(taxa_sums(ps.neg.ctrl)))*100 #of sequences (abundance of ASVs)

#Correct phyla names ####
#Sra to SRA
tax_table(ps.rare)[,2][tax_table(ps.rare)[,2] == "Sar"] = "SAR"

#Explore denoised data ####
subset_samples(ps.rare, sample_data(ps.rare)$neonic == "Y") #neonic-treated
subset_samples(ps.rare, sample_data(ps.rare)$neonic == "N") #control (non-treated)
#number of ASVs per sample
asv.rare.rich = estimate_richness(ps.rare, measures = "Observed") #ASV per sample (richness)
summary(asv.rare.rich)
sd(asv.rare.rich$Observed, na.rm=TRUE) /  
  sqrt(length(asv.rare.rich$Observed[!is.na(asv.rare.rich$Observed)])) #SE 

#save
saveRDS(ps2, "18S_nema_aca_soil.rds")
saveRDS(ps.rare, "18S_nema_aca_soil1000.rds")
