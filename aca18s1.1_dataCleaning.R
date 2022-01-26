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
plot_bar(ps.neg, fill="Kingdom")
#kingdom pie chart 
ps.neg %>%
  tax_glom(taxrank = "Kingdom", NArm=FALSE) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt() %>%
  group_by(Kingdom) %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  ggplot(aes(x="Kingdom", y=Abundance/sum(Abundance), fill = reorder(Kingdom,Abundance))) + 
  theme_void() +
  geom_bar(width=1, size = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text_repel(aes(x = 1.2,label = paste0(round(Abundance/sum(Abundance)*100, 1), "%")), 
                  position = position_stack(vjust = 0.8),size = 4)

#1. Remove undefined kingdoms ####
ps.neg.noNa = subset_taxa(ps.neg, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned"))
ps.neg.noNa = prune_taxa(taxa_sums(ps.neg.noNa)>0, ps.neg.noNa)
ps.neg.noNa
100 - (ntaxa(ps.neg.noNa)/ntaxa(ps.neg))*100 #of ASVs
100 - (sum(taxa_sums(ps.neg.noNa))/sum(taxa_sums(ps.neg)))*100 #of sequences (abundance of ASVs)

#2. Remove algae and fungi kingdoms 
ps.neg_euk = subset_taxa(ps.neg.noNa, Kingdom == "Metazoa" | Kingdom =="Eukaryota")
ps.neg_euk = prune_taxa(taxa_sums(ps.neg_euk)>0, ps.neg_euk)
ps.neg_euk
(ntaxa(ps.neg_euk)/ntaxa(ps.neg))*100 #of ASVs
100 - (sum(taxa_sums(ps.neg_euk))/sum(taxa_sums(ps.neg)))*100 #of sequences (abundance of ASVs)

#3. Keep samples with at least 1,000 reads ####
ps.neg_1Krds = prune_samples(sample_sums(ps.neg_euk)>=1000, ps.neg_euk)
ps.neg_1Krds = prune_taxa(taxa_sums(ps.neg_1Krds)>0, ps.neg_1Krds)
ps.neg_1Krds
setdiff(as.vector(sample_data(ps.neg.noNa)$sampleid),as.vector(sample_data(ps.neg_1Krds)$sampleid)) 

#4. Filter ASVs w/ less than 10 reads ####
ps.neg_seq10 = prune_taxa(taxa_sums(ps.neg_1Krds) > 10, ps.neg_1Krds)
ps.neg_seq10 = prune_samples(sample_sums(ps.neg_seq10)>0, ps.neg_seq10) 
ps.neg_seq10
100 - (ntaxa(ps.neg_seq10)/ntaxa(ps.neg))*100 

#5.Remove problematic annotations ####
ps.neg_seq10 %>% 
  transform_sample_counts(function(otu) otu/sum(otu)) %>%
  psmelt() %>%
  ggplot(aes(x=Genus, y=Abundance/sum(Abundance), fill= Family, color=Family)) + 
  theme_bw() +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        legend.text=element_text(size=8)) +
  guides(fill=guide_legend(ncol=15))

#5.1.Remove ASVs assigned to a genus that belongs to two different families (based on the barplot stacks, one stack w/ two colors) ####
#5.1.1.bilateria is not really a genus belonging to the nematoda phyla, it is a clade that includes nematoda and some other animals
rmGen1 = tax_table(ps.neg_seq10)[which(tax_table(ps.neg_seq10)[,6] =="Bilateria"),]; dim(rmGen1)
ps.neg_rmGen1 = subset_taxa(ps.neg_seq10, !rownames(tax_table(ps.neg_seq10)) %in% rownames(rmGen1))
#prune_samples(sample_sums(ps.neg_rmGen1)>0, ps.neg_rmGen1)
ps.neg_rmGen1
#5.1.2.heteromita belongs to the heteromitidae family, but in ASV233 it has been assigned to the viridiraptoridae family
tax_table(ps.neg_rmGen1)[which(tax_table(ps.neg_rmGen1)[,6] =="Heteromita"),]
rmGen2 = tax_table(ps.neg_rmGen1)[which(tax_table(ps.neg_rmGen1)[,5] =="Viridiraptoridae" & 
                                             tax_table(ps.neg_rmGen1)[,6] =="Heteromita"),];dim(rmGen2)
ps.neg_rmGen2 = subset_taxa(ps.neg_rmGen1, !rownames(tax_table(ps.neg_rmGen1)) %in% rownames(rmGen2))
ps.neg_rmGen2
#5.2. Remove ASV assigned to human parasitic nematodes
rmGen3 = tax_table(ps.neg_rmGen2)[which(tax_table(ps.neg_rmGen2)[,6] =="Eucoleus"),]; dim(rmGen3)
ps.neg_rmGen3 = subset_taxa(ps.neg_rmGen2, !rownames(tax_table(ps.neg_rmGen2)) %in% rownames(rmGen3))
#prune_samples(sample_sums(ps.neg_rmGen3)>0, ps.neg_rmGen3)
ps.neg_rmGen3

#6. Remove outliers ####
#6.1 NMDS ####
ps.ra = transform_sample_counts(ps.neg_rmGen3, function(otu) otu/sum(otu))
nmds = ordinate(ps.neg_rmGen3, method = "NMDS", k = 2, try = 100, distance = "bray") 
plot_ordination(ps.ra, nmds, color = "neonic") + 
  theme_bw() +
  geom_point() + ggtitle("nMDS") + geom_point(size = 5) +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3, nudge_y = -0.1)  #nudge_y to seperate id from point
#6.2 Alpha diversity ####
shn = estimate_richness(ps.neg_rmGen3, split=TRUE, measures="Shannon") 
rownames(shn)[(which(shn$Shannon<1))]
plot_richness(ps.neg_rmGen3, "month","neonic", measures = "Shannon") +
  geom_text(aes(label = sampleid))
#remove outliers
out = c("SDA7661","CDA7711","SDA9842","CTL0004","CTL4001")
ps2 = prune_samples(!sample_data(ps.neg_rmGen3)$sampleid %in% out, ps.neg_rmGen3)
ps2 = prune_taxa(taxa_sums(ps2)>0, ps2)
ps2

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

#Phyla bars ####
#phyla - proportions in the whole dataset
phlm2.melt = ps2 %>% 
  tax_glom(taxrank = "Phylum", NArm=FALSE) %>%
  transform_sample_counts(function(otu) otu/sum(otu)) %>%
  psmelt()
#plot
#plot to evaluate if there are genera with the same name belonging to two different families
colourCount2 = length(unique(phlm2.melt$Phylum))
ggplot(phlm2.melt, aes(x=Phylum, y=Abundance/sum(Abundance), fill=Phylum, color=Phylum)) + 
  theme_bw() +
  geom_bar(stat = "identity") +
  #facet_wrap(~neonic) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(colourCount2)) +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(colourCount2)) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1)) 

#Rarefaction ####
#rarefaction curve
source("../mp/scripts/my_functions.R")

rare_curve = calculate_rarefaction_curves(ps2,c('Observed', 'Shannon','Simpson'), c(1000,2000,5000,10000,20000))
summary(rare_curve)
rare_curve_summary = ddply(rare_curve, 
                           c('Depth', 'Sample', 'Measure'), summarise, 
                           Alpha_diversity_mean = mean(Alpha_diversity))
obs = rare_curve_summary[rare_curve_summary$Measure=="Observed",]
#FIG.1S-A####
# Plot
p.euk.rare = ggplot(data = obs,mapping = aes(x = Depth,y = Alpha_diversity_mean,colour=Sample,group=Sample)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,20000,1000))+ geom_line() + theme_bw() +theme(legend.position = "none") + 
  labs(x = "\nNumber of sequences", y = "Observed number of ASVs\n") + geom_point(size = 1) +
  theme(axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 20, face = "bold"))
p.euk.rare 
library(cowplot); packageVersion("cowplot") #"‘1.1.1’"
save_plot("graphs_final/ParizadehM_Fig.S1A_rare.pdf",p.euk.rare, ncol = 2, nrow = 2)

ps.rare = rarefy_even_depth(ps2, sample.size = 4000, rngseed = 5006, trimOTUs = TRUE, replace = TRUE)
ps.rare
rarecurve(t(otu_table(ps.rare)), step=100,label=FALSE, col = "darkred",
          main = "Rarefy to 4,000 reads per sample")

nsamples(ps2) - nsamples(ps.rare)
ntaxa(ps2) - ntaxa(ps.rare)
#lost samples and ASVs after data denoising and rarefaction at 4,000
100-(nsamples(ps.rare)/nsamples(ps.neg))*100 
100-(ntaxa(ps.rare)/ntaxa(ps.neg))*100 
100-(sum(taxa_sums(ps.rare))/sum(taxa_sums(ps.neg)))*100 #of sequences (abundance of ASVs)

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

#***RESULTS 3.1 ####
#Subset Control ####
ps.ctl = subset_samples(ps.rare, sample_data(ps.rare)$neonic == "N")
ps.ctl = prune_taxa(taxa_sums(ps.ctl)>0, ps.ctl)
ps.ctl

#plotbar-one stack ####
p.bar_all = ps.ctl %>%
  tax_glom(taxrank = "Phylum", NArm=FALSE) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt() %>%
  group_by(Phylum) %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  arrange(desc(Abundance)) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100)  %>%
  ggplot(aes(x="Phylum", y=Abundance/sum(Abundance), fill = reorder(Phylum,Abundance))) + 
  geom_bar(width=0.5, size = 1, stat = "identity") +
  #geom_text_repel(aes(x = 1,label = paste0(round(Abundance/sum(Abundance)*100, 1), "%")), 
  #          position = position_stack(vjust = 0.5),size = 4) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(name = "Phylum",
                    values = c("darkblue","tomato3","darkorange","bisque1",
                               "darkred","mediumorchid2","azure3","chartreuse2","tan2",
                               "seagreen","darkgoldenrod1","darkolivegreen4","lightsalmon","mediumvioletred",
                               "cyan3","indianred1","cornflowerblue","aquamarine3")) +
  # geom_text_repel(aes(x = 1,label = paste0(round(Abundance/sum(Abundance)*100, 1), "%")), 
  #                 position = position_stack(vjust = 0.8),size = 4) +
  xlab("Phyla") +
  ylab("Relative Abundance") +
  theme_bw() +
  theme(aspect.ratio = 1, #width (thinner)
        legend.title = element_text(size=14, face="bold"),
        legend.text=element_text(size=12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold")) +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 20, face = "bold")) 
p.bar_all

#plotbar - all samples ####
p.bar_smp = ps.ctl %>%
  tax_glom(taxrank = "Phylum", NArm=FALSE) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt() %>%
  mutate(Sample,Abundance/sum(Abundance)) %>%
  ggplot(aes(x=Sample, y=Abundance, fill = Phylum)) + 
  theme_bw() +
  geom_bar(width=0.7, size = 1, stat = "identity") +
  scale_fill_manual(name = "Phylum",
                    values = c("seagreen","cyan3","darkred","indianred1","bisque1",
                               "aquamarine3","mediumorchid2","darkorange","darkblue",
                               "darkolivegreen4","cornflowerblue","azure3","tomato3",
                               "darkgoldenrod1","mediumvioletred","lightsalmon","tan2","chartreuse2")) +
  xlab("Samples") +
  ylab("Relative Abundance") +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 0, vjust=0.5),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold")) +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 20, face = "bold")) 
p.bar_smp

#save
saveRDS(ps2, "18S_all_aca_soil.rds")
saveRDS(ps.rare, "18S_all_aca_soil4000.rds")
save.image("a2_1_qc1_aca_18S_all.RData")


