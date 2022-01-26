###############################################################
# Cleaning and denoising 16S data                             #
# Data: Miseq-16S - all Runs - Subset L'Acadie (ACA)          #
# Mona Parizadeh - 2020-2021                                  #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(vegan); packageVersion("vegan") #‘2.5.7’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(tidyverse); packageVersion("tidyverse") #‘1.3.0’

# Import data ####  
setwd("../mp/aca_16s/files/")
ps = readRDS("ps_t.rds") 
ps

#subset data (remove leaf, leaf bag control, May, other sites data) 
ps = subset_samples(ps, sample_data(ps)$habitat != "leaf" & sample_data(ps)$month != "May")
ps = prune_samples(!sample_data(ps)$sampleid %in% "CTL0001", ps) #leaf bag control
ctrl = sample_names(ps)[grep("CTL....",sample_names(ps))] 
ps = subset_samples(ps, sample_data(ps)$site == "ACA" | sample_data(ps)$sampleid %in% ctrl)
ps = prune_taxa(taxa_sums(ps)>0, ps)
sample_data(ps)$habitat = NULL
sample_data(ps)$site = NULL
sample_data(ps)$stage = NULL

#Explore positive controls ####
#subset Positive controls 
ct.posID = sample_names(ps)[grep("CTL..00",sample_names(ps))] 
ps.pos.ctl = subset_samples(ps, sample_data(ps)$sampleid %in% ct.posID)
ps.pos.ctl = prune_taxa(taxa_sums(ps.pos.ctl)>0, ps.pos.ctl)
#plot 
ps.pos.ctl %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #remove NAs
  plot_bar(fill="Genus") +
  labs(title = "Bacteria genera present in the positive controls") +
  xlab("Positive Control IDs") + ylab("Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "left")
# Remove positive controls
ps.neg_ctrl = subset_samples(ps, !sample_data(ps)$sampleid %in% ct.posID)
ps.neg_ctrl = prune_taxa(taxa_sums(ps.neg_ctrl)>0, ps.neg_ctrl) 
ps.neg_ctrl

# Explore data ####
subset_samples(ps.neg_ctrl, sample_data(ps.neg_ctrl)$neonic == "Y") #neonic-treated
subset_samples(ps.neg_ctrl, sample_data(ps.neg_ctrl)$neonic == "N") #control (non-treated)
subset_samples(ps.neg_ctrl, sample_data(ps.neg_ctrl)$sample_or_control == "control") #neg ctrl
#number of seq per sample
summary(sample_sums(ps.neg_ctrl))  #or: summary(apply(comm,1,sum))
sd(sample_sums(ps.neg_ctrl), na.rm=TRUE)/sqrt(length(sample_sums(ps.neg_ctrl)[!is.na(sample_sums(ps.neg_ctrl))])) #SE
head(sort(sample_sums(ps.neg_ctrl),TRUE))
hist(sample_sums(ps.neg_ctrl))
#ASV richness 
summary(estimate_richness(ps.neg_ctrl, measures = "Observed"))
#distribution of ASVs
hist(log10(taxa_sums(ps.neg_ctrl))) 
# ASV (OTU) table 
taxa_names(ps.neg_ctrl) = paste0("ASV", seq(ntaxa(ps.neg_ctrl))) #replace sequence w/ ASV
otu_mat = function(ps) as(otu_table(ps), "matrix")
otu_mat(ps.neg_ctrl)[1:5,1:5] 
# taxonomic table 
tax_mat = function(ps) as(tax_table(ps), "matrix")
tax_mat(ps.neg_ctrl)[1052:1062,]

#Explore contaminants ####
#Find Archaea 
Archaea = subset_taxa(ps.neg_ctrl, Kingdom == "Archaea")
(ntaxa(Archaea)/ntaxa(ps.neg_ctrl))*100 
# Find Cyanobacteria 
Cyanobacteria = subset_taxa(ps.neg_ctrl, Phylum=="Cyanobacteria")
(ntaxa(Cyanobacteria)/ntaxa(ps.neg_ctrl))*100 
# Find Chloroplast 
Chloroplast = subset_taxa(ps.neg_ctrl, Order=="Chloroplast")  ####keep it for later
(ntaxa(Chloroplast)/ntaxa(ps.neg_ctrl))*100 
# Find Mitochondria 
Mitochondria = subset_taxa(ps.neg_ctrl, Family=="Mitochondria")
(ntaxa(Mitochondria)/ntaxa(ps.neg_ctrl))*100
# Note: we still keep chloroplast and mitochondria for decontam purposes 
((ntaxa(Chloroplast)+ntaxa(Mitochondria))/ntaxa(ps.neg_ctrl))*100 #0.07%

#Denoising & Decontaminating data ####
ps.neg_ctrl.phlm = ps.neg_ctrl %>%
  tax_glom(taxrank = "Phylum", NArm=FALSE) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt() %>%
  group_by(Phylum) %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  arrange(desc(Abundance)) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100) 
ps.neg_ctrl.phlm$Phylum #find NA
ps.neg_ctrl.phlm[13,]
#1. Remove undefined phyla ####
ps.neg_ctrl.noNa = subset_taxa(ps.neg_ctrl, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps.neg_ctrl.noNa = prune_samples(sample_sums(ps.neg_ctrl.noNa)>0, ps.neg_ctrl.noNa)
ps.neg_ctrl.noNa
100-(ntaxa(ps.neg_ctrl.noNa)/ntaxa(ps.neg_ctrl))*100 #0.63% 
100-(sum(taxa_sums(ps.neg_ctrl.noNa))/sum(taxa_sums(ps.neg_ctrl)))*100 #of sequences (abundance of ASVs)

#2. Keep samples with at least 1000 reads ####
ps.neg_ctrl.1krds=prune_samples(sample_sums(ps.neg_ctrl.noNa)>=1000, ps.neg_ctrl.noNa)
ps.neg_ctrl.1krds=prune_taxa(taxa_sums(ps.neg_ctrl.1krds)>0, ps.neg_ctrl.1krds)
ps.neg_ctrl.1krds
#lost samples
setdiff(as.vector(sample_data(ps.neg_ctrl.noNa)$sampleid),as.vector(sample_data(ps.neg_ctrl.1krds)$sampleid)) #sample
ps.lost_1kreads = prune_samples(sample_sums(ps.neg_ctrl.noNa)<1000, ps.neg_ctrl.noNa)
ps.lost_1kreads
sample_names(ps.lost_1kreads)[grep("CTL..0.",sample_names(ps.lost_1kreads))] #lost
sample_names(ps.neg_ctrl.1krds)[grep("CTL..0.",sample_names(ps.neg_ctrl.1krds))] #left

#3. Remove contaminants ####
library(decontam); packageVersion("decontam") #‘1.1.2’
# 3.1. Inspect library sizes ####
# Assign the full sample_data
sample_df = function(ps) as(sample_data(ps), "data.frame")
sam.df = sample_df(ps.neg_ctrl.1krds)  # Put sample_data into a ggplot-friendly data.frame
sam.df$LibrarySize = sample_sums(ps.neg_ctrl.1krds)
sample_data(ps.neg_ctrl.1krds) = sample_data(sam.df)
# The shorthand way, assign one column
sample_data(ps.neg_ctrl.1krds)$LibrarySize = sample_sums(ps.neg_ctrl.1krds)
df = sam.df[order(sam.df$LibrarySize),]
df$Index = seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point()
# 3.2. Prevalence - Identify Contaminants ####
sample_data(ps.neg_ctrl.1krds)$is.neg = sample_data(ps.neg_ctrl.1krds)$sample_or_control == "control"
contamdf.prev = isContaminant(ps.neg_ctrl.1krds, method="prevalence", neg="is.neg")
# the default threshold for a contaminant is that it reaches a probability of 0.1
table(contamdf.prev$contaminant) #ASVs being recognized as contaminant (TRUE) or not (FALSE)
which(contamdf.prev$contaminant) 
# getting vector holding the identified contaminant IDs
contam_asvs = row.names(contamdf.prev[contamdf.prev$contaminant == TRUE, ])
length(contam_asvs)
contam_taxa = as.data.frame(tax_mat(ps.neg_ctrl.1krds))
head(contam_taxa)

#In the prevalence test there is a special value worth knowing,  threshold=0.5
#This will identify as contaminants all sequences that are more prevalent in negative controls than in samples.
contamdf.prev05 = isContaminant(ps.neg_ctrl.1krds, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
which(contamdf.prev05$contaminant)
# getting vector holding the identified contaminant IDs
contam05_asvs = row.names(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])
length(contam05_asvs)
contam05_taxa = as.data.frame(tax_mat(ps.neg_ctrl.noNa))
head(contam05_taxa)

# Visualize abundance of potential contaminant ASVs - prevalence cutoff 0.1
ps.pa = transform_sample_counts(ps.neg_ctrl.1krds, function(abund) 1*(abund>0))
ps.pa_ctl = prune_samples(sample_data(ps.pa)$sample_or_control == "control", ps.pa) # negative samples
sample_data(ps.pa_ctl)$sampleid
ps.pa_smp = prune_samples(sample_data(ps.pa)$sample_or_control == "sample", ps.pa)  # the real samples
nsamples(ps.pa_smp)
nsamples(ps.pa_ctl)
# Make data.frame of prevalence in positive and negative samples
df.pa = data.frame(pa.smp=taxa_sums(ps.pa_smp), pa.ctl=taxa_sums(ps.pa_ctl),
                   contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.ctl, y=pa.smp, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  ggtitle("Abundance of potential contaminant ASVs - prevalence cutoff 0.1")

# Visualize abundance of potential contaminant ASVs - prevalence cutoff 0.5
#threshold (cutoff): The probability threshold below which the null-hypothesis (not a contaminant) should be rejected in favor of the alternate hypothesis (contaminant)
# Make data.frame of prevalence in positive and negative samples
df.pa05 = data.frame(pa.smp=taxa_sums(ps.pa_smp), pa.ctl=taxa_sums(ps.pa_ctl),
                     contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa05, aes(x=pa.ctl, y=pa.smp, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + 
  ggtitle("Abundance of potential contaminant ASVs - prevalence cutoff 0.5") 
#Samples (red points) seem to split pretty cleanly into two branches that show up mostly in positive samples and another that shows up mostly in negative controls, and the contaminant assignment (at default probability threshold) has done a good job of identifying those mostly in negative controls.

# taxonomic identity and abundance of the potential contaminant ASVs - prevalence cutoff 0.1
subset_contam = subset_taxa(ps.neg_ctrl.1krds, contamdf.prev$contaminant==TRUE)
plot_bar(subset_contam, fill="Genus") + ggtitle("Taxonomic identity and abundance of the potential contaminant ASVs \nat Genus level - prevalence cutoff 0.1")
# Abundance of the contaminant ASVs by taxonomy
taxo.contam = as.data.frame(as(tax_table(subset_contam),"matrix"))
taxo.contam$abund = apply(otu_table(subset_contam),2,sum)
rownames(taxo.contam) = NULL
head(taxo.contam)

# taxonomic identity and abundance of the potential contaminant ASVs - prevalence cutoff 0.5
subset_contam05 = subset_taxa(ps.neg_ctrl.1krds, contamdf.prev05$contaminant==TRUE)
plot_bar(subset_contam05, fill="Genus") + ggtitle("Taxonomic identity and abundance of the potential contaminant ASVs \nat Genus level - prevalence cutoff 0.5")
# Abundance of the contaminant ASVs by taxonomy
taxo.contam05 = as.data.frame(as(tax_table(subset_contam05),"matrix"))
taxo.contam05$abund = apply(otu_table(subset_contam05),2,sum)
rownames(taxo.contam05) = NULL
head(taxo.contam)
# 3.2.a. Prevalence 0.1 ####
ps.notcontam = subset_taxa(ps.neg_ctrl.1krds, contamdf.prev$contaminant==FALSE)
ps.notcontam = prune_samples(sample_sums(ps.notcontam)>0, ps.notcontam)
ps.notcontam
# 3.2.b. Agressive Prevalence 0.5 ####
ps.notcontam05 = subset_taxa(ps.neg_ctrl.1krds, contamdf.prev05$contaminant==FALSE)
ps.notcontam05 = prune_samples(sample_sums(ps.notcontam05)>0, ps.notcontam05)
ps.notcontam05
#lost ASVs
ntaxa(ps.neg_ctrl.1krds) - ntaxa(ps.notcontam05) 
100-(ntaxa(ps.notcontam05)/ntaxa(ps.neg_ctrl.1krds))*100 #0.12% 
100 - (sum(taxa_sums(ps.notcontam05))/sum(taxa_sums(ps.neg_ctrl.1krds)))*100 #of sequences (abundance of ASVs)
# 3.2.c. Remove Chloroplast ####
ps.notcontamChlp = subset_taxa(ps.notcontam05, (Order!="Chloroplast" | is.na(Order)))
# 3.2.d. Remove Mitochondria ####
ps.notcontamChlpMitc = subset_taxa(ps.notcontamChlp, (Family!="Mitochondria") | is.na(Family))
ps.notcontamChlpMitc = prune_samples(sample_sums(ps.notcontamChlpMitc) > 0, ps.notcontamChlpMitc)
ps.notcontamChlpMitc

#4. Filter ASVs w/ less than 10 reads ####
ps.neg.ctrl_seq10 = prune_taxa(taxa_sums(ps.notcontamChlpMitc) > 10, ps.notcontamChlpMitc)
ps.neg.ctrl_seq10 = prune_samples(sample_sums(ps.neg.ctrl_seq10)>0, ps.neg.ctrl_seq10) 
ps.neg.ctrl_seq10
100 - (ntaxa(ps.neg.ctrl_seq10)/ntaxa(ps.neg_ctrl))*100 
100 - (sum(taxa_sums(ps.neg.ctrl_seq10))/sum(taxa_sums(ps.neg_ctrl)))*100 #of sequences (abundance of ASVs)

#5. Remove outliers ####
#5.1 NMDS ####
ps.ra = transform_sample_counts(ps.neg.ctrl_seq10, function(otu) otu/sum(otu))
nmds = ordinate(ps.neg.ctrl_seq10, method = "NMDS", k = 2, try = 100, distance = "bray") 
plot_ordination(ps.ra, nmds, color = "neonic") + 
  theme_bw() +
  geom_point() + ggtitle("nMDS") + geom_point(size = 5) +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3, nudge_y = -0.01)  #nudge_y to seperate id from point
#5.2 Alpha diversity ####
plot_richness(ps.neg.ctrl_seq10, "month","neonic", measures = "Shannon") +
  geom_text(aes(label = sampleid))
#remove outliers
out = c("CTL0002", "SDA7831", "SDA7852", "SDA8841")
ps2 = prune_samples(!sample_data(ps.neg.ctrl_seq10)$sampleid %in% out, ps.neg.ctrl_seq10)
ps2 = prune_taxa(taxa_sums(ps2)>0, ps2)
ps2
#nsamples(ps.neg.ctrl_seq10) - nsamples(ps2)
ntaxa(ps.neg.ctrl_seq10) - ntaxa(ps2)

#Rarefaction ####
#rarefaction curve
source("/data/users/mona/useful_scripts/my_functions.R")
rare_curve = calculate_rarefaction_curves(ps2,c('Observed', 'Shannon','Simpson'), c(1000,2000,5000,10000,20000))
summary(rare_curve)
rare_curve_summary = ddply(rare_curve, 
                           c('Depth', 'Sample', 'Measure'), summarise, 
                           Alpha_diversity_mean = mean(Alpha_diversity))
obs = rare_curve_summary[rare_curve_summary$Measure=="Observed",]
#FIG.1S-B####
# Plot
p.bac.rare = ggplot(data = obs,mapping = aes(x = Depth,y = Alpha_diversity_mean,colour=Sample,group=Sample)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,20000,1000))+ geom_line() + theme_bw() +theme(legend.position = "none") + 
  labs(x = "\nNumber of sequences", y = "Observed number of ASVs\n") + geom_point(size = 1) +
  theme(axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 18, face = "bold")) 
p.bac.rare

ps.rare10 = rarefy_even_depth(ps2, sample.size = 10000, rngseed = 8306, trimOTUs = TRUE, replace = TRUE)
ps.rare10
rarecurve(otu_table(ps.rare10), step=100,label=FALSE, col = "darkred",
          main = "Rarefy to 10,000 reads per sample")

nsamples(ps2) - nsamples(ps.rare10)
ntaxa(ps2) - ntaxa(ps.rare10)
#lost samples and ASVs after data denoising and rarefaction at 10,000
100-(nsamples(ps.rare10)/nsamples(ps.neg_ctrl))*100
100-(ntaxa(ps.rare10)/ntaxa(ps.neg_ctrl))*100 

#Explore denoised data ####
subset_samples(ps.rare10, sample_data(ps.rare10)$neonic == "Y") #neonic-treated
subset_samples(ps.rare10, sample_data(ps.rare10)$neonic == "N") #control (non-treated)
#number of ASVs per sample
asv.rare.rich = estimate_richness(ps.rare10, measures = "Observed") #ASV per sample (richness)
summary(asv.rare.rich)
sd(asv.rare.rich$Observed, na.rm=TRUE) /  
  sqrt(length(asv.rare.rich$Observed[!is.na(asv.rare.rich$Observed)])) #SE 

saveRDS(ps2, "16S_aca_soil.rds")
saveRDS(ps.rare10, "16S_aca_soil10000.rds")

