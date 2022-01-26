############################################################################################
# Explanatory analysis of soil samples                                                     #
# Studying soil nematode community variations in response to neonicotinoid seed treatment  #
# Data: Miseq-18S - L'Acadie (ACA)                                                         #
# Mona Parizadeh - 2020-2021                                                               #
############################################################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(vegan); packageVersion("vegan") #‘2.5.7’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(dplyr); packageVersion("dplyr") #‘1.0.4’

setwd("../mp/aca_18s/files/")
ps.nema = readRDS("18S_nema_aca_soil1000.rds")
ps.nema

#ASV richness 
asv.rich = estimate_richness(ps.nema, measures = "Observed")
summary(asv.rich)
sd(asv.rich$Observed, na.rm=TRUE) /  
  sqrt(length(asv.rich$Observed[!is.na(asv.rich$Observed)])) #SE
#distribution of ASVs
hist(log10(taxa_sums(ps.nema))) 
# ASV (OTU) table 
#taxa_names(ps.nema) = paste0("ASV", seq(ntaxa(ps.nema))) #replace sequence w/ ASV
otu_mat = function(ps) as(otu_table(ps), "matrix")
otu_mat(ps.nema)[1:5,1:5] 
# taxonomic table 
tax_mat = function(ps) as(tax_table(ps), "matrix")
tax_mat(ps.nema)[1:5,]

#ASV richness ####
asv.rich.nema = estimate_richness(ps.nema, measures = "Observed")
summary(asv.rich.nema)
sd(asv.rich.nema$Observed, na.rm=TRUE) /  
  sqrt(length(asv.rich.nema$Observed[!is.na(asv.rich.nema$Observed)])) #SE
# Treatments ####
ps.nema.neo = subset_samples(ps.nema, sample_data(ps.nema)$neonic == "Y") #neonic-treated
ps.nema.neo = prune_taxa(taxa_sums(ps.nema.neo)>0,ps.nema.neo); ps.nema.neo
ps.nema.ctl = subset_samples(ps.nema, sample_data(ps.nema)$neonic == "N") #control (non-treated)
ps.nema.ctl = prune_taxa(taxa_sums(ps.nema.ctl)>0,ps.nema.ctl); ps.nema.ctl

#Look for outliers ####
#NMDS ordination ####
ps.ra = transform_sample_counts(ps.nema, function(otu) otu/sum(otu))
nmds = ordinate(ps.nema, method = "NMDS", k = 2, try = 100, distance = "bray") 
plot_ordination(ps.ra, nmds, color = "year", shape = "neonic") + 
  theme_bw() +
  geom_point() + ggtitle("nMDS") + geom_point(size = 5) +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3, nudge_y = -0.05)  #nudge_y to seperate id from point
#Richness - Alpha diversity ####
plot_richness(ps.nema, "year","neonic", measures = "Shannon") +
  geom_text(aes(label = sampleid))

#PERMANOVA ####
#make dataframe
df = as(sample_data(ps.nema), "data.frame")
#bray-curtis distance
dis = phyloseq::distance(ps.nema,  method = "bray")
set.seed(901)
adns = adonis2(dis ~ neonic * host * year * month, df) #distance = bray
adns

#Ordination ####
pcoa = ordinate(ps.nema, method = "MDS", k = 2, try = 100, distance = "bray")
#ordinate
pcoa1 = plot_ordination(ps.nema, pcoa) 
pcoa1$layers
pcoa1$layers = pcoa1$layers[-1] #remove the original points to add the desired colors and shapes
#Fig.1 ####
#neonic,month,facet:year
#plot
fig.pcoa_neo_yr = pcoa1 + 
  theme_bw() +
  #group by neonic
  stat_ellipse(aes(fill=neonic, group=neonic), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.3, show.legend=TRUE) +
  scale_fill_manual(name="Treatement ellipses",values = c("wheat4","tomato3"),
                    labels = c(N="Control", Y="Neonicotinoid-treated")) + 
  geom_point(aes(color = neonic), size=4, alpha=0.85) +
  facet_wrap(~ year, ncol = 2,
             labeller=labeller(year = c("2016"="Soybean 2016","2017"="Corn 2017","2018"="Soybean 2018"))) + 
  #change facet font size
  theme(strip.text.x = element_text(size=16, face="bold")) +
  #change legend labels 
  scale_color_manual(name = "Treatment",values = c("cornflowerblue","mediumvioletred"),
                     labels = c(N="Control", Y="Neonicotinoid-treated")) +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.position = "right",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text( size = 16, face = "bold")) +
  guides(col = guide_legend(nrow = 3, title.position = "top"))  
fig.pcoa_neo_yr

#Richness ####
shn.rich = cbind(estimate_richness(ps.nema,measures = c('shannon','Observed','Simpson')),
                 sample_data(ps.nema))
summary(shn.rich$Shannon)
sd(shn.rich$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich$Shannon[!is.na(shn.rich$Shannon)])) #SE
summary(shn.rich$Observed)

# Wilcoxon rank-sum test (Mann-Whitney)
library(ggpubr); packageVersion("ggpubr") #‘0.3.0’
#shannon is sensitive to rare species
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "none") 
#simpson weights common species
compare_means(Simpson ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "none") 

compare_means(Observed ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "none")

#or
pairwise.wilcox.test(shn.rich$Shannon, sample_data(shn.rich)$neonic)
ggplot(shn.rich, aes(neonic, y = Shannon, alpha = 0.1)) + #alpha: color intensity)) 
  theme_bw() +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = 0.5, size = 14, label.y = 2.7) +# Add global p-value; #holm
  geom_jitter(alpha = 0.1)
#Richness / treatment ####
#control
shn.rich.ctl = cbind(estimate_richness(ps.nema.ctl,measures = c('shannon','Observed')),
                     sample_data(ps.nema.ctl))
summary(shn.rich.ctl$Shannon)
sd(shn.rich.ctl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.ctl$Shannon[!is.na(shn.rich.ctl$Shannon)])) #SE
summary(shn.rich.ctl$Observed)
sd(shn.rich.ctl$Observed, na.rm=TRUE) /  
  sqrt(length(shn.rich.ctl$Observed[!is.na(shn.rich.ctl$Observed)])) #SE
#neonic-treated
shn.rich.neo = cbind(estimate_richness(ps.nema.neo,measures = c('shannon','Observed')),
                     sample_data(ps.nema.neo))
summary(shn.rich.neo$Shannon)
sd(shn.rich.neo$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.neo$Shannon[!is.na(shn.rich.neo$Shannon)])) #SE
summary(shn.rich.neo$Observed)
sd(shn.rich.neo$Observed, na.rm=TRUE) /  
  sqrt(length(shn.rich.neo$Observed[!is.na(shn.rich.neo$Observed)])) #SE

#Nematode functions ####
#create tables of relative abundance of families for each treatment
rel_abund_neo = ps.nema.neo %>%
  tax_glom(taxrank = "Family", NArm=FALSE) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt() %>%
  group_by(Family) %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  arrange(desc(Abundance)) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100, neonic = "Y")
rel_abund_ctl = ps.nema.ctl %>%
  tax_glom(taxrank = "Family", NArm=FALSE) %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  psmelt() %>%
  group_by(Family) %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  arrange(desc(Abundance)) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100, neonic = "N") #add the rel_abund and the neonic columns

dim(left_join(rel_abund_ctl,rel_abund_neo,by = "Family"))
#put control and neonic together in one table
rel_abund_neo.ctl = rbind(rel_abund_neo, rel_abund_ctl)
View(rel_abund_neo.ctl)



