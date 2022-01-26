###########################################################################################################
# Explanatory analysis of soil samples                                                                    #
# Studying soil nematode family and trophic group variations in response to neonicotinoid seed treatment  #
# Data: Miseq-18S - L'Acadie (ACA)                                                                        #
# Mona Parizadeh - 2020-2021                                                                              #
###########################################################################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.34.0’
library(vegan); packageVersion("vegan") #‘2.5.7’
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(dplyr); packageVersion("dplyr") #‘1.0.4’
library(picante); packageVersion("picante") #'1.8.2'

# Import data #### 
setwd("../mp/aca_18s/files/")
ps.nema = readRDS("18S_nema_aca_soil1000.rds")
ps.nema

#Make families dataframe ####
#agglomerate at family level
ps.nema.fml = tax_glom(ps.nema,taxrank = "Family", NArm=FALSE)
#relative abundance
ps.ra = transform_sample_counts(ps.nema.fml, function(otu) otu/sum(otu))
#melt
nema.melt = psmelt(ps.ra) 
#sum(nema.melt$Abundance/length(levels(nema.melt$sampleid))) #verify that relative abundance of each sample equals to 1
colnames(nema.melt)
nema.melt_neo.fml = nema.melt[,-c(1,2,5:7,9:13)] #remove unnecessary columns
nema.melt_neo.fml.fac = factor(nema.melt_neo.fml$Family); levels(nema.melt_neo.fml.fac) #list of families

#Add trophic groups to dataframe ####
Trophic_group = unlist(lapply(nema.melt_neo.fml$Family, FUN = function(x){
  if(is.na(x) == TRUE){
    return("NA")
  } else if (x == "Belondiridae" || x == "Criconematidae" || x == "Merliniidae" || x == "Pratylenchidae" ||
             x == "Tylenchidae" || x == "Tylenchulidae"){
    #Belondiridae: (Yeates & Nathe: plant-feeder and omnivore)
    #check genera:
    # tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Belondiridae"]: Dorylaimellus, Oxydirus (genera, both herbivores)
    return("Herbivores") #/Plant-parasites
  } else if (x == "Anguinidae" || x == "Aphelenchidae" || x == "Aphelenchoididae" || x == "Diphtherophoridae"){
    #hyphal-feeder = fungivore?
    #Anguinidae, Aphelenchidae, Aphelenchoididae: (Yeates & Nath: hyphal-feeder; Nath: plant feeder)
    #check genera:
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Anguinidae"]: Pseudhalenchus (fungivore)
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Aphelenchidae"]: Aphelenchus, Paraphelenchus (both fungivores)
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Aphelenchoididae"]: Aphelenchoides (fungivore), Schistonchus (-)
    return("Fungivores")
  } else if (x == "Alaimidae" || x == "Cephalobidae" || x == "Diplopeltidae" || x == "Monhysteridae" ||
             x == "Neodiplogasteridae" || x == "Panagrolaimidae" || x == "Plectidae" || x == "Prismatolaimidae" ||
             x == "Rhabditidae"){
    #Neodiplogasteridae: (Yeates & Nath: bacterial-feeder, animal-predator, etc.) 
    #check genera:
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Neodiplogasteridae"]: Pristionchus (bacterivore, also predator)
    return("Bacterivores")
  } else if (x == "Mononchidae" || x == "Mylonchulidae" || x == "Trischistomatidae" || x == "Mermithidae"){
    #Mononchidae: (Yeates: animal-predator; Nath: plant-feeder, animal-predator)
    #check genera:
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Mononchidae"]: Prionchulus, Clarkus, Mononchus (all predators)
    #Trischistomatidae: (Nemaplex: generalist predators of small aquatic and soil organisms (like Tripylidae); neither Yeates nor Nath)
    #Mermithidae: insect-predator
    return("Predators")     #predacious - animal-predators
  # } else if (x == "Cyatholaimidae"){
  #   return("Algal-feeders") #unicellular eucaryote feeder,
  # } else if (x == "Mermithidae"){
  #   #Mermithidae: insect-parasite
  #   return("Insect-parasites")
  } else if (x == "Dorylaimidae" || x == "Nordiidae" || x == "Qudsianematidae" || x == "Leptonchidae" || x == "Cyatholaimidae"){
    return("Omnivores") #omnivorous
    #Nordiidae: (Yeates & Nath: omnivore, plant-feeder)
    #check genera:
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Nordiidae"]: (-)
    #Proleptonchus: (Yeates & Nath: hyphal-feeder and omnivore)
    #check genera:
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Leptonchidae"]: Proleptonchus (omnivore)
    #Qudsianematidae: (Yeates: - ; Nath: omnivore, animal-predator)
    #check genera:
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Qudsianematidae"]->(Dorylaimidae (omnivore family) in Yeates): Ecumenicus (omnivore), Labronema and Thonus (both),
    #Leptonchidae: (Yeates: hyphal-feeder, Nath: hyphal-feeder and omnivore)
    #check genera:
    #tax_table(ps.nema)[,5:6][tax_table(ps.nema)[,5] == "Leptonchidae"]: Proleptonchus (omnivore)
    #Cyatholaimidae: omnivore, algal-feeder
    } else {
    return("undefined")
  }
}))
nema.melt_neo.fml.trpc = nema.melt_neo.fml
nema.melt_neo.fml.trpc$Trophic_group <- Trophic_group
nema.melt_neo.fml.trpc

#Order relative abundance ####
#trophic groups
nema_trpc.ord = nema.melt_neo.fml.trpc %>%
  group_by(neonic,Trophic_group)  %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  arrange(dplyr::desc(Abundance)) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100)  #add the rel_abund and the neonic columns
nema_trpc.ord
nema_trpc.ord[nema_trpc.ord$neonic == "Y",]
nema_trpc.ord[nema_trpc.ord$neonic == "N",]

#Table S1 ####
#families
nema_fml.ord = nema.melt_neo.fml.trpc %>%
  group_by(neonic,Family)  %>%
  # Add abundances within each phylum 
  summarize_at("Abundance", sum)  %>%
  arrange(dplyr::desc(Abundance)) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100) %>%
  na.omit() #toral rel_abund is ~ 160 (1.NAs removed, 2.both N and Y together)

nema_fml.ord[nema_fml.ord$neonic == "Y",]
nema_fml.ord[nema_fml.ord$neonic == "N",]

#Trophic Groups ####
#remove NAs
all_trpc_noNA = nema.melt_neo.fml.trpc %>% na.omit() 
#define variables as factors
all_tp_gr = get_variable(all_trpc_noNA, "Trophic_group")
all_neo_gr = get_variable(all_trpc_noNA, "neonic")
#group
group.neo.all_tp_gr = paste(all_tp_gr,all_neo_gr, sep = "")
#order as factor
group.neo.all_tp_gr.fac = factor(group.neo.all_tp_gr); levels(group.neo.all_tp_gr.fac)

library(ggpubr)
stat.test_neo = compare_means(Abundance ~ neonic, nema.melt_neo.fml.trpc, method = "wilcox.test", paired = FALSE, 
                                  p.adjust.method = "holm", group.by = "Trophic_group") %>% na.omit() 
stat.test_neo[stat.test_neo$p.adj<0.05,] #no significance

#Fig.2A ####
labels = c("Control Bacterivores","Neonicotinoid-treated Bacterivores",
           "Control Fungivores","Neonicotinoid-treated Fungivores",
           "Control Herbivores","Neonicotinoid-treated Herbivores",
           "Control Predators", "Neonicotinoid-treated Predators",
           "Control Omnivores","Neonicotinoid-treated Omnivores")
gg.neo.trpc = nema.melt_neo.fml.trpc %>%
  na.omit() %>%
  ggplot(aes(x = Trophic_group, y = Abundance, fill = group.neo.all_tp_gr.fac, alpha = group.neo.all_tp_gr.fac)) + #alpha = factor(neonic), 
  geom_boxplot(aes(fill = group.neo.all_tp_gr.fac), #color = group.neo.tp_gr.fac,color = "black"
               position = position_dodge2(0.9,preserve = "single")) + #outlier.colour = NA,
  scale_alpha_manual(name = "Treatments & Trophic Groups", labels = labels, 
                     values = c(0.3,0.9,0.3,0.9,0.3,0.9,0.3,0.9,0.3,0.9)) +
  scale_fill_manual(name = "Treatments & Trophic Groups",
                    labels = labels,
                    values = c("cornflowerblue","cornflowerblue",
                               "darkgoldenrod1","darkgoldenrod1","indianred3","indianred3",
                               "#A29B32","#A29B32","mediumvioletred","mediumvioletred")) +
  scale_color_manual(name = "Treatments & Trophic Groups",
                     labels = labels,
                     values = c("cornflowerblue","cornflowerblue",
                                "darkgoldenrod1","darkgoldenrod1","indianred3","indianred3",
                                "#A29B32","#A29B32","mediumvioletred","mediumvioletred")) +
  scale_y_sqrt() +
  theme_bw() +
  xlab("\nTrophic Group") +
  #ylab("Relative Abundance") +
  ylab("Relative Abundance (square-root transformed)") +
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title  = element_text(size = 14, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size=14, face="bold"),legend.text=element_text(size=12),
        #strip.text = element_text(size=13, face="bold"),
        strip.background = element_rect(color = "white", size = 1)) +
  guides(fill = guide_legend(override.aes = aes(label = ""))) + #remove the "a" from legend boxes
  labs(tag = "A)") + theme(plot.tag = element_text(size = 18, face = "bold"))
gg.neo.trpc

#Families ####
#remove NAs
all_fml_trpc_noNA = nema.melt_neo.fml.trpc %>% na.omit() 
#define variables as factors
tp_gr = get_variable(all_fml_trpc_noNA, "Trophic_group")
neo_gr = get_variable(all_fml_trpc_noNA, "neonic")
#group
group.neo.tp_gr = paste(tp_gr,neo_gr, sep = "")
#order as factor
group.neo.tp_gr.fac = factor(group.neo.tp_gr)
#wilcoxon test
library(ggpubr)
stat.test_neo.fml = compare_means(Abundance ~ neonic, nema.melt_neo.fml.trpc, method = "wilcox.test", paired = FALSE, 
                                  p.adjust.method = "holm", group.by = "Family") %>% na.omit() 

stat.test <- nema.melt_neo.fml.trpc %>%
  group_by(Family) %>%
  wilcox_test(Abundance ~ neonic) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>% na.omit() 
stat.test[stat.test$p.adj<0.05,]

#Table S1 ####
stat.test_neo.fml[stat.test_neo.fml$p.adj<0.05,]

#Fig.2B####
gg.neo.fml.trpc = nema.melt_neo.fml.trpc %>%
  na.omit() %>%
  ggplot(aes(x = Family, y = Abundance, fill = group.neo.tp_gr.fac, alpha = group.neo.tp_gr.fac)) + #alpha = factor(neonic), 
  geom_boxplot(aes(fill = group.neo.tp_gr.fac), #color = group.neo.tp_gr.fac,color = "black"
               position = position_dodge2(0.9,preserve = "single")) + #outlier.colour = NA,
  scale_alpha_manual(name = "Treatments & Trophic Groups", labels = labels, 
                     values = c(0.3,0.9,0.3,0.9,0.3,0.9,0.3,0.9,0.3,0.9)) +
  facet_grid(~ Trophic_group, scales='free_x', space = "free") +
  scale_fill_manual(name = "Treatments & Trophic Groups",
                    labels = labels,
                    values = c("cornflowerblue","cornflowerblue",
                               "darkgoldenrod1","darkgoldenrod1","indianred3","indianred3",
                               "#A29B32","#A29B32","mediumvioletred","mediumvioletred")) +
  scale_color_manual(name = "Treatments & Trophic Groups",
                    labels = labels,
                    values = c("cornflowerblue","cornflowerblue",
                               "darkgoldenrod1","darkgoldenrod1","indianred3","indianred3",
                               "#A29B32","#A29B32","mediumvioletred","mediumvioletred")) +
  scale_y_sqrt() +
  theme_bw() +
  xlab("Family") +
  #ylab("Relative Abundance") +
  ylab("Relative Abundance (square-root transformed)") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title  = element_text(size = 14, face = "bold"),
        legend.position = "none",
        legend.title = element_text(size=14, face="bold"),legend.text=element_text(size=12),
        strip.text = element_text(size=13, face="bold"),
        strip.background = element_rect(color = "white", size = 1)) +
  guides(fill = guide_legend(override.aes = aes(label = ""))) + #remove the "a" from legend boxes
  labs(tag = "B)") + theme(plot.tag = element_text(size = 18, face = "bold"))
gg.neo.fml.trpc

