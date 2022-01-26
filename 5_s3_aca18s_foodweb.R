###############################################################
# Faunal analysis                                             #
# Data: Miseq-18S - L'Acadie (ACA)                            #
# Mona Parizadeh - 2020-2021                                  #
###############################################################

# Load libraries ####
library(ggplot2); packageVersion("ggplot2") #‘3.3.3’
library(tidyverse); packageVersion("tidyverse") #‘1.3.0’

#Import the table created by NINJA based on food web analysis
setwd("../mp/aca_18s/files/")
fw <- read.csv("NINJA_data.csv", sep=";"); dim(fw)

#subset data
fw$year = factor(as.character(fw$year))
fw_y <- subset.data.frame(fw, dataset = "year")
fw_m <- subset.data.frame(fw, dataset = "month")

#plot
ggplot(fw_y, aes(x=Structure_Index, y=Enrichment_Index, colour=year, shape=treatment)) +
  geom_point(size=5) +
  geom_hline(yintercept=50) +
  geom_vline(xintercept=50) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = c("#308005", "#7f5d00", "#a81a0d")) +
   labs(
    title = "",
    y="Enrichment Index",
    x="Structure Index",
    fill="Year"
  ) +
  theme_classic() 

#Fig.S2 ####
#% variables
group.yr.trt = paste(fw_y$year,fw_y$treatment, sep = "");group.yr.trt
#plot
fd.wb_yr = ggplot(fw_y, aes(x=Structure_Index, y=Enrichment_Index, colour=group.yr.trt, shape=group.yr.trt)) +
  #facet_wrap(facet = "month", nrow=2) +
  theme_bw() +
  geom_point(size=8) +
  geom_hline(yintercept=50) +
  geom_vline(xintercept=50) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_shape_manual(name="Treatement & Year", values = c(19,19,17,17,15,15),
                     labels = c("2016Control"="Control 2016", "2016Neonic"="Neonicotinoid-treated 2016",
                                "2017Control"="Control 2017", "2017Neonic"="Neonicotinoid-treated 2017",
                                "2018Control"="Control 2018", "2018Neonic"="Neonicotinoid-treated 2018")) +
  scale_colour_manual(name="Treatement & Year", values = c("cornflowerblue","mediumvioletred",
                                                           "cornflowerblue","mediumvioletred",
                                                           "cornflowerblue","mediumvioletred"),
                      labels = c("2016Control"="Control 2016", "2016Neonic"="Neonicotinoid-treated 2016",
                                 "2017Control"="Control 2017", "2017Neonic"="Neonicotinoid-treated 2017",
                                 "2018Control"="Control 2018", "2018Neonic"="Neonicotinoid-treated 2018")) +
  labs(title = "", y="Enrichment Index", x="Structure Index") +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.position = "right",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text( size = 16, face = "bold"))  
fd.wb_yr 

#plot
ggplot(fw_m, aes(x=Structure_Index, y=Enrichment_Index, colour=year, shape=month)) + 
  facet_wrap(facet = "treatment", nrow=2) +
  geom_point(size=5) +
  geom_hline(yintercept=50) +
  geom_vline(xintercept=50) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = c("#308005", "#7f5d00", "#a81a0d")) +
  labs(
    title = "",
    y="Enrichment Index",
    x="Structure Index",
    fill="Year"
  ) +
  theme_classic() 
  
  Stat <- fw %>% 
  group_by(treatment) %>%
  dplyr::summarise(across(
    .cols = is.numeric, 
    .fns = list(Mean = mean, SD = sd), na.rm = TRUE, 
    .names = "{col}_{fn}"
  ))
View(Stat)

#Maturity index
anova_MI <- aov(Maturity_Index ~ treatment*year/month, data = fw)
summary(anova_MI)

#Omnivores (J'ai test? les autres et pas de diff.)
anova_Omni <- aov(Omnivore_footprint ~ treatment*year/month, data = fw)
summary(anova_Omni)


