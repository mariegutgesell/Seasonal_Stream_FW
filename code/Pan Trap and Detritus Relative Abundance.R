##Seasonal Stream Food Webs -- Pan trap and organic matter/woody debris (terrestrial subsidy)
##By: Marie Gutgesell
##Date: June 9, 2022


##Set Working Directory
setwd("~/Desktop/Seasonal Stream Food Webs/Data/")

##Load Packages 
library(ggplot2)
library(Hmisc)
library (dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
library(ggpmisc)
library(reshape)
library(stringr)
library(ggpubr)

pan_data <- read.csv("M3 Insect Data/Pan Trap ID_M3 2018_MASTER.csv") %>%
  filter(Replicate == "Rep_1") %>%
  filter(Aq_Terr == "Terrestrial")
pan_data$newdate <-strptime(as.character(pan_data$Date), "%Y-%m-%d")

##Import csv file containing cleaned stream conditions data  
data <- read.csv("Stream Condition Data/2x-weekly Stream Conditions_cleaned_MASTERCSV.csv", na.strings=c("","NA"))

data$newdate <-strptime(as.character(data$Date), "%d/%m/%Y")
data$newdate <- as.Date(data$newdate)

##Calculate total woody debris/organic matter 
detritus <- c("Woody_Debris_substrate", "Woody_Debris_floating", "Decaying_OM_substrate", "Decaying_OM_floating")
data[detritus] <- sapply(data[detritus],as.numeric)
data$detritus_total <- rowSums(data[,detritus])



##Pan-Trap Data --------
##Subsampling date across range of 5 sampling points before e-fishing took place -- NOTE: 4 samples still to be ID: July 30, all 3 sites and Oct 15, 2018 for Taylor
data_spring <- pan_data[pan_data$newdate >= "2018-05-14" & pan_data$newdate <= "2018-05-28",]
data_summer <- pan_data[pan_data$newdate >= "2018-07-16" & pan_data$newdate <= "2018-07-30",]
data_fall <- pan_data[pan_data$newdate >= "2018-10-08" & pan_data$newdate <= "2018-10-22",]

pan_data_subset <- rbind(data_spring, data_summer, data_fall)
rm(data_spring, data_summer, data_fall)
pan_data_subset$Abundance <- as.numeric(pan_data_subset$Abundance) ##need to go through and sort out NAs

##Calculate total abundance of invertebrates
pan_data_subset <- na.omit(pan_data_subset)

pan_total_abundance <-  aggregate(pan_data_subset$Abundance, by = list(Site_ID = pan_data_subset$Site_ID, newdate = pan_data_subset$newdate), FUN = sum)
pan_total_abundance$newdate <- as.character(pan_total_abundance$newdate)
pan_total_abundance <- pan_total_abundance %>%
  mutate(Season = case_when(
    startsWith(newdate, "2018-05") ~ "Spring",
    startsWith(newdate, "2018-07") ~ "Summer", 
    startsWith(newdate, "2018-10") ~ "Fall"
  ) )
pan_total_abundance$Site_ID <- ordered(pan_total_abundance$Site_ID,
                                       levels = c("HC", "EP4", "AT"))
pan_total_abundance$Site_ID <- as.character(pan_total_abundance$Site_ID)
pan_total_abundance <- pan_total_abundance %>%
  mutate(Site = case_when(
    startsWith(Site_ID, "HC") ~ "Low Impact",
    startsWith(Site_ID, "EP") ~ "Mid Impact", 
    startsWith(Site_ID, "AT") ~ "High Impact"
  ) )

pan_total_abundance$Season <- ordered(pan_total_abundance$Season,
                                      levels = c("Spring", "Summer", "Fall"))

pan_total_abundance$Site <- ordered(pan_total_abundance$Site,
                                    levels = c("Low Impact", "Mid Impact", "High Impact"))

seasonal_pan_plot <- ggplot(pan_total_abundance, aes(x = Season, y = x, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Abundance") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 

seasonal_pan_plot


##Detritus

seasonal_det_plot <- ggplot(data_subset, aes(x = Season, y = detritus_total, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Relative Abundance Score") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

seasonal_det_plot
