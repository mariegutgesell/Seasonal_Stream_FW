##Code: Chapter 2 M3 2x-week pan trap and emergent tent exploratory analysis 
##By: Marie Gutgesell
##Date: January 8, 2021
##Updated: May 5, 2022

##Where left off on May 5 -- want to go through data files/plots and sort out any missing data, identifying benthic data that overlaps w/ e-fishing dates
##How do i deal w/ aquatic emergence falling back in? -- not really terrestrial subsidy.. 

##Set Working Directory
setwd("~/Desktop/Seasonal Stream Food Webs/Data/M3 Insect Data/")


##Load Packages 
library(ggplot2)
library(corrplot) 
library(Hmisc)
library ( ggbiplot )
library (dplyr)
library(tidyr)

library(tidyverse)
library(lubridate)
library(dplyr)
library(ggpmisc)
library(reshape)

##Import csv file containing ID data from emergent, pan and benthic traps 
emergent_data <- read.csv ("Emergent Trap ID_M3 2018_MASTERCSV.csv")

pan_data <- read.csv ("Pan Trap ID_M3 2018_MASTER.csv")

benthic_data <- read.csv ("Benthic Monitoring ID Sheet_MASTER.csv")

##Emergent Trap Exploratory Time Series 
emergent_data$Abundance <- as.numeric(as.character(emergent_data$Abundance))

emergent_data <- emergent_data %>%
  na.omit()

##Create subset of only 1st replicate data, unless data was not recorded, then using second replicate  
emergent_R1 <- emergent_data %>%
  filter(Replicate == "R_1" | X == "USE")

hc_ET <- emergent_R1 %>%
  filter(Site_Code == "HC")
as.character(hc_ET$Date)
hc_seasonal_abundance <- aggregate(hc_ET$Abundance, by = list(hc_ET$Date), FUN = sum, ra.rm = TRUE)
colnames(hc_seasonal_abundance)[1] <- "Date"
hc_seasonal_ET <- ggplot(hc_seasonal_abundance, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(hc_seasonal_ET)

ep4_ET <- emergent_R1 %>%
  filter(Site_Code == "EP4")
ep4_seasonal_abundance <- aggregate(ep4_ET$Abundance, by = list(ep4_ET$Date), FUN = sum, ra.rm = TRUE)
colnames(ep4_seasonal_abundance)[1] <- "Date"
ep4_seasonal_ET <- ggplot(ep4_seasonal_abundance, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(ep4_seasonal_ET)

at_ET <- emergent_R1 %>%
  filter(Site_Code == "AT")
at_seasonal_abundance <- aggregate(at_ET$Abundance, by = list(at_ET$Date), FUN = sum, ra.rm = TRUE)
colnames(at_seasonal_abundance)[1] <- "Date"
at_seasonal_ET <- ggplot(at_seasonal_abundance, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(at_seasonal_ET)

seasonal_abundance_total <- left_join(hc_seasonal_abundance, ep4_seasonal_abundance, by = "Date") %>%
  left_join(at_seasonal_abundance)
colnames(seasonal_abundance_total)[2] <- "HC"
colnames(seasonal_abundance_total)[3] <- "EP4"
colnames(seasonal_abundance_total)[4] <- "AT"
seasonal_abundance_total <- melt(seasonal_abundance_total, id.vars = "Date")

seasonal_ET <- ggplot(seasonal_abundance_total, aes(
  x = Date,
  y = value,
  group = variable, 
  col = variable
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(seasonal_ET)


###Pan Trap 
pan_data$Abundance <- as.numeric(as.character(pan_data$Abundance))

pan_data_rep1 <- pan_data %>%
  filter(Replicate == "Rep_1")

hc_PT <- pan_data_rep1 %>%
  filter(Site_ID == "HC")
as.character(hc_PT$Date)
hc_seasonal_abundance_pt <- aggregate(hc_PT$Abundance, by = list(hc_PT$Date), FUN = sum, ra.rm = TRUE)
colnames(hc_seasonal_abundance_pt)[1] <- "Date"
hc_seasonal_PT <- ggplot(hc_seasonal_abundance_pt, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(hc_seasonal_PT)

ep4_PT <- pan_data_rep1 %>%
  filter(Site_ID == "EP4")
ep4_seasonal_abundance_pt <- aggregate(ep4_PT$Abundance, by = list(ep4_PT$Date), FUN = sum, ra.rm = TRUE)
colnames(ep4_seasonal_abundance_pt)[1] <- "Date"
ep4_seasonal_PT <- ggplot(ep4_seasonal_abundance_pt, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(ep4_seasonal_PT)

at_PT <- pan_data_rep1 %>%
  filter(Site_ID == "AT")
at_seasonal_abundance_pt <- aggregate(at_PT$Abundance, by = list(at_PT$Date), FUN = sum, ra.rm = TRUE)
colnames(at_seasonal_abundance_pt)[1] <- "Date"
at_seasonal_PT <- ggplot(at_seasonal_abundance_pt, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(at_seasonal_PT)

pt_seasonal_abundance_total <- left_join(hc_seasonal_abundance_pt, ep4_seasonal_abundance_pt, by = "Date") %>%
  left_join(at_seasonal_abundance_pt)
colnames(pt_seasonal_abundance_total)[2] <- "HC"
colnames(pt_seasonal_abundance_total)[3] <- "EP4"
colnames(pt_seasonal_abundance_total)[4] <- "AT"
pt_seasonal_abundance_total <- melt(pt_seasonal_abundance_total, id.vars = "Date")

seasonal_PT <- ggplot(pt_seasonal_abundance_total, aes(
  x = Date,
  y = value,
  group = variable, 
  col = variable
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(seasonal_PT)

###Benthic 
benthic_data$abundance <- as.numeric(as.character(benthic_data$abundance))

benthic_data_rep1 <- benthic_data %>%
  filter(Replicate == "Rep_1")

hc_BM <- benthic_data_rep1 %>%
  filter(Site_ID == "HC")
as.character(hc_BM$Date)
hc_seasonal_abundance_bm <- aggregate(hc_BM$abundance, by = list(hc_BM$Date), FUN = sum, ra.rm = TRUE)
colnames(hc_seasonal_abundance_bm)[1] <- "Date"
hc_seasonal_BM <- ggplot(hc_seasonal_abundance_bm, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(hc_seasonal_BM)

ep4_BM <- benthic_data_rep1 %>%
  filter(Site_ID == "EP4")
ep4_seasonal_abundance_bm <- aggregate(ep4_BM$abundance, by = list(ep4_BM$Date), FUN = sum, ra.rm = TRUE)
colnames(ep4_seasonal_abundance_bm)[1] <- "Date"
ep4_seasonal_BM <- ggplot(ep4_seasonal_abundance_bm, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(ep4_seasonal_BM)

at_BM <- benthic_data_rep1 %>%
  filter(Site_ID == "AT")
at_seasonal_abundance_bm <- aggregate(at_BM$abundance, by = list(at_BM$Date), FUN = sum, ra.rm = TRUE)
colnames(at_seasonal_abundance_bm)[1] <- "Date"
at_seasonal_BM <- ggplot(at_seasonal_abundance_bm, aes(
  x = Date,
  y = x,
  group = 1
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(at_seasonal_BM)

bm_seasonal_abundance_total <- left_join(hc_seasonal_abundance_bm, ep4_seasonal_abundance_bm, by = "Date") %>%
  left_join(at_seasonal_abundance_bm)
colnames(bm_seasonal_abundance_total)[2] <- "HC"
colnames(bm_seasonal_abundance_total)[3] <- "EP4"
colnames(bm_seasonal_abundance_total)[4] <- "AT"
bm_seasonal_abundance_total <- melt(bm_seasonal_abundance_total, id.vars = "Date")

seasonal_BM <- ggplot(bm_seasonal_abundance_total, aes(
  x = Date,
  y = value,
  group = variable, 
  col = variable
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(seasonal_BM)

###Comparing all 3 traps within sites
##HawkCliff - Conservation Area

hc_seasonal_abundance_total <- left_join(hc_seasonal_abundance_bm, hc_seasonal_abundance_pt, by = "Date") %>%
  left_join(hc_seasonal_abundance)
colnames(hc_seasonal_abundance_total)[2] <- "In-stream Insect Abundance"
colnames(hc_seasonal_abundance_total)[3] <- "Terrestrial In-Fall Insect Abundance"
colnames(hc_seasonal_abundance_total)[4] <- "Emergent Insect Abundance"
hc_seasonal_abundance_total <- melt(hc_seasonal_abundance_total, id.vars = "Date")

seasonal_HC <- ggplot(hc_seasonal_abundance_total, aes(
  x = Date,
  y = value,
  group = variable, 
  col = variable
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(seasonal_HC)

##EP4 - Mid-Impact Farm

ep4_seasonal_abundance_total <- left_join(ep4_seasonal_abundance_bm, ep4_seasonal_abundance_pt, by = "Date") %>%
  left_join(ep4_seasonal_abundance)
colnames(ep4_seasonal_abundance_total)[2] <- "In-stream Insect Abundance"
colnames(ep4_seasonal_abundance_total)[3] <- "Terrestrial In-Fall Insect Abundance"
colnames(ep4_seasonal_abundance_total)[4] <- "Emergent Insect Abundance"
ep4_seasonal_abundance_total <- melt(ep4_seasonal_abundance_total, id.vars = "Date")

seasonal_EP4 <- ggplot(ep4_seasonal_abundance_total, aes(
  x = Date,
  y = value,
  group = variable, 
  col = variable
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(seasonal_EP4)

##Taylor - Conventional Farm

at_seasonal_abundance_total <- left_join(at_seasonal_abundance_bm, at_seasonal_abundance_pt, by = "Date") %>%
  left_join(at_seasonal_abundance)
colnames(at_seasonal_abundance_total)[2] <- "In-stream Insect Abundance"
colnames(at_seasonal_abundance_total)[3] <- "Terrestrial In-Fall Insect Abundance"
colnames(at_seasonal_abundance_total)[4] <- "Emergent Insect Abundance"
at_seasonal_abundance_total <- melt(at_seasonal_abundance_total, id.vars = "Date")

seasonal_AT <- ggplot(at_seasonal_abundance_total, aes(
  x = Date,
  y = value,
  group = variable, 
  col = variable
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(seasonal_AT)
