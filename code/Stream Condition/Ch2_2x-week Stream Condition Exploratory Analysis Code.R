##Code: Chapter 2 M3 2x-week stream condition exploratory analysis 
##By: Marie Gutgesell
##Date: January 8, 2021

##Set Working Directory
setwd("~/Desktop/Seasonal Stream Food Webs/Data/Stream Condition Data/")


##Load Packages 
library(ggplot2)
library(tidyverse)
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

##Import csv file containing cleaned stream conditions data  
data <- read.csv("2x-weekly Stream Conditions_cleaned_MASTERCSV.csv", na.strings=c("","NA"))

data$newdate <-strptime(as.character(data$Date), "%d/%m/%Y")
data$txtdate <- format(data$newdate, "%y-%m-%d")

algae <- c("Algae_floating", "Algae_filaments", "Algae_attached", "Algae_slim_crusts")
data$algae_total <- rowSums(data[,algae])

macrophyte <- c("Macrophyte_Emergent", "Macrophyte_rooted_floating", "Macrophyte_submerged", "Macrophyte_free_floating")
data$macrophyte_total <- rowSums(data[,macrophyte])

organic_matter <- c("Decaying_OM_substrate", "Decaying_OM_floating")
data$Decaying_OM_substrate <- as.numeric(data$Decaying_OM_substrate)
data$Decaying_OM_floating <- as.numeric(data$Decaying_OM_floating)

data$OM_total <- rowSums(data[,organic_matter])

woody_debris <- c("Woody_Debris_substrate", "Woody_Debris_floating")
data$Woody_Debris_floating <- as.numeric(data$Woody_Debris_floating)
data$Woody_Debris_substrate <- as.numeric(data$Woody_Debris_substrate)
data$WD_total <- rowSums(data[,woody_debris])

##Riparian Canopy Cover

hc <- data %>%
  filter(Name == "Hawk Cliff")

canopy_cover <- ggplot(data, aes(
  x = txtdate,
  y = Canopy_Cover_avg,
  group = Site_Code,
  col = Site_Code
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(canopy_cover)

canopy_cover_hc <- ggplot(hc, aes(
  x = txtdate,
  y = Canopy_Cover_avg,
  group = Name
)) +
  geom_line()+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(canopy_cover_hc)

##Hydraulic head -- need to convert this to discharge?
hydraulic_head<- ggplot(data, aes(
  x = txtdate,
  y = Hydraulic_head_cm,
  group = Site_Code,
  col = Site_Code
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(hydraulic_head)

#Add together totals, find ways to visual data 
##algae
algae_total_seasonal<- ggplot(data, aes(
  x = txtdate,
  y = algae_total,
  group = Site_Code,
  col = Site_Code
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(algae_total_seasonal)
##macrophyte
macrophyte_total_seasonal<- ggplot(data, aes(
  x = txtdate,
  y = macrophyte_total,
  group = Site_Code,
  col = Site_Code
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(macrophyte_total_seasonal)
##woody debris
woodydebris_total_seasonal<- ggplot(data, aes(
  x = txtdate,
  y = WD_total,
  group = Site_Code,
  col = Site_Code
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(woodydebris_total_seasonal)
##organic matter
OM_total_seasonal<- ggplot(data, aes(
  x = txtdate,
  y = OM_total,
  group = Site_Code,
  col = Site_Code
)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(OM_total_seasonal)
