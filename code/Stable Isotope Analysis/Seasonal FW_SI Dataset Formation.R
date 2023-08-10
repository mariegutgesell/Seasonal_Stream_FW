##Seasonal Stream Food Webs_SI Dataset Formation
##By: Marie Gutgesell
##Date: May 4, 2022
##Updated Code: 

##Set Working Directory
setwd("~/Desktop/Seasonal Stream Food Webs/Data/Stable Isotope Data")

##Load Packages -- 
library(tidyverse)
library (dplyr)
library(tidyr)
library(lubridate)
library(readxl)

##Import dH and CN data
dH_HC <- read.csv("HC_dH_clean.csv", stringsAsFactors = FALSE)
names(dH_HC)[2] <- "Sample_ID"
names(dH_HC)[5] <- "d2H"
names(dH_HC)[9] <- "Species"
dH_HC <- dH_HC %>%
  select(Sample_ID, d2H, Species)

dH_EP4 <- read.csv("EP4_dH_clean.csv", stringsAsFactors = FALSE)
names(dH_EP4)[2] <- "Sample_ID"
names(dH_EP4)[5] <- "d2H"
names(dH_EP4)[9] <- "Species"
dH_EP4 <- dH_EP4 %>%
  select(Sample_ID, d2H, Species)


dH_AT <- read.csv("AT_dH_clean.csv", stringsAsFactors = FALSE)
names(dH_AT)[2] <- "Sample_ID"
names(dH_AT)[5] <- "d2H"
names(dH_AT)[9] <- "Species"
dH_AT <- dH_AT %>%
  select(Sample_ID, d2H, Species)


dH_all <- rbind(dH_HC, dH_EP4, dH_AT) %>%
  group_by(Sample_ID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  na.omit()


##
dCN_all <- read.csv("SI data MASTER Marie Gutgesell 2019_clean.csv", stringsAsFactors = FALSE)
names(dCN_all)[1] <- "Sample_ID"
names(dCN_all)[2] <- "d13C"
names(dCN_all)[3] <- "per_C"
names(dCN_all)[4] <- "d15N"
names(dCN_all)[5] <- "per_N"
names(dCN_all)[6] <- "CN_ratio"
dCN_all <- dCN_all %>%
  select(Sample_ID, d13C, per_C, per_N, d15N, CN_ratio) %>%
  na.omit() %>%
  group_by(Sample_ID) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

SI_all <- left_join(dCN_all, dH_all, by = "Sample_ID")

na_test <- SI_all[is.na(SI_all$d2H),]
##15/18 are other fish species (i.e., not CC) that were not sent for dH analysis
##other missing: AT1-18-AQS-03 -- missing in dH dataset, looks like not done, maybe see if sent?
##other missing: AT2-19-CC-35 -- were dropped, not clear why not analyzed
##other missing: AT2-19-CC-37 -- were dropped, not clear why not analyzed



SI_all <- SI_all %>%
  mutate(Site = case_when(
    startsWith(Sample_ID, "HC") ~ "HawkCliff",
    startsWith(Sample_ID, "AT") ~ "Taylor",
    startsWith(Sample_ID, "EP4") ~ "Perl_4",
  ))

SI_all <- SI_all %>%
  mutate(Season = case_when(
    startsWith(Sample_ID, "HC1") ~ "Spring",
    startsWith(Sample_ID, "AT1") ~ "Spring",
    startsWith(Sample_ID, "EP41") ~ "Spring",
    startsWith(Sample_ID, "HC2") ~ "Summer",
    startsWith(Sample_ID, "AT2") ~ "Summer",
    startsWith(Sample_ID, "EP42") ~ "Summer",
    startsWith(Sample_ID, "HC3") ~ "Fall",
    startsWith(Sample_ID, "AT3") ~ "Fall",
    startsWith(Sample_ID, "EP43") ~ "Fall",
  ))


SI_all <- SI_all %>%
  mutate(Species = case_when(
    grepl("CC", Sample_ID) ~ "Creek Chub",
    grepl("AL", Sample_ID) ~ "Algae",
    grepl("DT", Sample_ID) ~ "Detritus",
    grepl("DET", Sample_ID) ~ "Detritus",
    grepl("AQS", Sample_ID) ~ "Aquatic Snail",
    grepl("SC", Sample_ID) ~ "Goldenrod",
    grepl("TS", Sample_ID) ~ "Terrestrial Snail",
    grepl("GH", Sample_ID) ~ "Grasshopper",
    grepl("SN", Sample_ID) ~ "Aquatic Snail",
    grepl("CAT", Sample_ID) ~ "Caterpillar",
    grepl("WBM", Sample_ID) ~ "Waterboatman",
    grepl("EP", Sample_ID) ~ "Mayfly",
    grepl("HE", Sample_ID) ~ "Terrestrial Snail",
    grepl("SU", Sample_ID) ~ "Terrestrial Snail",
    grepl("SNAIL", Sample_ID) ~ "Terrestrial Snail",
    grepl("EGS", Sample_ID) ~ "Terrestrial Snail",
    grepl("CL", Sample_ID) ~ "Clam",
    grepl("CA", Sample_ID) ~ "Caterpillar",
    grepl("TR", Sample_ID) ~ "Caddisfly",
    grepl("AM", Sample_ID) ~ "Amphipod",
    grepl("CA", Sample_ID) ~ "Caterpillar",
  ))

SI_all <- SI_all %>%
  mutate(Trophic_Group = case_when(
    grepl("CC", Sample_ID) ~ "Creek Chub",
    grepl("AL", Sample_ID) ~ "Aquatic Primary Producer",
    grepl("DT", Sample_ID) ~ "Terrestrial Primary Producer",
    grepl("DET", Sample_ID) ~ "Terrestrial Primary Producer",
    grepl("AQS", Sample_ID) ~ "Aquatic Invert",
    grepl("SC", Sample_ID) ~ "Terrestrial Primary Producer",
    grepl("TS", Sample_ID) ~ "Terrestrial Invert",
    grepl("GH", Sample_ID) ~ "Terrestrial Invert",
    grepl("SN", Sample_ID) ~ "Aquatic Invert",
    grepl("CAT", Sample_ID) ~ "Terrestrial Invert",
    grepl("WBM", Sample_ID) ~ "Aquatic Invert",
    grepl("EP", Sample_ID) ~ "Aquatic Invert",
    grepl("HE", Sample_ID) ~ "Terrestrial Invert",
    grepl("SU", Sample_ID) ~ "Terrestrial Invert",
    grepl("SNAIL", Sample_ID) ~ "Terrestrial Invert",
    grepl("EGS", Sample_ID) ~ "Terrestrial Invert",
    grepl("CL", Sample_ID) ~ "Aquatic Invert",
    grepl("CA", Sample_ID) ~ "Terrestrial Invert",
    grepl("TR", Sample_ID) ~ "Aquatic Invert",
    grepl("AM", Sample_ID) ~ "Aquatic Invert",
    grepl("CA", Sample_ID) ~ "Terrestrial Invert",
  ))

####select only CC, algae and det
SI_cc_al_dt <- SI_all %>%
  na.omit() %>%
  filter(Species == "Creek Chub" | Species == "Algae" | Species == "Detritus")
  
  
rm(dCN_all, dH_all, dH_AT, dH_EP4, dH_HC, na_test)






