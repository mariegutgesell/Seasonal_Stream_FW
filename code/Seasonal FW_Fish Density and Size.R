##Seasonal Stream Food Webs_Fish Density & Size Plots
##By: Marie Gutgesell
##Date: May 5, 2022
##Updated Code: Oct 24, 2023



library(readxl)
library(dbplyr)
library(tidyverse)
fish_meta <- read.csv("data/Fish_Collection_Data/2018_Fish_Marie_CollectionData.csv")

##okay, so somehow in cleaning the datasheet there are a whole bunch of new bulk samples.. for hc fall, why?? need to go pick up raw data sheets, not sure wtf is happening 

##Looking at total fish density per site across season

##there are 2 unknown #s for bulk samples, need to filter those out first
fish_meta <- fish_meta %>%
  filter(Number_of_Individuals != "unknown")

fish_meta$Number_of_Individuals <- as.numeric(fish_meta$Number_of_Individuals)
fish_meta$Season <- ordered(fish_meta$Season, 
                              levels = c("Spring", "Summer", "Fall"))

##total fish
fish_density_sp <- aggregate(Number_of_Individuals ~ Season + Site + Species, fish_meta, sum)
  

fish_density_all <- aggregate(Number_of_Individuals ~ Season + Site, fish_meta, sum)

hc_fall <- fish_meta %>%
  filter(Season == "Fall") %>%
  filter(Site == "Hawk Cliff") %>%
  summarise(vars(Number_of_Individuals), sum)

aggregate(hc_fall$Number_of_Individuals, sum)

fish_density_all <- fish_density_all %>%
  mutate(Site_Code = case_when(
    startsWith(Site, "H") ~ "Low Impact",
    startsWith(Site, "T") ~ "High Impact", 
    startsWith(Site, "P") ~ "Mid Impact"
  ) )

fish_density_all$Site_Code <- ordered(fish_density_all$Site_Code,
                                          levels = c("Low Impact", "Mid Impact", "High Impact"))

group_colours <- c(Spring = "chartreuse4", Summer = "darkgoldenrod1", Fall = "cornflowerblue")

fish_density_all_plot <- ggplot(fish_density_all, aes(x = Season, y = Number_of_Individuals, fill = Site_Code)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
   theme_classic() +
  facet_wrap(~Site_Code) +
  ylab("Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 

fish_density_all_plot

fish_density_sp_plot <- ggplot(fish_density_sp, aes(x= Season, y = Number_of_Individuals, fill = Species)) +
  geom_col(position = "stack") +
  facet_wrap(~Site) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "right", text = element_text(family = "Times New Roman")) 


fish_density_sp_plot


