##Fish Biomass standardized to m2 

library(readxl)
library(dbplyr)
library(lubridate)
library(tidyverse)

##Calculate Reach Area
sm<- read.csv("data/Stream Condition Data/2x-weekly Stream Conditions_cleaned_MASTERCSV.csv", na.strings=c("","NA"))

sm$newdate <-strptime(as.character(sm$Date), "%d/%m/%Y")
sm$txtdate <- format(sm$newdate, "%y-%m-%d")

sm$Wetted_width_cm <- as.numeric(sm$Wetted_width_cm)

sm_sampling <- sm %>%
  select(newdate, txtdate, Site_Code, Water_depth_cm, Wetted_width_cm, Hydraulic_head_cm) %>%
  filter(newdate %in% c("2018-05-28", "2018-07-30", "2018-10-22")) %>%
  mutate(reach_length_m = 40) %>%
  mutate(reach_area = reach_length_m*(Wetted_width_cm/100)) %>%
  mutate(Season = case_when(
    newdate == "2018-05-28" ~ "Spring",
    newdate == "2018-07-30" ~ "Summer",
    newdate == "2018-10-22" ~ "Fall"
  ))


##Calculate Fish Biomass
fish_meta <- read.csv("data/Fish_Collection_Data/2018_Fish_Marie_CollectionData.csv")

fish_meta$Weight_g <- as.numeric(fish_meta$Weight_g)
bm <- fish_meta %>%
  select(Season, Site, Species, Weight_g, Number_of_Individuals) %>%
  filter(!is.na(Weight_g)) %>%
  group_by(Season, Site) %>%
  summarise_at(vars(Weight_g, Number_of_Individuals), list(total = sum, mean = mean, sd = sd)) %>%
  mutate(Site_Code = case_when(
    Site == "Hawk Cliff" ~ "HC",
    Site == "Perl 4" ~ "EP4",
    Site == "Taylor" ~ "AT"
  ))



##looking at size of fish for SI tissue taken
si_cc <- fish_meta %>%
  filter(SI_tissue_taken == "Y" & Species == "CC")

ggplot(si_cc, aes(x = Site, y = Weight_g, color = Season)) +
  geom_boxplot()

si$Fork_Length_per_individual_mm <- as.numeric(si$Fork_Length_per_individual_mm)
ggplot(si, aes(x = Site, y = Fork_Length_per_individual_mm, color = Season)) +
  geom_boxplot()

##join together
df <- left_join(bm, sm_sampling, by = c("Season", "Site_Code")) %>%
  mutate(biomass_m2 = Weight_g_total/reach_area)

df$Season <- ordered(df$Season, 
                            levels = c("Spring", "Summer", "Fall"))
df$Site_Code <- ordered(df$Site_Code,
                        levels = c("HC", "EP4", "AT"))

bm_plot <- ggplot(df, aes(x = Season, y = biomass_m2, fill = Site_Code)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic() +
  facet_wrap(~Site_Code) +
  ylab("Biomass (g/m2)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 

bm_plot


bm_plot_2 <- ggplot(df, aes(x = Season, y = Weight_g_total, fill = Site_Code)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic() +
  facet_wrap(~Site_Code) +
  ylab("Total Biomass (g)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 

bm_plot_2

ggplot(df, aes(x = Season, y = Number_of_Individuals_total, fill = Site_Code)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic() +
  facet_wrap(~Site_Code) +
  ylab("Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 





###
hc <- fish_meta %>% filter(Site == "Hawk Cliff")
hist(hc$Weight_g)
  
ep <- fish_meta %>% filter(Site == "Perl 4")
hist(ep$Weight_g)

at <- fish_meta %>% filter(Site == "Taylor")
hist(at$Weight_g)





##look at size spectrum rather than biomass pyramid, and look at body size of fish used for isotope data 


test <- fish_meta %>%
  filter(Site == "Taylor" & Season == "Fall")
##Looking at total fish biomass per site across season, standardized to m2 of sampling reach 

##there are 2 unknown #s for bulk samples, need to filter those out first
fish_meta <- fish_meta %>%
  filter(Number_of_Individuals != "unknown")

