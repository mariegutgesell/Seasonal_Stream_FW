##Stomach Content Data for Seasonal Stream Food Webs -- Creating summary figure



##Load Packages -- 
library(tidyverse)
library (dplyr)
library(tidyr)
library(lubridate)
library(readxl)
library(ggpubr)

fish_meta <- read.csv("data/Fish_Collection_Data/2018_Fish_Marie_CollectionData.csv")

fish_meta_cc <- fish_meta %>%
  filter(Species == "CC") %>%
  filter(Kept == "Y") %>%
  select(Season, Site, Species, Sample_ID, Stomach_Contents)

taylor_spring_sc <- read_excel("data/Fish_Collection_Data/Taylor_Spring_Minnow Trap_Stomach Contents.xlsx") %>%
  select(Season, Site, Species, Sample_ID, Stomach_Contents)

fish_meta_cc <- rbind(fish_meta_cc, taylor_spring_sc)

##need to import taylor spring cc samples -- these came from minnow trap
##Stomach content groups:
##Empty = empty
##Not taken = unknown
##Undetermined = unknown
##Algae = Algae
##Detritus = Detritus
##Insects = Insect (also: beetles, worms,)
##Fish = Fish
##Algae + insect = Omnivore (also: plant matter + insect)
##not including parasites as part of stomach contents
##brown goo? 

stomach_content_types <- as.data.frame(unique(fish_meta_cc$Stomach_Contents))

fish_meta_cc <- fish_meta_cc %>%
  mutate(stomach_content_group = case_when(
    startsWith(Stomach_Contents, "algae, insect") ~ "Omnivore",
    startsWith(Stomach_Contents, "insects, algae") ~ "Omnivore",
    startsWith(Stomach_Contents, "algae/bug") ~ "Omnivore",
    startsWith(Stomach_Contents, "plant matter/insect") ~ "Omnivore",
    startsWith(Stomach_Contents, "gastropod/algae") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus/insects") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus, insects") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus, fish, insects") ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus/bug") ~ "Omnivore",
    startsWith(Stomach_Contents, "worms") ~ "Insect",
    grepl("empty", Stomach_Contents) ~ "Empty",
    grepl("not taken", Stomach_Contents) ~ "Unknown",
    grepl("unid", Stomach_Contents) ~ "Unknown",
    grepl("Unid", Stomach_Contents) ~ "Unknown",
    grepl("worm/algae", Stomach_Contents) ~ "Omnivore",
    startsWith(Stomach_Contents, "detritus") ~ "Detritus",
    startsWith(Stomach_Contents, "insect") ~ "Insect",
    startsWith(Stomach_Contents, "insects, algae") ~ "Omnivore",
    grepl("beetle", Stomach_Contents) ~ "Insect", 
    grepl("insects/worm", Stomach_Contents) ~ "Insect",
    startsWith(Stomach_Contents, "algae") ~ "Algae",
    grepl("UNK", Stomach_Contents) ~"Unknown",
    grepl("damselfly", Stomach_Contents) ~"Insect",
    grepl("chironomid", Stomach_Contents) ~"Insect",
    grepl("brown goo", Stomach_Contents) ~ "Unknown"
  ))


fish_meta_cc_test <- fish_meta_cc %>%
  group_by(Site, Season, stomach_content_group) %>%
  summarise(count = n())
  
fish_meta_cc_test$Season <- ordered(fish_meta_cc_test$Season, 
                            levels = c("Spring", "Summer", "Fall"))

stomach_content_plot <- ggplot(fish_meta_cc_test, aes(x = Season, y = count, fill = stomach_content_group)) +
  geom_col(position = "stack") +
  facet_wrap(~Site) +
  theme_classic()
stomach_content_plot

##Calculate proportion of fish w/ stomach content groups
stomach_content_prop_function <- function(x){
  df_prep <- x
  n <- nrow(df_prep)
  df_prep_count <- df_prep %>%
    group_by(stomach_content_group) %>%
    summarise(count = n())
  prop <- df_prep_count$count/n
  prop <- as.data.frame(cbind(df_prep_count, prop)) 
  prop <- merge(df_prep, prop)
}

stomach_prop <- split(fish_meta_cc, paste0(fish_meta_cc$Season, fish_meta_cc$Site)) %>%
  map(stomach_content_prop_function) %>%
  bind_rows()

stomach_prop <- stomach_prop %>%
  select(Season, Site, stomach_content_group, count, prop) %>%
  unique()

stomach_prop$Season <- ordered(stomach_prop$Season, 
                               levels = c("Spring", "Summer", "Fall"))

stomach_prop <- stomach_prop %>%
  mutate(Site_Code = case_when(
    startsWith(Site, "H") ~ "Low Impact",
    startsWith(Site, "T") ~ "High Impact", 
    startsWith(Site, "P") ~ "Mid Impact"
  ) )

stomach_prop$Site_Code <- ordered(stomach_prop$Site_Code, 
                               levels = c("Low Impact", "Mid Impact", "High Impact"))

stomach_content_prop_plot <- ggplot(stomach_prop, aes(x = Season, y = prop, fill = stomach_content_group)) +
  geom_col(position = "stack", colour = "black")  +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2")) +
  facet_wrap(~Site_Code)+
  ylab("Proportion of Stomachs") +
  labs(fill = "Stomach Content") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(),  text = element_text(family = "Times New Roman"), strip.background = element_blank()) 

stomach_content_prop_plot

##plotting out % omnivory separately
stomach_prop_omni <- stomach_prop %>%
  filter(stomach_content_group == "Omnivore") %>%
  select(Site_Code, Season, stomach_content_group, prop) %>%
  complete(Site_Code, Season, stomach_content_group, fill = list(prop=0)) %>%
  group_by(Site_Code, stomach_content_group) %>%
  summarise_at(vars(prop), list(prop_mean = mean))


omni_prop_plot <- ggplot(stomach_prop_omni, aes(x = Site_Code, y = prop_mean, fill = stomach_content_group )) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#F0E442")) +
  ylab("Mean Proportion of Stomachs") +
  labs(fill = "Stomach Content") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(),  text = element_text(family = "Times New Roman"), strip.background = element_blank()) 

omni_prop_plot

##plotting out % insects separately
stomach_prop_insect <- stomach_prop %>%
  filter(stomach_content_group == "Insect") %>%
  select(Site_Code, Season, stomach_content_group, prop) %>%
  complete(Site_Code, Season, stomach_content_group, fill = list(prop=0)) %>%
  group_by(Site_Code, stomach_content_group) %>%
  summarise_at(vars(prop), list(prop_mean = mean))

##calculate magnitude of difference between low and high impact proportion of predators
mag_diff_pred <- 0.6500000/0.2623188

insect_prop_plot <- ggplot(stomach_prop_insect, aes(x = Site_Code, y = prop_mean, fill = stomach_content_group)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#009E73")) +
  ylab("Mean Proportion of Stomachs") +
  labs(fill = "Stomach Content") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(),  text = element_text(family = "Times New Roman"), strip.background = element_blank()) 

insect_prop_plot

stomach_content_fig <- ggarrange(stomach_content_prop_plot, 
                                 ggarrange(omni_prop_plot, insect_prop_plot, nrow = 2, labels = c("b)", "c)"), font.label = list(colour = "black", size = 12, family = "Times New Roman")),
                           ncol = 2,  labels = "a)",font.label = list(colour = "black", size = 12, family = "Times New Roman"))
stomach_content_fig


