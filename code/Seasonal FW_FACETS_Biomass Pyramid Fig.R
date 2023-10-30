##FACETS_Biomass Pyramid Code -- Final Figure 4 
##Using untransformed and no outlier removed data, I think this makes most sense as not doing any statistics on it .. so don't have a reason to transform, and not enough samples to remove outliers

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

##Read in datafiles
benthic_data <- read.csv ("data/M3 Insect Data/Benthic Monitoring ID Sheet_MASTER.csv")
benthic_data$newdate <-strptime(as.character(benthic_data$Date), "%Y-%m-%d")
#benthic_data$newdate <- str_trim(benthic_data$newdate)

benthic_data$newdate
#benthic_data$txtdate <- format(benthic_data$newdate, "%y-%m-%d")
benthic_data$newdate <- as.Date(benthic_data$newdate)



##Import csv file containing cleaned stream conditions data  
data <- read.csv("data/Stream Condition Data/2x-weekly Stream Conditions_cleaned_MASTERCSV.csv", na.strings=c("","NA"))

data$newdate <-strptime(as.character(data$Date), "%d/%m/%Y")
data$newdate <- as.Date(data$newdate)

#data$txtdate <- format(data$newdate, "%y-%m-%d")

##Calculate total algae and macrophyte availability
algae <- c("Algae_floating", "Algae_filaments", "Algae_attached", "Algae_slim_crusts")
data$algae_total <- rowSums(data[,algae])

macrophyte <- c("Macrophyte_Emergent", "Macrophyte_rooted_floating", "Macrophyte_submerged", "Macrophyte_free_floating")
data$macrophyte_total <- rowSums(data[,macrophyte])


##Benthic Data --------
##Subsampling date across range of 5 sampling points before e-fishing took place -- NOTE: 4 samples still to be ID: July 30, all 3 sites and Oct 15, 2018 for Taylor
data_spring <- benthic_data[benthic_data$newdate >= "2018-05-14" & benthic_data$newdate <= "2018-05-28",]
data_summer <- benthic_data[benthic_data$newdate >= "2018-07-16" & benthic_data$newdate <= "2018-07-30",]
data_fall <- benthic_data[benthic_data$newdate >= "2018-10-08" & benthic_data$newdate <= "2018-10-22",]

benthic_data_subset <- rbind(data_spring, data_summer, data_fall)
rm(data_spring, data_summer, data_fall)
benthic_data_subset$abundance <- as.numeric(benthic_data_subset$abundance)

##Calculate total abundance of invertebrates
benthic_total_abundance <- aggregate(benthic_data_subset$abundance, by = list(Site_ID = benthic_data_subset$Site_ID, newdate = benthic_data_subset$newdate), FUN = sum)
benthic_total_abundance$Site_ID <- ordered(benthic_total_abundance$Site_ID,
                                           levels = c("HC", "EP4", "AT"))


##Create Boxplot of Benthic Abundances
benthic_total_abundance$newdate <- as.character(benthic_total_abundance$newdate)
benthic_total_abundance <- benthic_total_abundance %>%
  mutate(Season = case_when(
    startsWith(newdate, "2018-05") ~ "Spring",
    startsWith(newdate, "2018-07") ~ "Summer", 
    startsWith(newdate, "2018-10") ~ "Fall"
  ) )

benthic_total_abundance$Site_ID <- as.character(benthic_total_abundance$Site_ID)
benthic_total_abundance <- benthic_total_abundance %>%
  mutate(Site = case_when(
    startsWith(Site_ID, "HC") ~ "Low Impact",
    startsWith(Site_ID, "EP") ~ "Mid Impact", 
    startsWith(Site_ID, "AT") ~ "High Impact"
  ) )


benthic_total_abundance$Season <- ordered(benthic_total_abundance$Season,
                                          levels = c("Spring", "Summer", "Fall"))
benthic_total_abundance$Site <- ordered(benthic_total_abundance$Site,
                                        levels = c("Low Impact", "Mid Impact", "High Impact"))

##Seasonal Benthic Invert Plot
##no outliers removed, not square root transformed
seasonal_benthic_plot <- ggplot(benthic_total_abundance, aes(x = Season, y = x, fill = Site)) +
  stat_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Abundance") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 

seasonal_benthic_plot


##Aquatic vegetation --------
##Subsampling date across range of 5 sampling points before e-fishing took place --
data_spring <- data[data$newdate >= "2018-05-14" & data$newdate <= "2018-05-28",]
data_summer <- data[data$newdate >= "2018-07-16" & data$newdate <= "2018-07-30",]
data_fall <- data[data$newdate >= "2018-10-08" & data$newdate <= "2018-10-22",]

aqveg_data_subset <- rbind(data_spring, data_summer, data_fall)
#rm(data_spring, data_summer, data_fall)

aqveg_data_subset <- aqveg_data_subset %>%
  select(Site_Code, newdate, algae_total, macrophyte_total) %>%
  mutate(aqveg_total = algae_total + macrophyte_total)

#aqveg_total_abundance <- aggregate(aqveg_data_subset$abundance, by = list(Site_ID = aqveg_data_subset$Site_ID, newdate = aqveg_data_subset$newdate), FUN = sum)
aqveg_data_subset$Site_Code <- ordered(aqveg_data_subset$Site_Code,
                                       levels = c("HC", "EP4", "AT"))

##still missing data for: EP4 and HC - 19-07-2018 -- not in datasheet, need to find raw sheet and input info

aqveg_data_subset$newdate <- as.character(aqveg_data_subset$newdate)
aqveg_data_subset <- aqveg_data_subset %>%
  mutate(Season = case_when(
    startsWith(newdate, "2018-05") ~ "Spring",
    startsWith(newdate, "2018-07") ~ "Summer", 
    startsWith(newdate, "2018-10") ~ "Fall"
  ) )

aqveg_data_subset$Site_Code <- as.character(aqveg_data_subset$Site_Code)
aqveg_data_subset <- aqveg_data_subset %>%
  mutate(Site = case_when(
    startsWith(Site_Code, "HC") ~ "Low Impact",
    startsWith(Site_Code, "EP") ~ "Mid Impact", 
    startsWith(Site_Code, "AT") ~ "High Impact"
  ) )

aqveg_data_subset$Season <- ordered(aqveg_data_subset$Season,
                                    levels = c("Spring", "Summer", "Fall"))

aqveg_data_subset$Site <- ordered(aqveg_data_subset$Site,
                                  levels = c("Low Impact", "Mid Impact", "High Impact"))

##aquatic veg plot, not transformed and no outliers removed
seasonal_aqveg_plot <- ggplot(aqveg_data_subset, aes(x = Season, y = aqveg_total, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Relative Abundance Score") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

seasonal_aqveg_plot

##Figure S4: Combining Aq_Veg and Benthic Invert Seasonal Abundance Plots
##Create figure legend for 3 sites  -- come back to this, likely want to change legend to be same as in main fig 
data <- data.frame(
  Xdata = rnorm(3),
  Ydata = rnorm(3),
  LegendData = c("Low Impact", "Mid Impact", "High Impact")
)
data$LegendData <- factor(data$LegendData, levels = c("Low Impact", "Mid Impact", "High Impact"))
gplot <- ggplot(data, aes(Xdata, fill = LegendData)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  theme_classic()+
  guides(colour = guide_legend(title = "Site")) +
  theme(legend.title = element_text(family = "Times New Roman"), legend.text = element_text(size = 12, family = "Times New Roman"), legend.position = "bottom")+
  guides(fill = guide_legend(title="Site"))
gplot

leg_fig <- get_legend(gplot)

benthic_aqveg <- ggarrange(seasonal_benthic_plot, seasonal_aqveg_plot, legend = "bottom", common.legend = TRUE, legend.grob = leg_fig,
                           labels = c("a)", "b)"),
                           ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
benthic_aqveg

##Import fish density plot

source("code/Seasonal FW_Fish Density and Size.R")

##Putting together mean densities --- SEASONAL
seasonal_benthic_abundance_mean <- benthic_total_abundance %>%
  group_by(Site, Season) %>%
  summarise_at(vars(x), list(Abundance = mean))

aqveg_data_subset <- na.omit(aqveg_data_subset)
seasonal_aqveg_abundance_mean <- aqveg_data_subset %>%
  group_by(Site, Season) %>%
  summarise_at(vars(aqveg_total), list(Abundance = mean))

names(fish_density_all)[4] <- "Site"
names(fish_density_all)[2] <- "Site_Name"

seasonal_mean_abundance <- left_join(seasonal_benthic_abundance_mean, seasonal_aqveg_abundance_mean, by = c("Site", "Season")) %>%
  left_join(fish_density_all[,c(1,3:4)], by = c("Site", "Season"))
names(seasonal_mean_abundance)[3] <- "Benthic_Inverts"
names(seasonal_mean_abundance)[4] <- "Aquatic_Vegetation"
names(seasonal_mean_abundance)[5] <- "Fish"

##Calculating R:C ratio mean
seasonal_r_c_ratio <- seasonal_mean_abundance %>%
  mutate(R_C = Aquatic_Vegetation/Benthic_Inverts)

##Change in R:C ratio from lowest season to highest season
##maximum difference in R:C ratio from spring to fall at low-impact site
hc_max_diff <- 0.14285714/0.03149606

##maximum difference in R:C ratio from spring to fall at mid-impact site
ep_max_diff <- 0.14018692/0.03571429

##maximum difference in R:C ratio from spring to fall at high-impact site
at_max_diff <- 3.70000000/0.16556291

r_c_season_plot <- ggplot(seasonal_r_c_ratio, aes(x = Season, y = R_C, fill = Site)) +
  geom_col(colour = "black") +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Mean Seasonal R:C Ratio") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

r_c_season_plot

##Mean annual R:C ratio  
annual_benthic_abundance_mean <- benthic_total_abundance %>%
  group_by(Site) %>%
  summarise_at(vars(x), list(Abundance = mean))

aqveg_data_subset <- na.omit(aqveg_data_subset)
annual_aqveg_abundance_mean <- aqveg_data_subset %>%
  group_by(Site) %>%
  summarise_at(vars(aqveg_total), list(Abundance = mean))

annual_fish_abundance_mean <- fish_density_all%>%
  group_by(Site) %>%
  summarise_at(vars(Number_of_Individuals), list(name = mean))


annual_mean_abundance <- left_join(annual_benthic_abundance_mean, annual_aqveg_abundance_mean, by = "Site") %>%
  left_join(annual_fish_abundance_mean, by = "Site")

names(annual_mean_abundance)[2] <- "Benthic_Inverts"
names(annual_mean_abundance)[3] <- "Aquatic_Vegetation"
names(annual_mean_abundance)[4] <- "Fish"

annual_r_c_ratio <- annual_mean_abundance %>%
  mutate(R_C = Aquatic_Vegetation/Benthic_Inverts)

annual_r_c_plot <- ggplot(annual_r_c_ratio, aes(x = Site, y = R_C, fill = Site)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Mean R:C Ratio")
annual_r_c_plot

mean_fish_plot <- ggplot(annual_mean_abundance, aes(x = Site, y = Fish, fill = Site)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Mean Abundance")
mean_fish_plot



##Joining seasonal and annual mean R:C figs together
##create figure legend
library(ggpubr)
data <- data.frame(
  Xdata = rnorm(3),
  Ydata = rnorm(3),
  LegendData = c("Low Impact", "Mid Impact", "High Impact")
)
data$LegendData <- factor(data$LegendData, levels = c("Low Impact", "Mid Impact", "High Impact"))
gplot <- ggplot(data, aes(Xdata, fill = LegendData)) +
  geom_bar(colour = "black") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  theme_classic()+
  guides(colour = guide_legend(title = "Site")) +
  theme(legend.title = element_text(family = "Times New Roman"), legend.text = element_text(size = 12, family = "Times New Roman"), legend.position = "bottom")+
  guides(fill = guide_legend(title="Site"))
gplot

leg_fig <- get_legend(gplot)


biomass_pyramid_seasonal_annual <- ggarrange(fish_density_all_plot, mean_fish_plot, r_c_season_plot, annual_r_c_plot, legend = "bottom", common.legend = TRUE, legend.grob = leg_fig,
                                             labels = c("a)", "b)", "c)", "d)"),
                                             ncol = 2, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))


biomass_pyramid_seasonal_annual


