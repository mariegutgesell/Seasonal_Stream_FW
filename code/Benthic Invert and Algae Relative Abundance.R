##Seasonal Stream Food Webs -- Benthic invert & aquatic vegetation data (within stream resources)
##By: Marie Gutgesell
##Date: June 9, 2022
##Updated: Oct 24, 2023

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

##Function to remove outliers
##removing outliers
#' Replace outliers
#'
#' Replace outliers with NA. Outliers are defined as values that fall outside plus or minus
#' 1.5 * IQR.
#'
#' @return Numeric vector of same length.
#' @param x Numeric vector to replace.
#' @param na.rm Should NA values be replaced beforehand?
#'
#' @export
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  val <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  y
}

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

##do sqrt transformation to improve normality -- doesn't fix it perfectly but does improve it
benthic_total_abundance <- benthic_total_abundance %>%
  mutate(Abundance_sqrt = sqrt(x))

benthic_total_abundance_sub <- benthic_total_abundance %>% 
  group_by(Site, Season) %>% 
  mutate(Abundance_sqrt = remove_outliers(Abundance_sqrt)) %>% 
  ungroup() %>% 
  filter(!is.na(Abundance_sqrt))


##plots
##no outliers removed, not square root transformed
seasonal_benthic_plot <- ggplot(benthic_total_abundance, aes(x = Season, y = x, fill = Site)) +
  stat_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Abundance") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 
  
seasonal_benthic_plot

##no outliers removed, square root transformed
seasonal_benthic_plot_sqrt <- ggplot(benthic_total_abundance, aes(x = Season, y = Abundance_sqrt, fill = Site)) +
  stat_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Abundance") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 

seasonal_benthic_plot_sqrt

##outliers removed, square root transformed
seasonal_benthic_plot_sqrt_noout <- ggplot(benthic_total_abundance_sub, aes(x = Season, y = Abundance_sqrt, fill = Site)) +
  stat_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Abundance (Sqrt)") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank()) 

seasonal_benthic_plot_sqrt_noout

##i don't think removing outliers makes sense. there is such a low sample size. we are also not doing any statistics, so does it make any sense to sqrt transform? 
##The plot looks better, but that can not be reason enough, i don't think it makes sense.. in general this analysis is pretty sketchy

##to-do:  finish IDing last 4 samples - July 30th and Oct 15 for AT


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


##sqrt transformation and remove outliers
aqveg_data_subset <- aqveg_data_subset %>%
  mutate(aqveg_total_sqrt = sqrt(aqveg_total)) 

aqveg_data_subset_noout <- aqveg_data_subset %>% 
  group_by(Site, Season) %>% 
  mutate(aqveg_total_sqrt = remove_outliers(aqveg_total_sqrt)) %>% 
  ungroup() %>% 
  filter(!is.na(aqveg_total_sqrt))

##aquatic veg plot, not transformed and no outliers removed
seasonal_aqveg_plot <- ggplot(aqveg_data_subset, aes(x = Season, y = aqveg_total, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Relative Abundance Score") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())
  
seasonal_aqveg_plot

##aquatic veg plot, transformed and no outliers removed
seasonal_aqveg_plot_sqrt <- ggplot(aqveg_data_subset, aes(x = Season, y = aqveg_total_sqrt, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Relative Abundance Score") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

seasonal_aqveg_plot_sqrt

##square root transformed and outliers removed
seasonal_aqveg_plot_sqrt_noout <- ggplot(aqveg_data_subset_noout, aes(x = Season, y = aqveg_total_sqrt, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Relative Abundance Score (sqrt)") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

seasonal_aqveg_plot_sqrt_noout

##there are no outliers in the algae data -- so makes no difference, not sure square rooting makes sense 
###import fish density plot

source("code/Seasonal FW_Fish Density and Size.R")

##Final plot 
##Create figure legend for 3 sites 
data <- data.frame(
  Xdata = rnorm(3),
  Ydata = rnorm(3),
  LegendData = c("Low Impact", "Mid Impact", "High Impact")
)
data$LegendData <- factor(data$LegendData, levels = c("Low Impact", "Mid Impact", "High Impact"))
gplot <- ggplot(data, aes(Xdata, Ydata, color = LegendData)) +
  geom_point() +
  scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic()+
  guides(colour = guide_legend(title = "Site")) +
  theme(legend.title = element_text(family = "Times New Roman"), legend.text = element_text(size = 10, family = "Times New Roman"), legend.position = "right")
gplot

leg_fig <- get_legend(gplot)


fish_benthic_aqveg <- ggarrange(fish_density_all_plot, seasonal_benthic_plot, seasonal_aqveg_plot, legend = "right", common.legend = TRUE, legend.grob = leg_fig,
                                labels = c("a)", "b)", "c)"),
                                ncol = 1, nrow = 3, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
fish_benthic_aqveg

benthic_aqveg <- ggarrange(seasonal_benthic_plot, seasonal_aqveg_plot, legend = "right", common.legend = TRUE, legend.grob = leg_fig,
                                labels = c("a)", "b)"),
                                ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
benthic_aqveg

##Testing outliers and assumptions of anova -----------
#boxplot(x ~ Season + Site, benthic_total_abundance)

#hist(benthic_total_abundance$x)
#benthic_total_abundance[benthic_total_abundance<=0] <- 0.001


#benthic_total_abundance <- benthic_total_abundance %>%
#  mutate(Abundance_sqrt = sqrt(x)) %>%
#  mutate(Abundance_log = log(x))

#hist(benthic_total_abundance$Abundance_sqrt)##sqrt does improve normality but still not fully normally distributed
#hist(benthic_total_abundance$Abundance_log)
#qqnorm(benthic_total_abundance$x)
#qqnorm(benthic_total_abundance$Abundance_sqrt) ##i think sqrt transformation does better job of normalizing data
#qqnorm(benthic_total_abundance$Abundance_log)

#boxplot(Abundance_sqrt ~ Season + Site, benthic_total_abundance)
 

#boxplot(aqveg_total ~ Season +Site, aqveg_data_subset) ##no outliers
#hist(aqveg_data_subset$aqveg_total) 
#qqnorm(aqveg_data_subset$aqveg_total)

#aqveg_data_subset[aqveg_data_subset<=0] <- 0.001

#aqveg_data_subset <- aqveg_data_subset %>%
#  mutate(aqveg_total_sqrt = sqrt(aqveg_total)) %>%
#  mutate(aqveg_total_log = log(aqveg_total))

#hist(aqveg_data_subset$aqveg_total_sqrt) ##definitely improves normality 
#qqnorm(aqveg_data_subset$aqveg_total_sqrt)
#hist(aqveg_data_subset$aqveg_total_log) ##does not improve, 0's get separated
#qqnorm(aqveg_data_subset$aqveg_total_log)

##sqrt transformation helps improve normality most

###ANOVAs ---------------

benthic_anova <- aov(Abundance_sqrt ~ Site*Season, data = benthic_total_abundance_sub)
summary(benthic_anova)

benthic_tukey <- TukeyHSD(benthic_anova)
benthic_tukey_table <- benthic_tukey[["Site:Season"]]
write.csv(benthic_tukey_table, "outputs/benthic_inverts_tukey_table.csv")

aqveg_anova <- aov(aqveg_total_sqrt ~ Site*Season, data = aqveg_data_subset)
summary(aqveg_anova)

aqveg_tukey <- TukeyHSD(aqveg_anova)
aqveg_tukey_table <- aqveg_tukey[["Site:Season"]]
write.csv(aqveg_tukey_table, "outputs/aqveg_tukey_table.csv")


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

#mean_abundance_long <- gather(mean_abundance, group, abundance, Benthic_Inverts:Fish)

#mean_abundance_long$group <- ordered(mean_abundance_long$group,
#                                  levels = c("Fish", "Benthic_Inverts", "Aquatic_Vegetation"))

##not sure what this following series of plots is for, not using in MS
#mean_abundance_plot <- ggplot(mean_abundance_long, aes(y = abundance, x = group, fill = Site)) +
#  geom_col() +
#  facet_wrap(~Site) +
#  theme_classic() +
#  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
#  ylab("Mean Abundance")

#mean_abundance_plot

#mean_benthic_plot <- ggplot(benthic_abundance_mean, aes(x = Site, y = Abundance, fill = Site)) +
#  geom_col() +
#  theme_classic() +
#  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
#  ylab("Mean Abundance")
#mean_benthic_plot

#mean_aqveg_plot <- ggplot(aqveg_abundance_mean, aes(x = Site, y = Abundance, fill = Site)) +
#  geom_col() +
#  theme_classic() +
#  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
#  ylab("Mean Abundance Score")
#mean_aqveg_plot

#mean_fish_plot <- ggplot(fish_abundance_mean, aes(x = Site, y = name, fill = Site)) +
#  geom_col() +
#  theme_classic() +
#  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
#  ylab("Mean Abundance")
#mean_fish_plot


#fish_benthic_aqveg_mean <- ggarrange(mean_fish_plot, mean_benthic_plot, mean_aqveg_plot, legend = "none",
#                                labels = c("a)", "b)", "c)"),
#                                ncol = 1, nrow = 3, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
#fish_benthic_aqveg_mean

##Calculating R:C ratio mean
seasonal_r_c_ratio <- seasonal_mean_abundance %>%
  mutate(R_C = Aquatic_Vegetation/Benthic_Inverts)

r_c_season_plot <- ggplot(seasonal_r_c_ratio, aes(x = Season, y = R_C, fill = Site)) +
  geom_col(colour = "black") +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Mean Seasonal R:C Ratio") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

r_c_season_plot


biomass_pyramid_seasonal <- ggarrange(fish_density_all_plot, r_c_season_plot, legend = "none",
                             labels = c("a)", "b)"),
                             ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
biomass_pyramid_seasonal

##so i think it makes the most sense to use not transformed and no outlier removed, as not doing statistics so transformation does not make sense... and not enough samples to make outliers legit.. 


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

##Calculating P:C ratio - mean seasonal

seasonal_p_c_ratio <- seasonal_mean_abundance %>%
  mutate(P_C = Fish/Aquatic_Vegetation)

p_c_season_plot <- ggplot(seasonal_p_c_ratio, aes(x = Season, y = P_C, fill = Site)) +
  geom_col(colour = "black") +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Mean P:C Ratio") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

p_c_season_plot

annual_p_c_ratio <- annual_mean_abundance %>%
  mutate(P_C = Fish/Aquatic_Vegetation)

annual_p_c_plot <- ggplot(annual_p_c_ratio, aes(x = Site, y = P_C, fill = Site)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Mean P:C Ratio")
annual_p_c_plot

biomass_pyramid_seasonal <- ggarrange(fish_density_all_plot, r_c_season_plot, legend = "none",
                                      labels = c("a)", "b)"),
                                      ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
biomass_pyramid_seasonal

biomass_pyramid_seasonal_2 <- ggarrange(p_c_season_plot, r_c_season_plot, legend = "none",
                                      labels = c("a)", "b)"),
                                      ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
biomass_pyramid_seasonal_2


mean_fish_plot <- ggplot(annual_mean_abundance, aes(x = Site, y = Fish, fill = Site)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 12),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Mean Abundance")
mean_fish_plot

biomass_pyramid <- ggarrange(mean_fish_plot, annual_r_c_plot, legend = "none",
                             labels = c("a)", "b)"),
                             ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
biomass_pyramid

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

###NOT USING ------------------------
#####R:C ratio -- seasonal and annual -- using transformed and outlier removed data, not sure which is correct? first want to see how it is different

##Putting together mean densities --- SEASONAL
seasonal_benthic_abundance_mean <- benthic_total_abundance_sub %>%
  group_by(Site, Season) %>%
  summarise_at(vars(Abundance_sqrt), list(Abundance = mean))

aqveg_data_subset <- na.omit(aqveg_data_subset)
seasonal_aqveg_abundance_mean <- aqveg_data_subset %>%
  group_by(Site, Season) %>%
  summarise_at(vars(aqveg_total_sqrt), list(Abundance = mean))

names(fish_density_all)[4] <- "Site"
names(fish_density_all)[2] <- "Site_Name"

seasonal_mean_abundance <- left_join(seasonal_benthic_abundance_mean, seasonal_aqveg_abundance_mean, by = c("Site", "Season")) %>%
  left_join(fish_density_all[,c(1,3:4)], by = c("Site", "Season"))
names(seasonal_mean_abundance)[3] <- "Benthic_Inverts"
names(seasonal_mean_abundance)[4] <- "Aquatic_Vegetation"
names(seasonal_mean_abundance)[5] <- "Fish"

##Calculating R:C ratio seasonal mean
seasonal_r_c_ratio <- seasonal_mean_abundance %>%
  mutate(R_C = Aquatic_Vegetation/Benthic_Inverts)

3.7/0.16556

r_c_season_plot <- ggplot(seasonal_r_c_ratio, aes(x = Season, y = R_C, fill = Site)) +
  geom_col(colour = "black") +
  facet_wrap(~Site) +
  theme_classic() +
  ylab("Mean Seasonal R:C Ratio") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),axis.text.y = element_text(size = 14),axis.title.y=element_text(size = 14), axis.title.x = element_blank(), legend.position = "bottom", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())

r_c_season_plot

fish_density_all_plot

biomass_pyramid_seasonal <- ggarrange(fish_density_all_plot, r_c_season_plot, legend = "none",
                                      labels = c("a)", "b)"),
                                      ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
biomass_pyramid_seasonal

##Mean R:C ratio  -- ANNUAL
annual_benthic_abundance_mean <- benthic_total_abundance_sub %>%
  group_by(Site) %>%
  summarise_at(vars(Abundance_sqrt), list(Abundance = mean))

aqveg_data_subset <- na.omit(aqveg_data_subset)
annual_aqveg_abundance_mean <- aqveg_data_subset %>%
  group_by(Site) %>%
  summarise_at(vars(aqveg_total_sqrt), list(Abundance = mean))

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Mean Annual R:C Ratio")
annual_r_c_plot

mean_fish_plot <- ggplot(annual_mean_abundance, aes(x = Site, y = Fish, fill = Site)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Mean Annual Abundance")
mean_fish_plot

biomass_pyramid <- ggarrange(mean_fish_plot, annual_r_c_plot, legend = "none",
                             labels = c("a)", "b)"),
                             ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
biomass_pyramid

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
  

