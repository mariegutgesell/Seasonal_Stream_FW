###M3 Forest & Riparian Buffer Figure
##Date: Nov 3, 2022
##By: Marie Gutgesell

library(dplyr)
library(ggplot2)

setwd("~/Desktop/Seasonal Stream Food Webs/Data/")

landuse <- read.csv("Land-use Data/landuse_250_m3.csv")
buffer <- read.csv("Land-use Data/M3_2018_Buffer Width.csv") 

#%>%
 # filter(Measurement_Distance_Upstream == "0m") %>%
#  select(sitecode, Average_Buffer_Width_m)

buffer <- buffer %>%
  group_by(sitecode) %>%
  summarise(Avg_buff = mean(Average_Buffer_Width_m))

forest_buff <- left_join(landuse, buffer, by = "sitecode")

forest_buff <- forest_buff %>%
  mutate(Site = case_when(
    startsWith(sitecode, "HC") ~ "Low Impact",
    startsWith(sitecode, "AT") ~ "High Impact",
    startsWith(sitecode, "P4") ~ "Mid Impact",
  ))

forest_buff$Site <- ordered(forest_buff$Site,
                                      levels = c("Low Impact", "Mid Impact", "High Impact"))

forest_plot <- ggplot(forest_buff, aes(x = Site, y = forest, fill = Site)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("% Forest Cover_250m radius")
forest_plot

buffer_plot <- ggplot(forest_buff, aes(x = Site, y = Avg_buff, fill = Site)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),axis.text.y = element_text(size = 14),axis.title.y=element_text(size = 14), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Average Riparian Buffer Width (m)")
buffer_plot

forest_buff_plot <- ggarrange(forest_plot, buffer_plot, legend = "none",
                             labels = c("a)", "b)"),
                             ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
forest_buff_plot


####BELOW -- trying to figure out how to have two different y-axis, not working, come back to
#Scaling factor
sf <- max(forest_buff$forest)/max(forest_buff$Average_Buffer_Width_m)
#Transform
DF_long <- forest_buff %>%
  mutate(Average_Buffer_Width_m = Average_Buffer_Width_m*sf) %>%
  pivot_longer(names_to = "y_new", values_to = "val", forest:Average_Buffer_Width_m)

#Plot
ggplot(DF_long, aes(x=sitecode)) +
  geom_bar( aes(y = val, fill = y_new, group = y_new),
            stat="identity", position=position_dodge(),
            color="black", alpha=.6)  +
  scale_fill_manual(values = c("blue", "red")) +
  scale_y_continuous(name = "V2",labels = scales::comma,sec.axis = sec_axis(~./sf, name="V3",
                                                                            labels = scales::comma))+
  labs(fill='variable')+
  theme_bw()+
  theme(legend.position = 'top',
        plot.title = element_text(color='black',face='bold',hjust=0.5),
        axis.text = element_text(color='black',face='bold'),
        axis.title = element_text(color='black',face='bold'),
        legend.text = element_text(color='black',face='bold'),
        legend.title = element_text(color='black',face='bold'))+
  ggtitle('My barplot')


# Value used to transform the data
coeff <- 10

# A few constants
temperatureColor <- "#69b3a2"
priceColor <- rgb(0.2, 0.6, 0.9, 1)

ggplot(forest_buff, aes(x=sitecode)) +
  geom_bar(aes(y=forest), stat = "identity", colour = "black") + 
  geom_bar(aes(y=Average_Buffer_Width_m), stat = "identity", colour = "black") +
  scale_y_continuous(
    # Features of the first axis
    name = "% Forest",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Average Buffer Width (m)")
  ) + 
  theme_classic() 
