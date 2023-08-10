###Nutrient Data of M3 sites
library(ggpubr)
setwd("~/Desktop/Seasonal Stream Food Webs/Data/")

n_p <- read.csv("Nutrient Data/2018_M3_TN_TP.csv")


n_p <- n_p %>%
  mutate(Site_Type = case_when(
    startsWith(Site, "H") ~ "Low Impact",
    startsWith(Site, "T") ~ "High Impact",
    startsWith(Site, "P") ~ "Mid Impact",
  ))

n_p$Site_Type <- ordered(n_p$Site_Type,
                            levels = c("Low Impact", "Mid Impact", "High Impact"))


n_plot <- ggplot(n_p, aes(x = Site_Type, y = TN_mg_per_L, fill = Site_Type)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),axis.text.y = element_text(size = 14),axis.title.y=element_text(size = 14), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Total Nitrogen (mg/L)")
n_plot

p_plot <- ggplot(n_p, aes(x = Site_Type, y = TP_mg_per_L, fill = Site_Type)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 10), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Times New Roman"), strip.background = element_blank(),strip.text = element_blank())+
  ylab("Total Phosphorous (mg/L)")
p_plot

n_p_plot <- ggarrange(n_plot, p_plot, legend = "none",
                              labels = c("a)", "b)"),
                              ncol = 1, nrow = 2, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
n_p_plot


##ANOVA
##testing anova assumptions
##normality
hist(n_p$TN_mg_per_L)
qqPlot(n_p$TN_mg_per_L)

hist(n_p$TP_mg_per_L)
qqPlot(n_p$TP_mg_per_L)


n_p <- n_p%>%
  mutate(TN_log = log(TN_mg_per_L)) %>%
  mutate(TP_log = log(TP_mg_per_L)) %>%
  mutate(TN_sqrt = sqrt(TN_mg_per_L)) %>%
  mutate(TN_sq = '^'(TN_mg_per_L, 2)) %>%
  mutate(TN_exp = exp(TN_mg_per_L))
hist(n_p$TN_exp)
qqPlot(n_p$TN_log)


n_anova <- aov(TN_mg_per_L ~ Site_Type, data = n_p)
summary(n_anova)
TukeyHSD(n_anova)

p_anova <- aov(TP_mg_per_L ~ Site_Type, data = n_p)
summary(p_anova)
TukeyHSD(p_anova)
