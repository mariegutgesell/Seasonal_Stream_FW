##Isotopic Niche Space Figure
source("~/Desktop/Seasonal Stream Food Webs/Code/Stable Isotope Analysis/Seasonal FW_Isotopic Niche Width_Trophic Metrics_NicheRover.R")

##join ellipse.plot and SEA_plot

ellipse.plot
SEA_plot
##Final plot 
##Create figure legend for 3 sites 
data <- data.frame(
  Xdata = rnorm(3),
  Ydata = rnorm(3),
  LegendData = c("Low Impact", "Mid Impact", "High Impact"),
  LegendShape = c("Spring", "Summer", "Fall")
)

data$LegendData <- factor(data$LegendData, levels = c("Low Impact", "Mid Impact", "High Impact"))
data$LegendShape <- factor(data$LegendShape, levels = c("Spring", "Summer", "Fall"))
gplot <- ggplot(data, aes(Xdata, Ydata, color = LegendData, shape = LegendShape)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  scale_shape_manual(values = c( 1, 2, 0)) +
  theme_classic()+
  guides( shape = guide_legend(title = "Season"), colour = guide_legend(title = "Site")) +
  theme(legend.title = element_text(family = "Times New Roman", size = 12), legend.text = element_text(size = 12, family = "Times New Roman"), legend.position = "right")
gplot

leg_fig <- get_legend(gplot)


iso_niche_fig <- ggarrange(ellipse.plot, SEA_plot, legend = "right", common.legend = TRUE, legend.grob = leg_fig,
                                labels = c("a)", "b)"),
                                ncol = 2, nrow = 1, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
iso_niche_fig

iso_niche_fig_2 <- ggarrange(ellipse.plot, 
          ggarrange(SEA_plot, SEA_plot_mean, legend = "none", nrow = 2, labels = c("b)", "c)"), font.label = list(colour = "black", size = 12, family = "Times New Roman")),legend = "right", common.legend = TRUE, legend.grob = leg_fig,
          ncol = 2,  labels = "a)",font.label = list(colour = "black", size = 12, family = "Times New Roman"))

iso_niche_fig_2


