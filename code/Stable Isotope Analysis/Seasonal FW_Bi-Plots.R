##Code: Seasonal Stream Food Web Stable Isotope Bi-plots - exploring data across site
##By: Marie Gutgesell
##Date: May 4, 2022

source("code/Stable Isotope Analysis/Seasonal FW_SI Dataset Formation.R")

##Uncorrected dH biplots

##only CC, algae and det
SI_cc_al_dt$Season <- ordered(SI_cc_al_dt$Season, 
                              levels = c("Spring", "Summer", "Fall"))



seasonal_stream_si_plot_site <-  SI_cc_al_dt %>%
  ggplot(aes(x = d2H,y = d15N, color = Species)) +
  geom_point(size = 4) +
  theme_bw() +
  xlab("d2H") + 
  ylab("d15N") + 
  ggtitle("Seasonal Stream Food Web Data") +
  theme(text=element_text(size=10)) +
  theme(axis.line.x = element_line(color="black", size = 2),
        axis.line.y = element_line(color="black", size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Site + Season)

seasonal_stream_si_plot_site

##Creek chub and all baselines
###remove other fish species
SI_all <- SI_all %>%
  na.omit()

SI_all$Season <- ordered(SI_all$Season, 
                              levels = c("Spring", "Summer", "Fall"))

seasonal_stream_si_plot_site_all <-  SI_all %>%
  ggplot(aes(x = d2H,y = d15N, color = Species)) +
  geom_point(size = 4) +
  theme_bw() +
  xlab("d2H") + 
  ylab("d15N") + 
  ggtitle("Seasonal Stream Food Web Data") +
  theme(text=element_text(size=10)) +
  theme(axis.line.x = element_line(color="black", size = 2),
        axis.line.y = element_line(color="black", size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Site + Season)

seasonal_stream_si_plot_site_all

seasonal_stream_si_plot_site_all <-  SI_all %>%
  ggplot(aes(x = d2H,y = d15N, color = Trophic_Group)) +
  geom_point(size = 4) +
  theme_bw() +
  xlab("d2H") + 
  ylab("d15N") + 
  ggtitle("Seasonal Stream Food Web Data") +
  theme(text=element_text(size=10)) +
  theme(axis.line.x = element_line(color="black", size = 2),
        axis.line.y = element_line(color="black", size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Site + Season)

seasonal_stream_si_plot_site_all


##dH TRANSFORMATION  ---------------
##Equations for Aquatic Inverts
####dHadujustedendmember = [dHendmember - (w1 x dHwater)]/(1-w1)
##Equations For Creek Chub: 
#### w2 = 1 - (1-w1)^2 
####dHadjusted = [dHraw - (w2*dHwater)]/1-w2

##Equations from Page et al., 2017
##dH20 estimation taken as mean of dH20 values of two sites that span our study sites from Gibson et al., 2020 - dHwater = -62.7
##w1 estimates from Solomon et al., 2009 --> have tried both w = 0.12 (value for shredder insects) and w = 0.17 (avg. for all organisms) made no difference to regressions, in Solomon et al., 2009 there is high variation in invertebrate w, and within each order we have we have shredders and scrappers... not sure how much will change outcome, perhaps best to use average for all organisms (this is what Matt did for Em's MS)

w2 = 1 - (1-0.17)^2

creekchub <- SI_cc_al_dt %>%
  filter(Species == "Creek Chub")

output_dHadjust_cc <- matrix(nrow = 86, ncol = 1)
dH_adjust_cc <- for(i in 1:nrow(creekchub)){
  output_dHadjust_cc[i,] <- (creekchub$d2H[i] - (0.31 * (-62.7)))/(1-0.31)
}
output_dHadjust_cc
creekchub_adjusted <- cbind(creekchub, output_dHadjust_cc)

SI_cc_al_dt_adj <- left_join(SI_cc_al_dt, creekchub_adjusted, by = c("Sample_ID", "d13C", "d15N", "d2H", "CN_ratio", "per_C", "per_N" ,"Site", "Species", "Season", "Trophic_Group"))

seasonal_stream_si_plot_site_adj <- ggplot(SI_cc_al_dt_adj) +
  geom_point(aes(x = d2H, y = d15N, color = Species)) +
  geom_point(aes(x = output_dHadjust_cc, y = d15N, colour = "purple")) +
  xlab("d2H_adj") + 
  ylab("d15N") + 
  ggtitle("Seasonal Stream Food Web Data") +
  theme(text=element_text(size=10)) +
  theme(axis.line.x = element_line(color="black", size = 2),
        axis.line.y = element_line(color="black", size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Site + Season)

seasonal_stream_si_plot_site_adj

##Plotting only creek chub, all on one plot
cc_seasonal_stream_si_plot_site_adj <- ggplot(SI_cc_al_dt_adj) +
  geom_point(aes(x = output_dHadjust_cc, y = d15N, colour = Season, shape = Site)) +
  xlab("d2H_adj") + 
  ylab("d15N") + 
  ggtitle("Seasonal Stream Food Web Data") +
  theme(text=element_text(size=10)) +
  theme(axis.line.x = element_line(color="black", size = 2),
        axis.line.y = element_line(color="black", size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

cc_seasonal_stream_si_plot_site_adj
##still need to do dH adjustment of aquatic inverts 


rm(SI_cc_al_dt, creekchub, creekchub_adjusted, output_dHadjust_cc, seasonal_stream_si_plot_site, seasonal_stream_si_plot_site_adj, seasonal_stream_si_plot_site_all)
