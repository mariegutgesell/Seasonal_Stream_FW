##Seasonal Stream Food Web Isotope Analysis - Isotopic Niche Width using SIBER on trophic position and % coupling metrics

##

set.seed(1)
library(SIBER)


source("code/Stable Isotope Analysis/Seasonal FW_TP & Coupling.R")
rm(cc_seasonal_stream_si_plot_site_adj, coupling_boxplot, SI_all, SI_cc_al_dt_adj, TP_boxplot)
##Set up dataframe for SIBER ------------
seasonal_df_split_pred_tp_bound$Season <- as.character(seasonal_df_split_pred_tp_bound$Season)
siber_cc_tm <- seasonal_df_split_pred_tp_bound %>%
  mutate(group = case_when(
    startsWith(Season, "Spring") ~ "1",
    startsWith(Season, "Summer") ~ "2",
    startsWith(Season, "Fall") ~ "3",
  ))

siber_cc_tm$Site <- as.character(siber_cc_tm$Site)
siber_cc_tm <- siber_cc_tm %>%
  mutate(community = case_when(
    startsWith(Site, "Hawk") ~ "1",
    startsWith(Site, "Tayl") ~ "3",
    startsWith(Site, "Perl") ~ "2",
  ))

##Set up dataframe so formatted for SIBER
siber_cc_tm <- siber_cc_tm %>%
  select(TC_pred, TP_pred, group, community)

names(siber_cc_tm)[1] <- "iso1"  ##% terrestrial coupling
names(siber_cc_tm)[2] <- "iso2" ##trophic position
siber_cc_tm$group <- as.integer(siber_cc_tm$group)
siber_cc_tm$community <- as.integer(siber_cc_tm$community)
typeof(siber_cc_tm$iso1)

siber_cc_tm <- as.data.frame(siber_cc_tm) %>%
  arrange(community, group)
str(siber_cc_tm)

##Testing normal distribution
library(mvnormtest)
shapiro.test(siber_cc_tm$iso1)
shapiro.test(siber_cc_tm$iso2)


##Create siber object ----------
siber_cc_object <- createSiberObject(siber_cc_tm)

##Plotting SEA -------------
#Create lists of plotting arguments to be passed onwards to each of the three plotting functions
community.hulls.args <- list(col = 1, lty = 1, lwd = 1) 
group.ellipses.args <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2) ##p = NULL creates Maximum Liklihood Standard Ellipse (incorporates approx. 40% of data)
group.hulls.args <- list(lty = 2, col = "grey20") 

##Make plot -- this currently plots the maximum likelihood standard ellipse
par(mfrow = c(1,1))
plotSiberObject(siber_cc_object,
                ax.pad = 2,
                hulls = F, community.hulls.args = community.hulls.args,
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = F, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression(Terrestrial_Coupling),
                ylab = expression(Trophic_Position)
)


group.ML <- groupMetricsML(siber_cc_object)
print(group.ML)

legend("topright", colnames(group.ML), 
       pch = c(1,1,1,2,2,2,3,3,3), col = c(1:3, 1:3, 1:3), lty = 1)

#Make ellipse plot in ggplot
siber_cc_ggplot <- siber_cc_tm %>% mutate(group = factor(group), 
                                       community = factor(community),
                                       TC_pred = iso1, 
                                       TP_pred = iso2,
                                       .keep = "unused")

first.plot <- ggplot(data = siber_cc_ggplot, 
                     aes(x = TC_pred, 
                         y = TP_pred)) + 
  geom_point(aes(color = community, shape = group), size = 5) +
  ylab("Trophic Position") +
  xlab("% Terrestrial Coupling") + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) + 
  theme_classic()
first.plot

# decide how big an ellipse you want to draw
#p.ell <- pchisq(1,2) ##use this for a standard ellipse
p.ell <- pchisq(1,2)
ellipse.plot <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community), 
                   fill = community,
                   color = community), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman")) 

ellipse.plot





##Function to extract range of x and y of ellipses for each site and season
ellipse_range_function <- function(x){
  df_prep <- x
  group <- df_prep[1,3]
  community <- df_prep[1,4]
  m <- nrow(df_prep)
  mu <- colMeans(df_prep[,1:2])
  sigma <- cov(df_prep[,1:2])
  p <- pchisq(1,2)
  plot(df_prep[,2] ~ df_prep[,1])
  ellipse <- addEllipse(mu, sigma, m, p.interval = p, small.sample = TRUE)
  ellipse <- as.data.frame(ellipse)
  ellipse_ymax <- max(ellipse$V2)
  ellipse_ymin <- min(ellipse$V2)
  ellipse_xmax <- max(ellipse$V1)
  ellipse_xmin <- min(ellipse$V1)
  ellipse_yrange <- ellipse_ymax - ellipse_ymin
  ellipse_xrange <- ellipse_xmax - ellipse_xmin
  ellipse_metrics <- as.data.frame(cbind(ellipse_xrange, ellipse_yrange, group, community))
}


ellipse_range <- split(siber_cc_tm, paste0(siber_cc_tm$community, siber_cc_tm$group)) %>%
  map(ellipse_range_function) %>%
  bind_rows()
ellipse_range$community.group <- paste(ellipse_range$community, ".", ellipse_range$group)

xrange_plot <- ggplot(ellipse_range, aes(x = community.group, y = ellipse_xrange)) +
  geom_point() +
  theme_classic() +
  ylab("%Terrestrial Coupling Range")
xrange_plot

yrange_plot <- ggplot(ellipse_range, aes(x = community.group, y = ellipse_yrange)) +
  geom_point() +
  theme_classic() +
  ylab("Trophic Position Range")
yrange_plot

##Calculate group summary metrics -- TA == Total Area, SEA = Standard Ellipse Area, SEAc = corrected SEA ---------
group.ML <- groupMetricsML(siber_cc_object)
print(group.ML)

##plotting summary metrics
group_ML_transposed <- as.data.frame(t(group.ML))
group_ML_transposed$names <- rownames(group_ML_transposed)
group.ML_plot <- ggplot(group_ML_transposed, aes(x = names, y = SEAc)) +
  geom_point() +
  theme_classic()
group.ML_plot

###Fitting Bayesian models to the data ---------------
library(rjags) ##note this requires JAGS to be installed on computer (JAGS = Just Another Gibbs Sampler)
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber_cc_object, parms, priors)

##Comparing among groups using Bayesian Standard Ellipse Area -----------

SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

##Determining Ellipse overlap of maximum liklihood ellipses -------------
##note: the p = 0.95 is from the link below, how do i know if this is using the SEAc? is it b/c of the maxLikOverlap function? no, based on the r documentation, the p should be determined same as in addEllipse function, where if p = NULL that plots the SEA
##https://cran.r-project.org/web/packages/SIBER/vignettes/siber-comparing-populations.html
##When p = NULL, this uses area of SEA, how do i change it so it is area of SEAc? 
##HawkCliff
##overlap between the SEAc of groups 1.2 and 1.3
overlap.G1.2.G1.3 <- maxLikOverlap("1.2", "1.3", siber_cc_object, p = NULL, n =)
print(overlap.G1.2.G1.3)
##overlap between the SEAc of groups 1.1 and 1.3
overlap.G1.1.G1.3 <- maxLikOverlap("1.1", "1.3", siber_cc_object, p = NULL, n =)
print(overlap.G1.1.G1.3)
##overlap between the SEAc of groups 1.1 and 1.2
overlap.G1.1.G1.2 <- maxLikOverlap("1.1", "1.2", siber_cc_object, p = NULL, n =)
print(overlap.G1.1.G1.2)

##repeat this for EP4 and Taylor
##EP4
##overlap between the SEAc of groups 2.1 and 2.2
overlap.G2.1.G2.2 <- maxLikOverlap("2.1", "2.2", siber_cc_object, p = NULL, n =)
print(overlap.G2.1.G2.2)
##overlap between the SEAc of groups 2.2 and 2.3
overlap.G2.2.G2.3 <- maxLikOverlap("2.2", "2.3", siber_cc_object, p = NULL, n =)
print(overlap.G2.2.G2.3)
##overlap between the SEAc of groups 2.1 and 2.3
overlap.G2.1.G2.3 <- maxLikOverlap("2.1", "2.3", siber_cc_object, p = NULL, n =)
print(overlap.G2.1.G2.3)

##AT
##overlap between the SEAc of groups 3.1 and 3.2
overlap.G3.1.G3.2 <- maxLikOverlap("3.1", "3.2", siber_cc_object, p = NULL, n =)
print(overlap.G3.1.G3.2)
##overlap between the SEAc of groups 3.2 and 3.3
overlap.G3.2.G3.3 <- maxLikOverlap("3.2", "3.3", siber_cc_object, p = NULL, n =)
print(overlap.G3.2.G3.3)
##overlap between the SEAc of groups 3.1 and 3.3
overlap.G3.1.G3.3 <- maxLikOverlap("3.1", "3.3", siber_cc_object, p = NULL, n =)
print(overlap.G3.1.G3.3)

site_overlap_sea <- as.data.frame(rbind(overlap.G1.1.G1.2, overlap.G1.1.G1.3, overlap.G1.2.G1.3, overlap.G2.1.G2.2, overlap.G2.1.G2.3, overlap.G2.2.G2.3, overlap.G3.1.G3.2, overlap.G3.1.G3.3, overlap.G3.2.G3.3))
site_overlap_sea$comparison <- rownames(site_overlap_sea)

options(scipen = 999)
site_overlap_sea <- site_overlap_sea %>%
  mutate(prop_overlap_both = (overlap/(area.1+area.2))*100)
site_overlap_sea <- site_overlap_sea %>%
  mutate(prop_overlap_first = (overlap/(area.1))*100)
site_overlap_sea <- site_overlap_sea %>%
  mutate(prop_overlap_second = (overlap/(area.2))*100)

dev.off()
site_overlap_sea_plot <- ggplot(site_overlap_sea, aes(x = comparison)) +
  geom_point(aes(y = prop_overlap_both)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) 


site_overlap_sea_plot


###Comparing overlap of ellipses using Bayesian posterior distribution -----------------
#currently w/ 10 draws from bayesian distribution, how many should i use?
par(mfrow = c(1,1))
plotSiberObject(siber_cc_object,
                ax.pad = 2,
                hulls = F, community.hulls.args = community.hulls.args,
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = F, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = "%Terrestrial Coupling",
                ylab = "Trophic Position"
)
legend("topright", colnames(group.ML), 
       pch = c(1,1,1,2,2,2,3,3,3), col = c(1:3, 1:3, 1:3), lty = 1)


##HC  
##overlap b/w 1.1 and 1.2
bayes.overlap.G1.1.G1.2 <- bayesianOverlap("1.1", "1.2", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360, do.plot = TRUE)
bayes.overlap.G1.1.G1.2$comparison <- "G1.1.G1.2"
##overlap b/w 1.1 and 1.3
bayes.overlap.G1.1.G1.3 <- bayesianOverlap("1.1", "1.3", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360, do.plot = TRUE)
bayes.overlap.G1.1.G1.3$comparison <- "G1.1.G1.3"
##overlap b/w 1.2 and 1.3
bayes.overlap.G1.2.G1.3 <- bayesianOverlap("1.2", "1.3", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360)
bayes.overlap.G1.2.G1.3$comparison <- "G1.2.G1.3"

##EP4
##overlap b/w 2.1 and 2.2
bayes.overlap.G2.1.G2.2 <- bayesianOverlap("2.1", "2.2", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360)
bayes.overlap.G2.1.G2.2$comparison <- "G2.1.G2.2"
##overlap b/w 2.1 and 2.3
bayes.overlap.G2.1.G2.3 <- bayesianOverlap("2.1", "2.3", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360)
bayes.overlap.G2.1.G2.3$comparison <- "G2.1.G2.3"
##overlap b/w 2.2 and 2.3
bayes.overlap.G2.2.G2.3 <- bayesianOverlap("2.2", "2.3", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360)
bayes.overlap.G2.2.G2.3$comparison <- "G2.2.G2.3"

##AT
##overlap b/w 3.1 and 3.2
bayes.overlap.G3.1.G3.2 <- bayesianOverlap("3.1", "3.2", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360)
bayes.overlap.G3.1.G3.2$comparison <- "G3.1.G3.2"
##overlap b/w 3.1 and 3.3
bayes.overlap.G3.1.G3.3 <- bayesianOverlap("3.1", "3.3", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360)
bayes.overlap.G3.1.G3.3$comparison <- "G3.1.G3.3"
##overlap b/w 3.2 and 3.3
bayes.overlap.G3.2.G3.3 <- bayesianOverlap("3.2", "3.3", ellipses.posterior, 
                                           draws = 10, p.interval = NULL,
                                           n = 360)
bayes.overlap.G3.2.G3.3$comparison <- "G3.2.G3.3"

bayesian_overlap_all <- as.data.frame(rbind(bayes.overlap.G1.1.G1.2, bayes.overlap.G1.1.G1.3, bayes.overlap.G1.2.G1.3, bayes.overlap.G2.1.G2.2, bayes.overlap.G2.1.G2.3, bayes.overlap.G2.2.G2.3, bayes.overlap.G3.1.G3.2, bayes.overlap.G3.1.G3.3, bayes.overlap.G3.2.G3.3))
#bayesian_overlap_all$comparison <- rownames(bayesian_overlap_all)

bayesian_overlap_all$prop_overlap <- as.numeric(bayesian_overlap_all$overlap / (bayesian_overlap_all$area1 + bayesian_overlap_all$area2))

bayesian_overlap_all_plot <- ggplot(bayesian_overlap_all, aes(y = prop_overlap, x = comparison)) +
  geom_boxplot() +
  theme_classic()
bayesian_overlap_all_plot


