##Seasonal Stream Food Web Isotope Analysis - Isotopic Niche Width using nicheROVER on trophic position and % coupling metrics

##Using nicheRover

source("code/Stable Isotope Analysis/Seasonal FW_TP & Coupling.R")
#rm(seasonal_df_split_pred_tp_bound_noout, cc_seasonal_stream_si_plot_site_adj, coupling_boxplot, SI_all, SI_cc_al_dt_adj, TP_boxplot, leg_fig, tc_dh_plot, tp_dn, tp_dn_plot, gplot, data)

library(nicheROVER)
?niche.size
#data("fish")

###1) Format data for use in nicheRover (TP and %TC as columns, each row is an observation)
##Create column that combines site and season
seasonal_df_split_pred_tp_bound_noout <- seasonal_df_split_pred_tp_bound_noout%>%
  mutate(Group = case_when(
    startsWith(Sample_ID, "HC1") ~ "HC_spring",
    startsWith(Sample_ID, "HC2") ~ "HC_summer", 
    startsWith(Sample_ID, "HC3") ~ "HC_fall",
    startsWith(Sample_ID, "EP41") ~ "EP4_spring",
    startsWith(Sample_ID, "EP42") ~ "EP4_summer", 
    startsWith(Sample_ID, "EP43") ~ "EP4_fall",
    startsWith(Sample_ID, "AT1") ~ "AT_spring",
    startsWith(Sample_ID, "AT2") ~ "AT_summer", 
    startsWith(Sample_ID, "AT3") ~ "AT_fall"
  ) )

data <- seasonal_df_split_pred_tp_bound_noout %>%
  select(Group, TC_pred_logit, TP_pred) ####Note: have used data w/ outliers removed and logit transformed % TC

typeof(data$TC_pred_logit)
typeof(data$TP_pred)

str(data)
#data[, 2:3] <- sapply(data[, 2:3], as.numeric)

data$Group <- ordered(data$Group,
                                    levels = c("HC_spring", "HC_summer", "HC_fall", "EP4_spring", "EP4_summer", "EP4_fall", "AT_spring", "AT_summer", "AT_fall"))

aggregate(data[2:3], data[1], mean) ##calculates mean TP and %TC for each site:season

##Generate the posterior distributions of u (mean) and V (variance) for isotope values of each species with the default prior
# generate parameter draws from the "default" posteriors of each fish
nsamples <- 1e4
system.time({
  fish.par <- tapply(1:nrow(data), data$Group,
                     function(ii) niw.post(nsamples = nsamples, X = data[ii,2:3]))
})

##Plot mean and variance parameters
clrs <- c("black", "red", "blue", "orange", "green", "purple", "yellow", "grey", "pink") # colors for each species

# mu1 (%TC), mu2 (TP), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

##plot all mu and sigma
par(mar = c(4.2, 4.2, 2, 1)+.1)
niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = TRUE)
legend("topright", legend = names(fish.par), fill = clrs)

# 2-d projections of 10 niche regions
clrs <- c("black", "red", "blue", "orange", "green", "purple", "yellow", "grey", "pink") # colors for each species
nsamples <- 10
fish.par <- tapply(1:nrow(data), data$Group,
                   function(ii) niw.post(nsamples = nsamples, X = data[ii,2:3]))

# format data for plotting function
data_NR <- tapply(1:nrow(data), data$Group, function(ii) X = data[ii,2:3])

niche.plot(niche.par = fish.par, niche.data = data_NR, pfrac = .05,
           iso.names = expression("Terrestrial coupling", "Trophic position"),
           col = clrs, xlab = expression("Terrestrial Coupling"))

# niche overlap plots for 95% niche region sizes
nsamples <- 10000
fish.par <- tapply(1:nrow(data), data$Group,
                   function(ii) niw.post(nsamples = nsamples, X = data[ii,2:3]))

# Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher accuracy.
# the variable over.stat can be supplied directly to the overlap.plot function

over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, 0.40))

#The mean overlap metrics calculated across iteratations for both niche 
#region sizes (alpha = .95 and alpha = .99) can be calculated and displayed in an array.
over.mean <- apply(over.stat, c(1:2,4), mean)*100 
mean_overlap <- round(over.mean, 2)
mean_overlap<- as.data.frame(mean_overlap)
write.csv(mean_overlap, "outputs/mean_niche_overlap.csv")


##determine credible interval of overlap
over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
cred_overlap_95 <- round(over.cred[,,,1]) # display alpha = .95 niche region
cred_overlap_95 <- as.data.frame(cred_overlap_95)
write.csv(cred_overlap_95, "outputs/CI_niche_overlap_95.csv")

cred_overlap_40 <- round(over.cred[,,,2])
cred_overlap_40 <- as.data.frame(cred_overlap_40)
write.csv(cred_overlap_40, "outputs/CI_niche_overlap_40.csv")


##Plot 95% overlap
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e4, alpha = .95)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

##Determining size distribution of niche
# posterior distribution of (mu, Sigma) for each site*season
nsamples <- 10000
fish.par <- tapply(1:nrow(data), data$Group,
                   function(ii) niw.post(nsamples = nsamples, X = data[ii,2:3]))

# posterior distribution of niche size by species
p.ell <- pchisq(1,2) ##SEA 
p.ell2 <- 0.40 ##slightly larger mean niche area, which makes sense

fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = p.ell)
})

##have estimate of 10000 SEA for each site-season
##I guess could do SEAc calculation for each of these estimates, and get an estimate of SEAc.. 
# point estimate and standard error
rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

##mean SEA is smaller than SEA from SIBER -- this does incorporate uncertainty already
##SEA equation: SEAc = SEA x [(n-1)/(n-2)]

dev.off()

boxplot(fish.size, col = clrs,
        ylab = "Niche Size (SEA)", xlab = "Group")

##Set up fish.size dataframe to format figure for manuscript
fish.size <- as.data.frame(fish.size)
fish.size <- fish.size %>%
  gather(value)
colnames(fish.size)[1] <- "Group"

fish.size$Group <- ordered(fish.size$Group,
                      levels = c("HC_spring", "HC_summer", "HC_fall", "EP4_spring", "EP4_summer", "EP4_fall", "AT_spring", "AT_summer", "AT_fall"))

fish.size$Group <- as.character(fish.size$Group)
fish.size <- fish.size%>%
  mutate(Site = case_when(
    startsWith(Group, "HC") ~ "Low Impact",
    startsWith(Group, "EP4") ~ "Mid Impact", 
    startsWith(Group, "AT") ~ "High Impact"
  ))
fish.size <- fish.size%>%
  mutate(Season = case_when(
    grepl("spring", Group) ~ "Spring",
    grepl("summer", Group) ~ "Summer", 
    grepl("fall", Group) ~ "Fall"
  ))

fish.size$Season <- ordered(fish.size$Season,
                           levels = c("Spring", "Summer", "Fall"))
fish.size$Site <- ordered(fish.size$Site,
                            levels = c("Low Impact", "Mid Impact", "High Impact"))


SEA_plot <- ggplot(fish.size, aes(x = Season, y = value, fill = Site)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))  +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "right", text = element_text(family = "Times New Roman")) +
  ylab("SEA") 
  
SEA_plot
##Looking at overlap of of 3 seasons of each site separately 

SEA_plot_mean <- ggplot(fish.size, aes(x = Site, y = value, fill = Site)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"))  +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "right", text = element_text(family = "Times New Roman")) +
  ylab("SEA") 

SEA_plot_mean


##ANOVA of SEAc -- seasonal and mean
sea_anova_seasonal <- aov(value ~ Site * Season, data = fish.size)
summary(sea_anova_seasonal)
sea_tukey_seasonal <- TukeyHSD(sea_anova_seasonal)
sea_tukey_seasonal[["Site:Season"]]

sea_anova_mean <- aov(value ~ Site, data = fish.size)
summary(sea_anova_mean)
sea_tukey_mean <- TukeyHSD(sea_anova_mean)
sea_tukey_mean[["Site"]]

#Make ellipse plot in ggplot
seasonal_df_split_pred_tp_bound_noout <- seasonal_df_split_pred_tp_bound_noout%>%
  mutate(Site_Type = case_when(
    startsWith(Group, "HC") ~ "Low Impact",
    startsWith(Group, "EP4") ~ "Mid Impact", 
    startsWith(Group, "AT") ~ "High Impact"
  ))


data_2 <- seasonal_df_split_pred_tp_bound_noout %>%
  select(Site_Type, Season, TC_pred_logit, TP_pred)

data_2$Site_Type <- ordered(data_2$Site_Type,
                          levels = c("Low Impact", "Mid Impact", "High Impact"))


first.plot <- ggplot(data = data_2, 
                     aes(x = TC_pred_logit, 
                         y = TP_pred)) + 
  geom_point(aes(color = Site_Type, shape = Season), size = 5) +
  ylab("Trophic Position") +
  xlab("% Terrestrial Energy Use (logit)") + 
  theme(text = element_text(size=16)) + 
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) + 
  theme_classic()
first.plot

# decide how big an ellipse you want to draw
#p.ell <- pchisq(1,2) ##creates predictive ellipse that contains approx. p % of data
p.ell <- pchisq(1,2)
ellipse.plot <- first.plot + 
  stat_ellipse(aes(group = interaction(Season, Site_Type), 
                   fill = Site_Type,
                   color = Site_Type), 
               alpha = 0.25, 
               level = 0.40,
               type = "norm",
               geom = "polygon") +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.text.x = element_text(hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_text(size = 12), legend.position = "right", text = element_text(family = "Times New Roman")) 

ellipse.plot

##whether using pchisq(1,2) or p = 0.40 essentially the same just p 0.40 is slightly larger, because pchisq(1,2) = ~0.393

##then try and correct SEA for SEAc -- see if any different, and if possible to use SEAc for overlap? 
##or try using overlap again in siber... 
