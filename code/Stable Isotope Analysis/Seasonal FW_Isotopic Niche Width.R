##Seasonal Stream Food Web Isotope Analysis - Isotopic Niche Width using SIBER on raw isotope data

##

set.seed(1)
library(SIBER)

source("code/Stable Isotope Analysis/Seasonal FW_Bi-Plots.R")

##Set up dataframe to prep for SIBER ------------
##Remove outliers (see Seasonal FW_Baseline Outliers.R code to see how outliers were determined)

SI_cc_al_dt_adj <- SI_cc_al_dt_adj %>%
  filter(Sample_ID != "AT1-18-AL-01") %>%
  filter(Sample_ID != "AT1-18-DT-01")

SI_cc_al_dt_adj <- SI_cc_al_dt_adj %>%
  mutate(group = case_when(
    startsWith(Sample_ID, "HC1") ~ "1",
    startsWith(Sample_ID, "AT1") ~ "1",
    startsWith(Sample_ID, "EP41") ~ "1",
    startsWith(Sample_ID, "HC2") ~ "2",
    startsWith(Sample_ID, "AT2") ~ "2",
    startsWith(Sample_ID, "EP42") ~ "2",
    startsWith(Sample_ID, "HC3") ~ "3",
    startsWith(Sample_ID, "AT3") ~ "3",
    startsWith(Sample_ID, "EP43") ~ "3",
  ))

SI_cc_al_dt_adj <- SI_cc_al_dt_adj %>%
  mutate(community = case_when(
    startsWith(Sample_ID, "HC") ~ "1",
    startsWith(Sample_ID, "AT") ~ "3",
    startsWith(Sample_ID, "EP4") ~ "2",
  ))

##Set up dataframe so formatted for SIBER
siber_cc <- SI_cc_al_dt_adj %>%
  filter(Species == "Creek Chub") %>%
  select(output_dHadjust_cc, d15N, group, community)

names(siber_cc)[1] <- "iso1" 
names(siber_cc)[2] <- "iso2"
siber_cc$group <- as.integer(siber_cc$group)
siber_cc$community <- as.integer(siber_cc$community)
typeof(siber_cc$iso1)

siber_cc <- as.data.frame(siber_cc) %>%
  arrange(community)
str(siber_cc)
##order community in ascending order -- seeing if this fixes issue w/ ellipse list names... 
#siber_cc$community <- ordered(siber_cc$community, 
#                    levels = c("1", "2", "3"))

##Testing normal distribution
library(mvnormtest)
shapiro.test(siber_cc$iso1)
shapiro.test(siber_cc$iso2)

##Create siber object ----------
siber_cc_object <- createSiberObject(siber_cc)

##Plotting SEA -------------
#Create lists of plotting arguments to be passed onwards to each of the three plotting functions
community.hulls.args <- list(col = 1, lty = 1, lwd = 1) 
group.ellipses.args <- list(n = 100, p.interval = NULL, lty = 1, lwd = 2) ## = NULL creates Maximum Liklihood Standard Ellipse (incorporates approx. 40% of data)
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
                xlab = expression({delta}^2*H~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)


group.ML <- groupMetricsML(siber_cc_object)
print(group.ML)

legend("topright", colnames(group.ML), 
       pch = c(1,1,1,2,2,2,3,3,3), col = c(1:3, 1:3, 1:3), lty = 1)



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

hc1_test <- ellipse_range_function(hc1)

ellipse_range <- split(siber_cc, paste0(siber_cc$community, siber_cc$group)) %>%
  map(ellipse_range_function) %>%
  bind_rows()
ellipse_range$community.group <- paste(ellipse_range$community, ".", ellipse_range$group)

xrange_plot <- ggplot(ellipse_range, aes(x = community.group, y = ellipse_xrange)) +
  geom_point() +
  theme_classic() +
  ylab("d2H Range")
xrange_plot

yrange_plot <- ggplot(ellipse_range, aes(x = community.group, y = ellipse_yrange)) +
  geom_point() +
  theme_classic() +
  ylab("d15N Range")
yrange_plot

##Can add other types of ellipse
# You can add more ellipses by directly calling plot.group.ellipses()
# Add an additional p.interval % prediction ellilpse
plotGroupEllipses(siber_cc_object, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# or you can add the XX% confidence interval around the bivariate means
# by specifying ci.mean = T along with whatever p.interval you want.
plotGroupEllipses(siber_cc_object, n = 100, p.interval = 0.95, ci.mean = T,
                  lty = 1, lwd = 2)

#Make ellipse plot in ggplot
siber_cc_ggplot <- siber_cc %>% mutate(group = factor(group), 
                                        community = factor(community),
                                        d2H = iso1, 
                                        d15N = iso2,
                                        .keep = "unused")

first.plot <- ggplot(data = siber_cc_ggplot, 
                     aes(x = d2H, 
                         y = d15N)) + 
  geom_point(aes(color = group, shape = community), size = 5) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{2}, "H (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_viridis_d() + 
  theme(panel.grid.minor = element_line(colour="white", size=0.5)) +
  scale_y_continuous(minor_breaks = seq(0 , 20, 2.5), breaks = seq(0, 20, 5))
first.plot

# decide how big an ellipse you want to draw
p.ell <- 0.40 ##creates predictive ellipse that contains approx. p % of data
ellipse.plot <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community), 
                   fill = group, 
                   color = group), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_viridis_d()
ellipse.plot

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

##Calculating point Layman metrics on raw isotope data (i.e., not bayesian distribution) ---------
##Calculate layman metrics for community 
community.ML <- as.data.frame(communityMetricsML(siber_cc_object))
print(community.ML)
community.ML$metrics <- rownames(community.ML)
names(community.ML)[1] <- "HC"
names(community.ML)[2] <- "EP4"
names(community.ML)[3] <- "AT"
community.ML_plot <- ggplot(community.ML, aes(x = metrics)) +
  geom_point(aes(y = HC, colour = "HawkCliff")) +
  geom_point(aes(y = EP4, colour = "Perl4")) +
  geom_point(aes(y = AT, colour = "Taylor")) +
  theme_classic() +
  ylab("Value")
community.ML_plot

##Calculate layman metrics for each group --  
hc.1 <- siber_cc %>%
  filter(group == "1", community == "1")
hc.1_layman <-laymanMetrics(hc.1$iso1, hc.1$iso2)
hc.2 <- siber_cc %>%
  filter(group == "2", community == "1")
hc.2_layman <-laymanMetrics(hc.2$iso1, hc.2$iso2)
hc.3 <- siber_cc %>%
  filter(group == "3", community == "1")
hc.3_layman <-laymanMetrics(hc.3$iso1, hc.3$iso2)

ep.1 <- siber_cc %>%
  filter(group == "1", community == "2")
ep.1_layman <-laymanMetrics(ep.1$iso1, ep.1$iso2)
ep.2 <- siber_cc %>%
  filter(group == "2", community == "2")
ep.2_layman <-laymanMetrics(ep.2$iso1, ep.2$iso2)
ep.3 <- siber_cc %>%
  filter(group == "3", community == "2")
ep.3_layman <-laymanMetrics(ep.3$iso1, ep.3$iso2)

at.1 <- siber_cc %>%
  filter(group == "1", community == "3")
at.1_layman <-laymanMetrics(at.1$iso1, at.1$iso2)
at.2 <- siber_cc %>%
  filter(group == "2", community == "3")
at.2_layman <-laymanMetrics(at.2$iso1, at.2$iso2)
at.3 <- siber_cc %>%
  filter(group == "3", community == "3")
at.3_layman <-laymanMetrics(at.3$iso1, at.3$iso2)

all_layman <- rbind(hc.1_layman[["metrics"]], hc.2_layman[["metrics"]], hc.3_layman[["metrics"]], ep.1_layman[["metrics"]], ep.2_layman[["metrics"]], ep.3_layman[["metrics"]], at.1_layman[["metrics"]], at.2_layman[["metrics"]], at.3_layman[["metrics"]])
all_layman <- as.data.frame(t(all_layman))
all_layman$metrics <- rownames(all_layman)

group_layman_plot <- ggplot(all_layman, aes(x = metrics)) +
  geom_point(aes(y = V1, colour = "HawkCliff"))+
  geom_point(aes(y = V2, colour = "HawkCliff")) +
  geom_point(aes(y = V3, colour = "HawkCliff")) +
  geom_point(aes(y = V4, colour = "Perl4")) +
  geom_point(aes(y = V5, colour = "Perl4")) +
  geom_point(aes(y = V6, colour = "Perl4")) +
  geom_point(aes(y = V7, colour = "Taylor")) +
  geom_point(aes(y = V8, colour = "Taylor")) +
  geom_point(aes(y = V9, colour = "Taylor")) +
  theme_classic() +
  ylab("Value")

group_layman_plot

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

##Plot out ellipses ---
library(ellipse)
# how many of the posterior draws do you want?
n.posts <- 10

# decide how big an ellipse you want to draw
#p.ell <- 0.95

# for a standard ellipse use
p.ell <- pchisq(1,2)

# a list to store the results
all_ellipses <- list()

# loop over groups
for (i in 1:length(ellipses.posterior)){
  
  # a dummy variable to build in the loop
  ell <- NULL
  post.id <- NULL
  
  for ( j in 1:n.posts){
    
    # covariance matrix
    Sigma  <- matrix(ellipses.posterior[[i]][j,1:4], 2, 2)
    
    # mean
    mu     <- ellipses.posterior[[i]][j,5:6]
    
    # ellipse points
    
    out <- ellipse::ellipse(Sigma, centre = mu , level = p.ell)
    
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
    
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

ellipse_df <- bind_rows(all_ellipses, .id = "id")


# now we need the group and community names

# extract them from the ellipses.posterior list
group_comm_names <- names(ellipses.posterior)[as.numeric(ellipse_df$id)]

# split them and conver to a matrix, NB byrow = T
split_group_comm <- matrix(unlist(strsplit(group_comm_names, "[.]")),
                           nrow(ellipse_df), 2, byrow = TRUE)

ellipse_df$community <- split_group_comm[,1]
ellipse_df$group     <- split_group_comm[,2]

ellipse_df <- dplyr::rename(ellipse_df, iso1 = x, iso2 = y)

##First plot raw data 
first.plot <- ggplot(data = siber_cc, aes(iso1, iso2)) +
  geom_point(aes(color = factor(community):factor(group)), size = 2)+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{2}, "H (\u2030)"))) + 
  theme(text = element_text(size=15))
print(first.plot)

second.plot <- first.plot + facet_wrap(~factor(community):factor(group))
print(second.plot)

third.plot <- second.plot + 
  geom_polygon(data = ellipse_df,
               mapping = aes(iso1, iso2,
                             group = rep,
                             color = factor(community):factor(group),
                             fill = NULL),
               fill = NA,
               alpha = 0.2)
print(third.plot)
#Calculate some credible intervals
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles) ##get credible interval for each group

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes)


##SEAb size -- Comparing the posterior distributions of Bayesian Ellipses ---------------
##SEAb ellipse comparison within group
##HawkCliff
##comparing probability of group 1.1 being smaller than group 1.2
Pg1.1_lt_g1.2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.1_lt_g1.2)
##comparing probability of group 1.2 being smaller than group 1.3
Pg1.2_lt_g1.3 <- sum( SEA.B[,2] < SEA.B[,3] ) / nrow(SEA.B)
print(Pg1.2_lt_g1.3)
##comparing probability of group 1.1 being smaller than group 1.3
Pg1.1_lt_g1.3 <- sum( SEA.B[,1] < SEA.B[,3] ) / nrow(SEA.B)
print(Pg1.1_lt_g1.3)

##EP4
##comparing probability of group 2.1 being smaller than group 2.2
Pg2.1_lt_g2.2 <- sum( SEA.B[,4] < SEA.B[,5] ) / nrow(SEA.B)
print(Pg2.1_lt_g2.2)
##comparing probability of group 2.2 being smaller than group 2.3
Pg2.2_lt_g2.3 <- sum( SEA.B[,5] < SEA.B[,6] ) / nrow(SEA.B)
print(Pg2.2_lt_g2.3)
##comparing probability of group 2.1 being smaller than group 2.3
Pg2.1_lt_g2.3 <- sum( SEA.B[,4] < SEA.B[,6] ) / nrow(SEA.B)
print(Pg2.1_lt_g2.3)

##Taylor
##comparing probability of group 3.1 being smaller than group 3.2
Pg3.1_lt_g3.2 <- sum( SEA.B[,7] < SEA.B[,8] ) / nrow(SEA.B)
print(Pg3.1_lt_g3.2)
##comparing probability of group 3.2 being smaller than group 3.3
Pg3.2_lt_g3.3 <- sum( SEA.B[,8] < SEA.B[,9] ) / nrow(SEA.B)
print(Pg3.2_lt_g3.3)
##comparing probability of group 3.1 being smaller than group 3.3
Pg3.1_lt_g3.3 <- sum( SEA.B[,7] < SEA.B[,9] ) / nrow(SEA.B)
print(Pg3.1_lt_g3.3)

site_ellipse_size_diff <- as.data.frame(rbind(Pg1.1_lt_g1.2, Pg1.2_lt_g1.3, Pg1.1_lt_g1.3, Pg2.1_lt_g2.2, Pg2.1_lt_g2.3, Pg2.2_lt_g2.3, Pg3.1_lt_g3.2, Pg3.1_lt_g3.3, Pg3.2_lt_g3.3))
site_ellipse_size_diff$comparison <- rownames(site_ellipse_size_diff)


##SEAb ellipse comparison across group, within season
##Spring
##comparing probability of group 1.1 being smaller than group 2.1
Pg1.1_lt_g2.1 <- sum( SEA.B[,1] < SEA.B[,4] ) / nrow(SEA.B)
print(Pg1.1_lt_g2.1)
##comparing probability of group 2.1 being smaller than group 3.1
Pg2.1_lt_g3.1 <- sum( SEA.B[,4] < SEA.B[,7] ) / nrow(SEA.B)
print(Pg2.1_lt_g3.1)
##comparing probability of group 1.1 being smaller than group 3.1
Pg1.1_lt_g3.1 <- sum( SEA.B[,1] < SEA.B[,7] ) / nrow(SEA.B)
print(Pg1.1_lt_g3.1)

##Summer
##comparing probability of group 1.2 being smaller than group 2.2
Pg1.2_lt_g2.2 <- sum( SEA.B[,2] < SEA.B[,5] ) / nrow(SEA.B)
print(Pg1.2_lt_g2.2)
##comparing probability of group 2.2 being smaller than group 3.2
Pg2.2_lt_g3.2 <- sum( SEA.B[,5] < SEA.B[,8] ) / nrow(SEA.B)
print(Pg2.2_lt_g3.2)
##comparing probability of group 1.2 being smaller than group 3.2
Pg1.2_lt_g3.2 <- sum( SEA.B[,2] < SEA.B[,8] ) / nrow(SEA.B)
print(Pg1.2_lt_g3.2)

##Fall
##comparing probability of group 1.3 being smaller than group 2.3
Pg1.3_lt_g2.3 <- sum( SEA.B[,3] < SEA.B[,7] ) / nrow(SEA.B)
print(Pg1.3_lt_g2.3)
##comparing probability of group 2.3 being smaller than group 3.3
Pg2.3_lt_g3.3 <- sum( SEA.B[,7] < SEA.B[,9] ) / nrow(SEA.B)
print(Pg2.3_lt_g3.3)
##comparing probability of group 1.3 being smaller than group 3.3
Pg1.3_lt_g3.3 <- sum( SEA.B[,3] < SEA.B[,9] ) / nrow(SEA.B)
print(Pg1.3_lt_g3.3)

seasonal_ellipse_size_diff <- as.data.frame(rbind(Pg1.1_lt_g2.1, Pg1.1_lt_g3.1, Pg2.1_lt_g3.1, Pg1.2_lt_g2.2, Pg1.2_lt_g3.2, Pg2.2_lt_g3.2, Pg1.3_lt_g2.3, Pg1.3_lt_g3.3, Pg2.3_lt_g3.3))
seasonal_ellipse_size_diff$comparison <- rownames(seasonal_ellipse_size_diff)

seasonal_ellipse_size_diff_plot <- ggplot(seasonal_ellipse_size_diff, aes(x = comparison, y = V1)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Proportion of Group1 < Group2")

  
seasonal_ellipse_size_diff_plot ##note: this is not in seasonal/site order
##can you compare more than 2 ellipses at a time?

##Comparing communities using Layman's Metrics on Bayesian distributions --------
# extract the posterior means
mu.post <- extractPosteriorMeans(siber_cc_object, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

LB_1.1 <- layman.B[[1.1]]
LB_1.2 <- layman.B[[1.2]]

##Community 1 -- HC
par(mfrow = c(3,1))
siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,10))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber_cc_object$ML.mu[[1]][1,1,],
                                 siber_cc_object$ML.mu[[1]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)

##Community 2 -- EP4
siberDensityPlot(layman.B[[2]], xticklabels = colnames(layman.B[[2]]), 
                 bty="L", ylim = c(0,10))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber_cc_object$ML.mu[[2]][1,1,],
                                 siber_cc_object$ML.mu[[2]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)

##Community 3 -- AT
siberDensityPlot(layman.B[[3]], xticklabels = colnames(layman.B[[3]]), 
                 bty="L", ylim = c(0,10))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber_cc_object$ML.mu[[3]][1,1,],
                                 siber_cc_object$ML.mu[[3]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)

##can i calculate the bayesian layman's metrics on each group within the community?

##Determining Ellipse overlap of maximum liklihood ellipses -------------
##note: the p = 0.95 is from the link below, how do i know if this is using the SEAc? is it b/c of the maxLikOverlap function? no, based on the r documentation, the p should be determined same as in addEllipse function, where if p = NULL that plots the SEA
##https://cran.r-project.org/web/packages/SIBER/vignettes/siber-comparing-populations.html
##When p = NULL, this uses area of SEA, how do i change it so it is area of SEAc? 
##HawkCliff
##overlap between the SEAc of groups 1.2 and 1.3
overlap.G1.2.G1.3 <- maxLikOverlap("1.2", "1.3", siber_cc_object, p = NULL, n = 100, do.plot = TRUE)
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
##can i calculate overlap of 3 ellipses? i feel like i should be able to but haven't found code as to how


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
                xlab = expression({delta}^2*H~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)
legend("topright", colnames(group.ML), 
       pch = c(1,1,1,2,2,2,3,3,3), col = c(1:3, 1:3, 1:3), lty = 1)


##HC  
##overlap b/w 1.1 and 1.2
bayes.overlap.G1.1.G1.2 <- bayesianOverlap("1.1", "1.2", ellipses.posterior, 
                                       draws = 10, p.interval = NULL,
                                       n = 100, do.plot = TRUE)
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


###Calculate Layman metrics of Bayesian distributions for each group ---------------
##Calculate layman metrics for each group --  
bayesian_group_laymans <- bayesianLayman(ellipses.posterior)





hc.1 <- siber_cc %>%
  filter(group == "1", community == "1")
hc.1_layman <-laymanMetrics(hc.1$iso1, hc.1$iso2)
hc.2 <- siber_cc %>%
  filter(group == "2", community == "1")
hc.2_layman <-laymanMetrics(hc.2$iso1, hc.2$iso2)
hc.3 <- siber_cc %>%
  filter(group == "3", community == "1")
hc.3_layman <-laymanMetrics(hc.3$iso1, hc.3$iso2)

ep.1 <- siber_cc %>%
  filter(group == "1", community == "2")
ep.1_layman <-laymanMetrics(ep.1$iso1, ep.1$iso2)
ep.2 <- siber_cc %>%
  filter(group == "2", community == "2")
ep.2_layman <-laymanMetrics(ep.2$iso1, ep.2$iso2)
ep.3 <- siber_cc %>%
  filter(group == "3", community == "2")
ep.3_layman <-laymanMetrics(ep.3$iso1, ep.3$iso2)

at.1 <- siber_cc %>%
  filter(group == "1", community == "3")
at.1_layman <-laymanMetrics(at.1$iso1, at.1$iso2)
at.2 <- siber_cc %>%
  filter(group == "2", community == "3")
at.2_layman <-laymanMetrics(at.2$iso1, at.2$iso2)
at.3 <- siber_cc %>%
  filter(group == "3", community == "3")
at.3_layman <-laymanMetrics(at.3$iso1, at.3$iso2)

all_layman <- rbind(hc.1_layman[["metrics"]], hc.2_layman[["metrics"]], hc.3_layman[["metrics"]], ep.1_layman[["metrics"]], ep.2_layman[["metrics"]], ep.3_layman[["metrics"]], at.1_layman[["metrics"]], at.2_layman[["metrics"]], at.3_layman[["metrics"]])
all_layman <- as.data.frame(t(all_layman))
all_layman$metrics <- rownames(all_layman)

group_layman_plot <- ggplot(all_layman, aes(x = metrics)) +
  geom_point(aes(y = V1, colour = "HawkCliff"))+
  geom_point(aes(y = V2, colour = "HawkCliff")) +
  geom_point(aes(y = V3, colour = "HawkCliff")) +
  geom_point(aes(y = V4, colour = "Perl4")) +
  geom_point(aes(y = V5, colour = "Perl4")) +
  geom_point(aes(y = V6, colour = "Perl4")) +
  geom_point(aes(y = V7, colour = "Taylor")) +
  geom_point(aes(y = V8, colour = "Taylor")) +
  geom_point(aes(y = V9, colour = "Taylor")) +
  theme_classic() +
  ylab("Value")

group_layman_plot