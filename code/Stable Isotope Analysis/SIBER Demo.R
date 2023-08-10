##SIBER Demo
##url: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html
set.seed(1)
library(SIBER)
data("demo.siber.data")
##Create siber object - convert this raw data form into an object that contains information 
#summary statistics that are useful for plotting, and a z-score transformed version of the data 
#which is used in the model fitting process before being back-transformed using the summary statistic information.
siber_example <- createSiberObject(demo.siber.data) 


#Create lists of plotting arguments to be passed onwards to each of the three plotting functions
##Convex hulls are drawn between the centres of each group within a community with hulls = T
##Ellipses are drawn around each group independently w/ ellipses = T. These ellipses can be made to be maximum likelihood standard ellipses by setting p = NULL, or can be made to be prediction ellipses that contain approximately p proportion of data. For example, p = 0.95 will draw an ellipse that encompasses approximately 95% of the data. The parameter n determines how many points are used to make each ellipse and hence how smooth the curves are
##Convext hulls are drawn around each group independently w/ group.hulls = T
community.hulls.args <- list(col = 1, lty = 1, lwd = 1) 
group.ellipses.args <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args <- list(lty = 2, col = "grey20") 

par(mfrow = c(1,1))
plotSiberObject(siber_example,
                ax.pad = 2,
                hulls = F, community.hulls.args = community.hulls.args,
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = F, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
                )

##Calculating summary statistics: TA, SEA, SEAc
##This function loops over each group within each community and calculates the convex hull total area, Standard Ellipse Area (SEA) and its corresponding small sample size corrected version SEAc based on the maximum likelihood estimates of the means and covariance matrices of each group.
group.ML <- groupMetricsML(siber_example)
print(group.ML)

###Fitting Bayesian models to the data
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
ellipses.posterior <- siberMVN(siber_example, parms, priors)


##Comparing among groups using Standard Ellipse Area
#When comparing individual groups with each other, be it within a single community, or groups among communities, 
#the Standard Ellipse Area (SEA) is the recommended method. Since the multivariate normal distributions have already been fitted to each group, 
#it only remains to calculate the SEA on the posterior distribution of covariance matrix for each group, thereby yielding the Bayesian SEA or SEA-B.
# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
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

