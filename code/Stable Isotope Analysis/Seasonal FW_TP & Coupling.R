##Code: Seasonal Stream Food Web Stable Isotope Analysis -- Calculate Trophic Position and %coupling
##By: Marie Gutgesell
##Date: May 4, 2022

source("~/Desktop/Seasonal Stream Food Webs/Code/Stable Isotope Analysis/Seasonal FW_Bi-Plots.R")


##Remove baseline outliers (see Seasonal FW_Baseline Outliers.R code to see how outliers were determined)

SI_cc_al_dt_adj <- SI_cc_al_dt_adj %>%
 filter(Sample_ID != "AT1-18-AL-01") %>%
  filter(Sample_ID != "AT1-18-DT-01")


##In Champagne et al., 2022 -- calculated TP using 1-source (terrestrial baselines) as main result, and two-baseline method to corroborate results in supplement -- likely b/c weird algae values that gave 0 or negative TPs, which is not issue here
##Calculating TP & dN using bound %TC b/w 0-1 --------
si_mixing_model_pred_bound <- function(x){
  df_prep <- x
  AB <- df_prep %>%
    filter(Species == "Algae")
  AB_mean_dH <- mean(AB$d2H)
  AB_mean_dN <- mean(AB$d15N)
  TB <- df_prep %>%
    filter(Species == "Detritus")
  TB_mean_dH <- mean(TB$d2H)
  TB_mean_dN <- mean(TB$d15N)
  df_pred <- df_prep %>%
    filter(Species =="Creek Chub")
  for(i in 1:nrow(df_pred)){
    TC_pred = ((df_pred$output_dHadjust_cc - AB_mean_dH)/(TB_mean_dH - AB_mean_dH))
    TC_pred[TC_pred>=1] <- 0.999
    TC_pred[TC_pred<=0] <- 0.001
    TP_pred = 1 + (df_pred$d15N - ((AB_mean_dN*(1 - TC_pred)) + (TB_mean_dN*TC_pred)))/3.4
    dN_pred = df_pred$d15N - (((1-TC_pred)*AB_mean_dN)+(TC_pred*TB_mean_dN))
    TC_TP_df_pred <- cbind(TC_pred, TP_pred, dN_pred)
  }
  TC_TP_df_pred <- as.data.frame(TC_TP_df_pred)
  TC_TP_df_pred$Species <- df_pred$Species
  TC_TP_df_pred$Site <- df_pred$Site
  TC_TP_df_pred$Season <- df_pred$Season
  TC_TP_df_pred$Sample_ID <- df_pred$Sample_ID
  TC_TP_df_pred
}

seasonal_df_split_pred_tp_bound <- split(SI_cc_al_dt_adj, paste0(SI_cc_al_dt_adj$Site, SI_cc_al_dt_adj$Season)) %>%
  map(si_mixing_model_pred_bound) %>%
  bind_rows()


tp_dn <- left_join(seasonal_df_split_pred_tp_bound, SI_cc_al_dt_adj, by = "Sample_ID")

tp_dn_plot <- ggplot(tp_dn, aes(x = d15N, y = TP_pred, color = Site.x)) +
  geom_point()
tp_dn_plot


tc_dh_plot <- ggplot(tp_dn, aes(x = d2H, y = TC_pred, color = Site.x)) +
  geom_point()
tc_dh_plot


seasonal_df_split_pred_tp_bound$Season <- ordered(seasonal_df_split_pred_tp_bound$Season, 
                              levels = c("Spring", "Summer", "Fall"))

##TESTING ANOVA ASSUMPTIONS -----------------

##testing normality of TP & % TC for use in ANOVA and isotopic niche analysis
library(car)
library(rstatix)

seasonal_df_split_pred_tp_bound %>%
  group_by(Site, Season) %>%
  shapiro_test(TC_pred)

seasonal_df_split_pred_tp_bound %>%
  group_by(Site, Season) %>%
  shapiro_test(TP_pred) ##only Perl_4 spring is not normal -- p = 0.008


#ggqqplot(seasonal_df_split_pred_tp_bound$TC_pred)
hist(seasonal_df_split_pred_tp_bound$TC_pred)
#ggqqplot(seasonal_df_split_pred_tp_bound$TP_pred)
hist(seasonal_df_split_pred_tp_bound$TP_pred)


seasonal_df_split_pred_tp_bound <- seasonal_df_split_pred_tp_bound %>%
  mutate(TC_pred_log = log(TC_pred)) %>%
  mutate(TP_pred_log = log(TP_pred)) %>%
  mutate(TC_pred_logit = logit(TC_pred)) %>% ##logit transformation appropriate for percentages/proportions
  mutate(TC_pred_qlogis = qlogis(TC_pred)) %>%
  mutate(TP_pred_sqrt = sqrt(TP_pred))

hist(seasonal_df_split_pred_tp_bound$TP_pred_sqrt)
hist(seasonal_df_split_pred_tp_bound$TC_pred_logit)

##this is equivalent to the logit transformation

seasonal_df_split_pred_tp_bound %>%
  group_by(Site, Season) %>%
  shapiro_test(TC_pred_logit) ##does improve normality, only 2 sites non-normal hc spring and perl 4 fall

seasonal_df_split_pred_tp_bound %>%
  group_by(Site, Season) %>%
  shapiro_test(TP_pred_log) ###improves normality of EP4 spring -- almost normal (p = 0.04)


sd <- seasonal_df_split_pred_tp_bound %>%
  group_by(Site, Season) %>%
  summarise_each(funs(sd)) ##have different variances, but doesn't look like transformations improve it...

##Normality results: using logit transformed % coupling and not transformed TP 

##Testing for homogeneity of variance 
library(car)
leveneTest(TC_pred_logit ~ Site,data = seasonal_df_split_pred_tp_bound)
leveneTest(TP_pred ~ Site,data = seasonal_df_split_pred_tp_bound)
leveneTest(TP_pred_sqrt ~ Site,data = seasonal_df_split_pred_tp_bound)

qqPlot(seasonal_df_split_pred_tp_bound$TC_pred_logit)
qqPlot(seasonal_df_split_pred_tp_bound$TP_pred)


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



seasonal_df_split_pred_tp_bound_noout <- seasonal_df_split_pred_tp_bound %>% 
  group_by(Site, Season) %>% 
  mutate(TC_pred_logit = remove_outliers(TC_pred_logit)) %>% 
  mutate(TP_pred = remove_outliers(TP_pred)) %>%
  ungroup() %>% 
  filter(!is.na(TC_pred_logit)) %>%
  filter(!is.na(TP_pred))

seasonal_df_split_pred_tp_bound_noout <- seasonal_df_split_pred_tp_bound_noout%>%
  mutate(Site_Type = case_when(
    startsWith(Sample_ID, "HC") ~ "Low Impact",
    startsWith(Sample_ID, "EP4") ~ "Mid Impact",
    startsWith(Sample_ID, "AT") ~ "High Impact"
  ) )

seasonal_df_split_pred_tp_bound_noout$Site_Type <- ordered(seasonal_df_split_pred_tp_bound_noout$Site_Type, 
                                                  levels = c("Low Impact", "Mid Impact", "High Impact"))

mean_trophic_metrics <- seasonal_df_split_pred_tp_bound_noout %>%
  group_by(Site, Season) %>%
  summarise_all(mean)

coupling_boxplot <- ggplot(seasonal_df_split_pred_tp_bound_noout, aes(x = Season, y = TC_pred, fill = Site_Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme_classic() +
  ylab("% Terrestrial Energy Use") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),axis.text.y = element_text(size = 14),axis.title.y=element_text(size = 14), axis.title.x = element_blank(), legend.position = "right", text = element_text(family = "Times New Roman"))+
  guides(fill = guide_legend(title="Site"))


coupling_boxplot


TP_boxplot <- ggplot(seasonal_df_split_pred_tp_bound_noout, aes(x = Season, y = TP_pred, fill = Site_Type)) +
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  ylab("Trophic Position") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),axis.text.y = element_text(size = 10),axis.title.y=element_text(size = 12), axis.title.x = element_blank(), legend.position = "right", text = element_text(family = "Times New Roman"))+ 
  guides(fill = guide_legend(title="Site"))

TP_boxplot

##Create figure legend for 3 sites 
library(ggpubr)
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
  theme(legend.title = element_text(family = "Times New Roman"), legend.text = element_text(size = 12, family = "Times New Roman"), legend.position = "right")
gplot

leg_fig <- get_legend(gplot)

tp_coupling_plot <- ggarrange(coupling_boxplot, TP_boxplot,  legend = "right", common.legend = TRUE, legend.grob = leg_fig,
                                labels = c("a)", "b)"),
                                ncol = 2, nrow = 1, font.label = list(colour = "black", size = 12, family = "Times New Roman"))
tp_coupling_plot


##ANOVA of % coupling and trophic position

coupling_anova <- aov(TC_pred_logit ~ Site_Type + Season, seasonal_df_split_pred_tp_bound_noout)
summary(coupling_anova)

coupling_tukey <- TukeyHSD(coupling_anova)
coupling_tukey_table <- coupling_tukey[["Site_Type:Season"]]
write.csv(coupling_tukey_table, "~/Desktop/Seasonal Stream Food Webs/R output tables/coupling_tukey_table.csv")

tp_anova <- aov(TP_pred ~ Site_Type + Season, seasonal_df_split_pred_tp_bound_noout)
summary(tp_anova)
tp_tukey <- TukeyHSD(tp_anova)
tp_tukey_table <- tp_tukey[["Site_Type:Season"]]
write.csv(tp_tukey_table, "~/Desktop/Seasonal Stream Food Webs/R output tables/tp_tukey_table.csv")

