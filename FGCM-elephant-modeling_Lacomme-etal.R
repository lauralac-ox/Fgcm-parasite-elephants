# R Scirpt for Manuscript: "Conservation strategies shape stress-related
# glucocorticoid levels of managed elephant populations more than parasitic infection"

# Authors: Laura Lacomme, Hervé Fritz, André Ganswindt, Chloé Guerbois, Benjamin Rey
# Corresponding Authors: Laura Lacomme. Email: laura.lacomme.pro@gmail.com
# Date: 2026

# R version 4.4.2

# SET-UP ----
# Load required packages
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(grid) 
library(hrbrthemes)
library(ggpubr)
library(car)
library(lme4)
library(lmerTest)
library(tidyverse)
library(tidytext)
library(dplyr)
library(ggResidpanel)
library(performance)
library(see)
library(effects)
library(jtools)
library(interactions)
library(glmmTMB)
library(MASS)
library(ggpmisc)
library(plyr)
library(labelled)
library(gtsummary)
library(AICcmodavg)
library(flextable)

# Load and clean dataset
setwd(dir = "C:/Users/biol0277/OneDrive - Nexus365/3. PHD MANUSCRIPTS/5. STAT ANALYSES/Session 1/data/Repository")
dt <- read.csv("dataset_fgcm_modeling.csv")
dt <- lapply(dt, function(col) {
  if (is.character(col)) {
    as.factor(col)
  } else {
    col
  }
})
dt <- as.data.frame(dt)
dt$campaign <- as.factor(dt$campaign)
dt$epg <- as.numeric(dt$epg)


# CHECK DATA  ----

# Histogram and QQ plot of FGCMs
par(mar=c(5,4,4,2)+0.1)

# Histogram
x2 <- seq(min(dt$hf_ng), max(dt$hf_ng), length = 50)
fun <- dnorm(x2, mean = mean(dt$hf_ng), sd = sd(dt$hf_ng))
hist(dt$hf_ng,
     #breaks=50,
     prob=TRUE,
     col="white")
lines(x2, fun, col = 2, lwd = 2)
lines(density(dt$hf_ng), col = 4, lwd = 2)

# QQ plot
qqnorm(dt$hf_ng, datax=TRUE, main="FGCMs")
qqline(dt$hf_ng,datax=TRUE)

# Stats test
t.test(dt$hf_ng) #hyp null (mean=0) is rejected
shapiro.test(dt$hf_ng) #hypothesis of normality is rejected
ks.test(dt$hf_ng,"pnorm",mean(dt$hf_ng),sd(dt$hf_ng)) #hypothesis of normality is rejected

# Stats test
t.test(dt$hf_ng) #hyp null (mean=0) is rejected
shapiro.test(dt$hf_ng) #hypothesis of normality is rejected
ks.test(dt$hf_ng,"pnorm",mean(dt$hf_ng),sd(dt$hf_ng)) #hypothesis of normality is rejected

# We use log transformed FGCMs
qqnorm(dt$log_hf_ng, datax=TRUE, main="Log-FGCMs")
qqline(dt$log_hf_ng,datax=TRUE)

ks.test(dt$log_hf_ng,"pnorm",mean(dt$log_hf_ng),sd(dt$log_hf_ng)) #hypothesis of normality is tolerated

dt_mod <- dt

## Vizualisation of raw data ----

# New facet label names 
ses.labs <- c("Session 1", "Session 2")
names(ses.labs) <- c("1", "2")
reg.labs <- c("Northern sites", "Southern sites")
names(reg.labs) <- c("north", "south")

counts <- dplyr::count(dt_mod, campaign, region, name = "n")
counts # produces columns: campaign, region, n


# Plot fgcms by site (Supplementary figure S2):

label_y_below_zero <- -10

# 1) per-panel range (campaign x region)
range_df <- dt_mod %>%
  group_by(campaign, region) %>%
  dplyr::summarize(
    min_hf_ng = min(hf_ng, na.rm = TRUE),
    max_hf_ng = max(hf_ng, na.rm = TRUE),
    range_hf_ng = max_hf_ng - min_hf_ng,
    .groups = "drop"
  )

# 2) per-site median and counts (explicitly using dplyr to avoid plyr::summarize)
label_df <- dt_mod %>%
  group_by(campaign, region, site) %>%
  dplyr::summarise(
    median_hf_ng = median(hf_ng, na.rm = TRUE),
    n = sum(!is.na(hf_ng)),
    upper_whisker = ifelse(all(is.na(hf_ng)), NA_real_, boxplot.stats(hf_ng[!is.na(hf_ng)])$stats[5]),
    min_hf_ng = min(hf_ng, na.rm = TRUE),
    max_hf_ng = max(hf_ng, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("n=", n),
    y = label_y_below_zero
  )

# 3) plot using geom_text with label_df
dt_mod %>%
  ggplot(aes(x = reorder_within(site, hf_ng, region, median), y = hf_ng)) +
  geom_boxplot() +
  xlab("Sites") +
  ylab("FGCM concentration (ng/g)") +
  geom_text(
    data = label_df,
    aes(x = reorder_within(site, median_hf_ng, region, FUN = median), y = y, label = label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  geom_text(
    data = counts,
    aes(x = -Inf, y = Inf, label = paste0("N=", n)),
    inherit.aes = FALSE,
    hjust = -0.15,  # tweak: <- moves left/right
    vjust = 1.50,   # tweak: >1 moves down from the top
    size = 3.5
  ) +
  facet_grid(campaign ~ region,
             labeller = labeller(campaign = ses.labs, region = reg.labs),
             scales = "free") +
  scale_x_reordered() +
  coord_cartesian(ylim = c(label_y_below_zero, NA)) +   # make sure the negative y is visible
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

# Plot parasite load by site (epg; Supplementary figure S3):

label_y_below_zero <- -200

# 1) per-panel range (campaign x region)
range_df <- dt_mod %>%
  group_by(campaign, region) %>%
  dplyr::summarize(
    min_epg = min(epg, na.rm = TRUE),
    max_epg = max(epg, na.rm = TRUE),
    range_epg = max_epg - min_epg,
    .groups = "drop"
  )

# 2) per-site median and counts (explicitly using dplyr to avoid plyr::summarize)
label_df <- dt_mod %>%
  group_by(campaign, region, site) %>%
  dplyr::summarise(
    median_epg = median(epg, na.rm = TRUE),
    n = sum(!is.na(epg)),
    upper_whisker = ifelse(all(is.na(epg)), NA_real_, boxplot.stats(epg[!is.na(epg)])$stats[5]),
    min_epg = min(epg, na.rm = TRUE),
    max_epg = max(epg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("n=", n),
    y = label_y_below_zero
  )

# 3) plot using geom_text with label_df
dt_mod %>%
  ggplot(aes(x = reorder_within(site, epg, region, median), y = epg)) +
  geom_boxplot() +
  xlab("Sites") +
  ylab("Parasite load (epg)") +
  geom_text(
    data = label_df,
    aes(x = reorder_within(site, median_epg, region, FUN = median), y = y, label = label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  geom_text(
    data = counts,
    aes(x = -Inf, y = Inf, label = paste0("N=", n)),
    inherit.aes = FALSE,
    hjust = -0.15,  # tweak: <- moves left/right
    vjust = 1.50,   # tweak: >1 moves down from the top
    size = 3.5
  ) +
  facet_grid(campaign ~ region,
             labeller = labeller(campaign = ses.labs, region = reg.labs),
             scales = "free") +
  scale_x_reordered() +
  coord_cartesian(ylim = c(label_y_below_zero, NA)) +   # make sure the negative y is visible
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.7))

# Number of individuals with epg = 0
sum(dt_mod$epg == 0) # 1 individual with epg = 0
# Prevalence of 218/219 = 99.5%.

# MODEL SELECTION ----
# Full model
mfull <- lmer(log_hf_ng ~  
                age + sex + age:sex +
                density_s + 
                log_area_s +
                ndvi_1m_s +
                log_area_s:density_s +
                ndvi_1m_s:density_s +
                self_drive +
                authority +
                epg_s +
                (1|site:session) ,
              data = dt_mod,
              REML = FALSE)

summary(mfull)
check_singularity(mfull) #FALSE
# No singularity issues.

## Model assumptions ----
resid_panel(mfull, smoother = TRUE, qqbands = TRUE) #OK (Supplentary figure S1)
check_autocorrelation(mfull) #autocorr detected

# Full model without interaction terms to calculate Variation Inflation Factor (VIF)
mfull_fe <- lmer(log_hf_ng ~ epg_s + 
                   age + sex +
                   density_s + 
                   ndvi_1m_s +
                   log_area_s +
                   self_drive +
                   authority +
                   (1|site:session),
                 data = dt_mod)

vif <- check_collinearity(mfull_fe) #all VIF<3
vif # Supplementary Table S1.
# No collinearity issues.

## Bakcward elimination based on AICc ----
step <- lmerTest::step(mfull, reduce.random=FALSE)
get_model(step)
print(step) # Supplentary Table S2

# Final model after stepwise selection
modfound <- lmer(log_hf_ng ~ 
                   density_s + 
                   ndvi_1m_s + 
                   log_area_s + 
                   self_drive + 
                   density_s:ndvi_1m_s + 
                   (1|site:session),
                 data = dt_mod,
                 REML=TRUE)

summary(modfound)
resid_panel(modfound, smoother = TRUE, qqbands = TRUE)
check_autocorrelation(modfound) #OK: Residuals appear to be independent and not autocorrelated (p = 0.086)
check_singularity(modfound) #FALSE
# No singularity issues.

## Model selection table (Supplementary Table S3)
mnull <- lmer(log_hf_ng ~ 1 + (1|site:session),
               data = dt_mod,
               REML =FALSE)

mfull <- lmer(log_hf_ng ~  
                age + sex + age:sex +
                density_s + 
                log_area_s +
                ndvi_1m_s +
                log_area_s:density_s +
                ndvi_1m_s:density_s +
                self_drive +
                authority +
                epg_s +
                (1|site:session) ,
              data = dt_mod,
              REML = FALSE)

modfound <- lmer(log_hf_ng ~ 
                   density_s + 
                   ndvi_1m_s + 
                   log_area_s + 
                   self_drive + 
                   density_s:ndvi_1m_s + 
                   (1|site:session),
                 data = dt_mod,
                 REML=FALSE)

cand.set <- list(mfull, modfound, mnull)
names <- c("Full model", "Minimal model", "Null model")
aictab(cand.set, modnames = names)

# MODEL OUTPUTS ----

## Model summary table (Table 1) ----
mfull <- lmer(log_hf_ng ~  
                age + sex + age:sex +
                density_s + 
                log_area_s +
                ndvi_1m_s +
                log_area_s:density_s +
                ndvi_1m_s:density_s +
                self_drive +
                authority +
                epg_s +
                (1|site:session) ,
              data = dt_mod,
              REML = TRUE)

t1 <- modfound %>%
  tbl_regression(exponentiate = F,
                 label = list(
                   density_s ~ "Density",
                   ndvi_1m_s ~ "NDVI",
                   log_area_s ~ "Area size",
                   self_drive ~ "Self-drive"
                 )) %>%
  bold_labels()%>%
  bold_p() %>%
  modify_header(
    label = "**Variables**")

t2 <- mfull %>%
  tbl_regression(exponentiate = F,
                 label = list(
                   epg_s ~ "EPG",
                   age ~ "Age",
                   sex ~ "Sex",
                   density_s ~ "Density",
                   ndvi_1m_s ~ "NDVI",
                   log_area_s ~ "Area size",
                   self_drive ~ "Self-drive",
                   authority ~ "Authority"
                 )) %>%
  bold_labels()%>%
  bold_p() %>%
  modify_header(
    label = "**Variables**")

tbl_merge_ex1 <-
  tbl_merge(
    tbls = list(t2, t1),
    tab_spanner = c("**Full model**", "**Minimal model**")
  ) 
tbl_merge_ex1

# Export the table
ft <- as_flex_table(tbl_merge_ex1)
save_as_docx(ft, path = "model_summary_table.docx")

## Plot effects (Figure 1) ----

# Model found with unscaled:
modfound_ns <- lmer(log_hf_ng ~ density + 
                      ndvi_1m + 
                      log_area + 
                      self_drive + 
                      density:ndvi_1m +
                      (1|site:session),
                    data = dt_mod)

# Interaction plot density x ndvi_1m
pndvi3 <- interact_plot(modfound_ns, pred = ndvi_1m,
                        modx = density,
                        modx.values = c(2, 0.2), 
                        plot.points = F,
                        partial.residuals = T,
                        jitter = 0.005,
                        interval = T,
                        x.label = "NDVI (1 month mean)",
                        y.label = "FGCMs (log scale)",
                        legend.main = "Density of elephants (ind./km²)",
                        line.thickness = 1,
                        colors = c("#0072B2", "#D55E00"),
                        point.alpha = 0.7)+
  scale_linetype_manual(values = c("solid", "solid"), guide = "none")+
  theme(legend.position = c(0.81, 0.11),
        legend.background = element_rect(
          fill = "white",     # background color of the box
          colour = "darkgrey",   # border color
          size = 0.6          # border thickness
        ),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.spacing = unit(0.1, "lines"),
        legend.spacing.y = unit(0.05, "lines"),
        axis.text = element_text(size=10),
        legend.margin = margin(0.5,0.5,0.5,0.5, "mm"),
        panel.background = element_rect(colour = "darkgrey",
                                        linetype = "solid")) +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

pndvi3

# Area size effect plot
parea3 <- effect_plot(modfound_ns, pred = log_area,
                      plot.points = F,
                      partial.residuals = T,
                      jitter = 0.005,
                      interval = TRUE,
                      int.type = "confidence",
                      x.label = "Area size (log scale)",      
                      y.label = "FGCMs (log scale)",
                      line.thickness = 1,
                      colors = "gray8",
                      point.alpha = 0.6) +
  theme(panel.background = element_rect(colour = "darkgrey",
                                        linetype = "solid"),
        axis.text = element_text(size=10))
parea3

# Self-drive effect plot
pdrive3 <- effect_plot(modfound_ns, pred = self_drive,
                       plot.points = F,
                       partial.residuals = T,
                       #interval = TRUE,
                       cat.interval.geom = "linerange",
                       x.label = "Self drive",      
                       y.label = "FGCMs (log scale)",
                       line.thickness = 3,
                       line.colors = c("#009E73", "#DF536B"),
                       point.alpha = 0.4) +
  theme(panel.background = element_rect(colour = "darkgrey",
                                        linetype = "solid"),
        axis.text = element_text(size=10))
pdrive3

# Parasite load null effect plot
pepg <- dt_mod %>%
  ggplot(aes(x = epg, y = log_hf_ng)) +
  geom_point(color = "gray8",
             alpha = 0.6,
             size = 1.5) +
  xlab("Parasite load (epg)") +
  ylab("FGCMs (log scale)") +
  theme_nice() +
  theme(
    panel.background = element_rect(colour = "darkgrey", linetype = "solid"),
    axis.text = element_text(size = 10))
pepg

# Figure 1 for Manuscript:
fig4 <- ggarrange(ggarrange(parea3, pdrive3,
                            ncol = 2,
                            labels = c("A", "B"),
                            font.label=list(color="black",size=16),
                            widths = c(1.4,1)),
                  ggarrange(pndvi3, pepg,
                            ncol = 2,
                            labels = c("C","D"),
                            font.label=list(color="black",size=16),
                            widths = c(1.4,1)),
                  nrow = 2, 
                  heights = c(1,1)) 

fig4

