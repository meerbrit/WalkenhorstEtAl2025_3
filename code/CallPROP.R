#### Analysis script for 2nd GEN Flutamide study: CALL PROPORTIONS #####################
# Bayesian Multilevel models: 
# Treatment (T), age (A), sex (S), body mass offset (W), 
# also checking for (normalised) competition score (C = H/P ratio & pups & adults),
# group size and cumulative rainfall
#
# BWalkenhorst, 2024 

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(ggplot2) 
library(ggpubr) 
library(car)
library(tidyverse)
library(tidybayes) 
library(brms)
library(rstan)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # 
library(extrafont)# use font_import() on first use

set.seed(23)

PROP_data <- readRDS('../data/PROP_data.rds')

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal(base_family='Calibri') +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(1.75)),
      axis.text = element_text(face = "bold", size= rel(1.25)),
      strip.text = element_text(face = "bold", size = rel(1.5), color='white'),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold", size = rel(1.25)),
      legend.text = element_text(face = 'italic', size = rel(1)))
}

get_age_vars <- function(){
  age_vars <- c((30- mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE), 
                (75- mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE),
                (120- mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE))
  return(age_vars)
}

get_org_value_for_z <- function(z_value, column){
  return ((z_value*sd(column)) + mean(column))   
}

get_z_value_for_org <- function(org_value, column){
  return ((org_value - mean(column))/sd(column))   
}

########################## BAYES MULTIVARIATE MODEL########################################
# weakly informative priors based on the first gen model
priors <- c(
  # Intercept: Still broadly informed but allows variability
  set_prior("normal(0, 2)", class = "Intercept", resp = "SumBEG"),
  set_prior("normal(0, 2)", class = "Intercept", resp = "SumDIG"),
  set_prior("normal(0, 3)", class = "Intercept", resp = "SumCC"),  # More variation in CC
  
  # Fixed effect coefficients: Weakly informative, keeping response-specific scaling
  set_prior("normal(0, 1.5)", class = "b", resp = "SumBEG"),
  set_prior("normal(0, 1)", class = "b", resp = "SumDIG"),
  set_prior("normal(0, 1.5)", class = "b", resp = "SumCC"),
  
  # Random effect (only ID-level variation)
  set_prior("normal(0.5, 0.3)", class = "sd", group = "ID", resp = "SumBEG"),
  set_prior("normal(0.6, 0.3)", class = "sd", group = "ID", resp = "SumDIG"),
  set_prior("normal(0.7, 0.4)", class = "sd", group = "ID", resp = "SumCC"),  # More variability in CC
  
  # Overdispersion priors (phi) - adjusted for scale
  set_prior("gamma(3, 1)", class = "phi", resp = "SumBEG"),
  set_prior("gamma(4, 1)", class = "phi", resp = "SumDIG"),  
  set_prior("gamma(5, 1.5)", class = "phi", resp = "SumCC"),  # More variation in CC
  
  # Zero-inflation priors (zi) - slightly informed
  set_prior("beta(2, 2)", class = "zi", resp = "SumBEG"),
  set_prior("beta(2, 2)", class = "zi", resp = "SumDIG"),
  set_prior("beta(3, 2)", class = "zi", resp = "SumCC")  # More probability around 0.1-0.2
)


bf_REP <- bf(Sum_BEG | trials(Total_calls) ~ TREATMENT + AGE_z + MaternalRank + SEX +
               TREATMENT:AGE_z + TREATMENT:MaternalRank + AGE_z:MaternalRank + 
               TREATMENT:AGE_z:MaternalRank +
               MaternalRank:AGE_z + MaternalRank:SEX + AGE_z:SEX + MaternalRank:AGE_z:SEX +
               WEIGHT_z + COMP_NORM_z + GS_z +  RAIN_z + (1|ID))

bf_DIG <- bf(Sum_DIG | trials(Total_calls) ~ TREATMENT + AGE_z + I(AGE_z^2) + MaternalRank + SEX +
               TREATMENT:AGE_z + TREATMENT:I(AGE_z^2) + TREATMENT:MaternalRank + AGE_z:MaternalRank + I(AGE_z^2):MaternalRank +
               TREATMENT:AGE_z:MaternalRank + TREATMENT:I(AGE_z^2):MaternalRank +
               MaternalRank:AGE_z + MaternalRank:I(AGE_z^2) + MaternalRank:SEX + AGE_z:SEX + I(AGE_z^2):SEX +
               MaternalRank:AGE_z:SEX + MaternalRank:I(AGE_z^2):SEX +
               WEIGHT_z + COMP_NORM_z + GS_z +  RAIN_z + (1|ID))

bf_CC <- bf(Sum_CC | trials(Total_calls) ~ TREATMENT + AGE_z + MaternalRank + SEX +
              TREATMENT:AGE_z + TREATMENT:MaternalRank + AGE_z:MaternalRank + 
              TREATMENT:AGE_z:MaternalRank +
              MaternalRank:AGE_z + MaternalRank:SEX + AGE_z:SEX + MaternalRank:AGE_z:SEX +
              WEIGHT_z + COMP_NORM_z + GS_z +  RAIN_z + (1|ID))


multivar_formula <- mvbrmsformula(bf_REP, bf_DIG, bf_CC)

B_prop <- brms::brm(formula = multivar_formula,
                    data = PROP_data, family = zero_inflated_beta_binomial(link='logit'),
                    chains = 4, iter = 10000, warmup = 2000, seed = 23, control = list(max_treedepth = 20, adapt_delta=0.99),
                    save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                    threads = threading(4),
                    prior = priors,
                    file="B_prop")


#### RESULTS: MULTIVAR model ####
B_prop <- readRDS("../models/B_prop.rds")

summary(B_prop)

plot(B_prop)
pp_check(B_prop, ndraws = 100, resp='SumBEG')
pp_check(B_prop, ndraws = 100, resp='SumDIG')
pp_check(B_prop, ndraws = 100, resp='SumCC')

posterior <- describe_posterior(
  B_prop$fit,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = c(-0.18, 0.18),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)
(output <- posterior[1:80,])
saveRDS(output, file='posterior_desc_B_prop_prior.rds')

loo_R2(B_prop, moment_match=T)

bayes_R2(B_prop)


### REP: TAM ####
REP_ROPE <- c(-0.18, 0.18)
(mat_stat<- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank, var="AGE_z", resp='SumBEG'))

pd(mat_stat)

p_significance(mat_stat, threshold = REP_ROPE)

rm(mat_stat)

### REP at 30 days (intercept) ####
REP_ROPE <- c(-0.18, 0.18)
(REP_int <- emmeans(B_prop,
                           pairwise ~ TREATMENT * MaternalRank | AGE_z,
                           at = list(AGE_z = get_z_value_for_org(30, PROP_data$REC_AGE)),
                           resp='SumBEG') )
                   # ,type='response') )

pd(REP_int)

p_significance(REP_int, threshold = REP_ROPE)

rm(REP_int)

### DIG slopes ####
DIG_ROPE <- c(-0.18, 0.18)
(treat_mat_stat <- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank, var = "AGE_z", max.degree = 2, resp = "SumDIG"))

pd(treat_mat_stat)

p_significance(treat_mat_stat, threshold = DIG_ROPE)

### DIG: aged 75 days ####
DIG_ROPE <- c(-0.18, 0.18)
(treat_mat_stat <- emmeans(B_prop,
                          pairwise ~ TREATMENT * MaternalRank | AGE_z,
                           at = list(AGE_z = get_z_value_for_org(75, PROP_data$REC_AGE)),
                           resp='SumDIG'))
                         # , type='response'))

pd(treat_mat_stat)


p_significance(treat_mat_stat, threshold = DIG_ROPE)


### CC TAM ####
CC_ROPE <- c(-0.18, 0.18)
(treat_mat <- emtrends(B_prop, pairwise ~ TREATMENT:MaternalRank, var="AGE_z", resp='SumCC'))

pd(treat_mat)

p_significance(treat_mat, threshold = CC_ROPE)

### CC at 120 days ####
CC_ROPE <- c(-0.18, 0.18)
(treat_mat_stat <- emmeans(B_prop,
                           pairwise ~ TREATMENT * MaternalRank | AGE_z,
                           at = list(AGE_z = get_z_value_for_org(120, PROP_data$REC_AGE)),
                           resp='SumCC' ) )
                          # , type='response'))

pd(treat_mat_stat)

p_significance(treat_mat_stat, threshold = CC_ROPE)


rm(treat_sex, mat_sex, treat_mat, treat_mat_stat,  CC_ROPE, DIG_ROPE, REP_ROPE)

#### Coefficient plots ####
posterior_desc <- readRDS('posterior_desc_B_prop_prior.rds')

REP_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumBEG"))
DIG_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumDIG"))
CC_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumCC"))

REP_desc <- REP_desc[(2:20),]
DIG_desc <- DIG_desc[(2:28),]
CC_desc <- CC_desc[c(2:20),]

# REP Coeff
# clean up labels:
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'b_SumBEG_', '')
REP_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", REP_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", REP_desc$Parameter), "DT", "DC"))
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTDT', 'DT')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTSC', 'SC')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'GS_z', 'Group size')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'IGroup sizeE2', 'Group size^2')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'AGE_z', 'Age')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'COMP_NORM_z', 'Competition')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'SEXM', 'Male')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'MaternalRankSUB', 'SUB mother')

custom_order <- c('DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'Age:SUB mother:Male',
                  'SUB mother:Male', 
                  'DT:Monthly rainfall','SC:Monthly rainfall','Monthly rainfall', 
                  'DT:Group size^2','SC:Group size^2','Group size^2', 
                  'DT:Group size','SC:Group size','Group size', 
                  'DT:Competition','SC:Competition','Competition', 
                  'DT:Body mass offset','SC:Body mass offset','Body mass offset', 
                  'DT:Male','SC:Male',"Male", 
                  'DT:Age','SC:Age',"Age", 
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT','SC') 

# Update the order of TREATMENT factor levels
REP_desc$TREATMENT <- factor(REP_desc$TREATMENT, levels = c("DC", "SC", "DT"))
REP_desc$Parameter <- factor(REP_desc$Parameter, levels = custom_order)
# Coeff_REP 700 * 800
ggplot(REP_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# DIG coeff
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'b_SumDIG_', '')
DIG_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", DIG_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", DIG_desc$Parameter), "DT", "DC"))
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTDT', 'DT')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTSC', 'SC')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IAGE_zE2', 'Age^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'AGE_z', 'Age')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'GS_z', 'Group size')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IGroup sizeE2', 'Group size^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'COMP_NORM_z', 'Competition')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'SEXM', 'Male')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'MaternalRankSUB', 'SUB mother')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c(
                  'DT:Age^2:SUB mother','SC:Age^2:SUB mother',"Age^2:SUB mother", 
                  'DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:Age^2:Male','SC:Age^2:Male','Age^2:Male',
                  'DT:Age:Male','SC:Age:Male','Age:Male',
                  "Age^2:SUB mother:Male", 
                  "Age:SUB mother:Male" ,
                  'DT:SUB mother:Male','SC:SUB mother:Male',"SUB mother:Male",
                  'Monthly rainfall',
                  'Group size^2', 
                  'Group size', 
                  'Competition', 
                  'Body mass offset', 
                  'DT:Male','SC:Male',"Male", 
                  'DT:Age^2','SC:Age^2',"Age^2", 
                  'DT:Age','SC:Age',"Age",
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT','SC') 

DIG_desc$Parameter <- factor(DIG_desc$Parameter, levels = custom_order)
DIG_desc$TREATMENT <- factor(DIG_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG 700*900
ggplot(DIG_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# CC Coeff
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'b_SumCC_', '')
CC_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", CC_desc$Parameter), "SC",
                            ifelse(grepl("TREATMENTDT", CC_desc$Parameter), "DT", "DC"))
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTDT', 'DT')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTSC', 'SC')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'IGS_zE2', 'Group size^2')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'GS_z', 'Group size')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'AGE_z', 'Age')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'COMP_NORM_z', 'Competition')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'SEXM', 'Male')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'MaternalRankSUB', 'SUB mother')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'RAIN_z', 'Monthly rainfall')

custom_order <- c('DT:Age:SUB mother','SC:Age:SUB mother','Age:SUB mother',
                  'DT:Age:Male','SC:Age:Male','Age:Male', 
                  'Age:SUB mother:Male',
                  'SUB mother:Male', 
                  'DT:Monthly rainfall','SC:Monthly rainfall','Monthly rainfall', 
                  'DT:Group size^2','SC:Group size^2','Group size^2', 
                  'DT:Group size','SC:Group size','Group size', 
                  'DT:Competition','SC:Competition','Competition', 
                  'DT:Body mass offset','SC:Body mass offset','Body mass offset', 
                  'DT:Male','SC:Male',"Male", 
                  'DT:Age','SC:Age',"Age", 
                  'DT:SUB mother','SC:SUB mother','SUB mother', 
                  'DT','SC') 

CC_desc$Parameter <- factor(CC_desc$Parameter, levels = custom_order)
CC_desc$TREATMENT <- factor(CC_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_CC 700*800
ggplot(CC_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

rm(REP_desc, DIG_desc, CC_desc, posterior_desc, custom_order)

#### MODEL ONTOGENY PLOTS ####
B_prop <- readRDS('B_prop.rds')

# get all needed values
sd_age <- sd(PROP_data$REC_AGE) 
mean_age <- mean(PROP_data$REC_AGE) 

range(PROP_data$REC_AGE)

rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE)
# 1 day steps not needed here 
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

#predictions based on mean values
PROP_pred <- B_prop %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(PROP_data$TREATMENT),
                                    SEX = levels(PROP_data$SEX),
                                    MaternalRank = levels(PROP_data$MaternalRank),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(PROP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(PROP_data$COMP_NORM_z),
                                    GS_z = mean(PROP_data$GS_z),
                                    RAIN_z = mean(PROP_data$RAIN_z),
                                    Total_calls=1), 
              re_formula = NA,  robust = T)


#unscale AGE_z values:
PROP_pred$REC_AGE <- PROP_pred$AGE_z * sd_age + mean_age
# ensure right format
PROP_pred$Call_prop <- PROP_pred$.epred
PROP_pred$Call_type <- as.factor(PROP_pred$.category)
PROP_pred$SEX <- as.factor(PROP_pred$SEX)
PROP_pred$MaternalRank <- as.factor(PROP_pred$MaternalRank)
PROP_pred$TREATMENT <- factor(PROP_pred$TREATMENT, levels = c("DC", "SC", "DT"))

# plot raw proportions
ggplot(PROP_pred, aes(x = REC_AGE, y = Call_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +  
  stat_lineribbon(.width = .95) + # shows uncertainty
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), name = 'Call type', labels = c('REP', 'DIG', 'CC'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# just REP ###
REP_pred <- subset(PROP_pred, Call_type == 'SumBEG')
REP_pred$REP_prop <- REP_pred$Call_prop

# TAM
ggplot(REP_pred, aes(x = REC_AGE, y = REP_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_wrap(~MaternalRank)
#MAS
ggplot(REP_pred, aes(x = REC_AGE, y = REP_prop, color = MaternalRank, fill = MaternalRank)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  labs(x = "Age (days)", y = "Repeat call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

rm(REP_pred)

# # just DIG
DIG_pred <- subset(PROP_pred, Call_type == 'SumDIG')
DIG_pred$DIG_prop <- DIG_pred$Call_prop

# MAS 800*500
ggplot(DIG_pred, aes(x = REC_AGE, y = DIG_prop, color = MaternalRank, fill = MaternalRank)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  labs(x = "Age (days)", y = "Digging call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()#+
  facet_wrap(~SEX)

# TAM 800*500
ggplot(DIG_pred, aes(x = REC_AGE, y = DIG_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() #+
  facet_wrap(~MaternalRank)

rm(DIG_pred)

# just CC
CC_pred <- subset(PROP_pred, Call_type == 'SumCC')
CC_pred$CC_prop <- CC_pred$Call_prop

# TAM 800*500
ggplot(CC_pred, aes(x = REC_AGE, y = CC_prop, color = TREATMENT, fill = TREATMENT)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean() +
  facet_wrap(~MaternalRank)

# MAS
ggplot(CC_pred, aes(x = REC_AGE, y = CC_prop, color = MaternalRank, fill = MaternalRank)) +
  stat_lineribbon(.width = .95) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('dominant', 'subordinate')) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n",
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

rm(CC_pred)



