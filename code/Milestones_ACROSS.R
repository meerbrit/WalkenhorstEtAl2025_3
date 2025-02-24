### Examine the transition ages for FLUT study subjects (F1 & F2) with their sample means ####
### BWalkenhorst, 2024 #######################################################################

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(readxl) 
library(writexl)
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

GEN_data <- readRDS('../data/GEN_milestone_data.rds')

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

################################################################################
################################################################################
### RESULTS ####
################################################################################
################################################################################
# milestone_data <- readRDS('../data/GEN_milestone_data.rds')
# 
# SEMI_data <- milestone_data %>%
#   select(GEN, MAT_STAT, ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, SEMI_AGE) %>%
#   drop_na(SEMI_AGE) %>%
#   rename(Milestone_Age = SEMI_AGE) %>%
#   mutate(Milestone_Type = "SEMI")
# # apply weights because of gen ratio mismatch
# gen_counts <- table(SEMI_data$GEN)
# # Compute weights so that each generation contributes equally
# SEMI_data$gen_weights <- ifelse(SEMI_data$GEN == 'F1', gen_counts[2] / gen_counts[1], 1)
# saveRDS(SEMI_data, 'SEMI_data.rds')
# 
# PEAK_data <- milestone_data %>%
#   select(GEN, MAT_STAT, ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, PEAK_DIG_AGE) %>%
#   drop_na(PEAK_DIG_AGE) %>%
#   rename(Milestone_Age = PEAK_DIG_AGE) %>%
#   mutate(Milestone_Type = "PEAK")
# gen_counts <- table(PEAK_data$GEN)
# PEAK_data$gen_weights <- ifelse(PEAK_data$GEN == 'F1', gen_counts[2] / gen_counts[1], 1)
# saveRDS(PEAK_data, 'PEAK_data.rds')
# 
# FULL_data <- milestone_data %>%
#   select(GEN, MAT_STAT, ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, FULL_AGE) %>%
#   drop_na(FULL_AGE) %>%
#   rename(Milestone_Age = FULL_AGE) %>%
#   mutate(Milestone_Type = "FULL")
# gen_counts <- table(FULL_data$GEN)
# FULL_data$gen_weights <- ifelse(FULL_data$GEN == 'F1', gen_counts[2] / gen_counts[1], 1)
# saveRDS(FULL_data, 'FULL_data.rds')

#### ANOVAs ####
 # REP -> DIG: SEMI ####
SEMI_data <- readRDS('..data/SEMI_data.rds')

mean(SEMI_data$Milestone_Age) #78.09201
sd(SEMI_data$Milestone_Age)# 16.31937

priors_SEMI <- c(
  # General prior for all slopes
  set_prior("normal(0, 15)", class = "b"), # wide range for Rx effects

  # Prior for the intercept
  set_prior("normal(80, 15)", class = "Intercept")
)

SEMI_anova <- brm(formula = Milestone_Age | weights(gen_weights) ~ 
                    (TREATMENT * GEN * MAT_STAT) + (MAT_STAT * GEN * SEX) + (1|ID),
                  data = SEMI_data,
                  family = student(link='identity'), 
                  chains = 4, iter = 3000, warmup = 750, seed = 42234223, control = list(max_treedepth = 20),
                  cores=4, backend = 'cmdstanr', init= 'random',
                  prior = priors_SEMI, threads = threading(4),
                  file="SEMI_anova_ACROSS"
)

# Results: SEMI ####
SEMI_anova <- readRDS("../models/SEMI_anova_ACROSS.rds")

summary(SEMI_anova)


plot(SEMI_anova)
pp_check(SEMI_anova, ndraws=100)

loo_R2(SEMI_anova, moment_match = T)

bayes_R2(SEMI_anova)

describe_posterior(
  SEMI_anova,
  effects = "all", 
  component = "all",
  rope_range = rope_range(SEMI_anova),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

### TG ####
(treat_gen<- emmeans(SEMI_anova, pairwise ~ TREATMENT:GEN, nesting = c("MAT_STAT %in% GEN")))

pd(treat_gen)

p_significance(treat_gen, threshold = rope_range(SEMI_anova))

# Plot TG
treat_gen <- emmeans(SEMI_anova, ~ TREATMENT:GEN, nesting = c("MAT_STAT %in% GEN"))
emm_data <- as.data.frame(treat_gen)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))
emm_data$GEN <- as.factor(emm_data$GEN)
levels(emm_data$GEN) <- c('F1', 'F2')
SEMI_data$GEN <- as.factor(SEMI_data$GEN)
levels(SEMI_data$GEN) <- c('F1', 'F2')

milestones_org <- SEMI_data %>%
  group_by(ID, TREATMENT, SEX, GEN, MAT_STAT) %>%
  summarise(TRANS_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")
  
# 700* 500: SEMI_TG
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = GEN)) +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = TRANS_age, shape = GEN, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "F0 status", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  scale_shape_discrete(name = 'Generation', labels=c('F1', 'F2'))+
  labs(x = "F0 status", y = "Age (days)\n") +
  theme_clean()

rm(treat_gen)

### F2 only: TM ####
(F2_only <- emmeans(SEMI_anova, pairwise ~ TREATMENT:MAT_STAT, at = list(GEN = 'F2', MAT_STAT=c('DOM', 'SUB')), nesting = c("MAT_STAT %in% GEN")))

pd(F2_only)

p_significance(F2_only, threshold = rope_range(SEMI_anova))

# Plot TM
treat_mat <- emmeans(SEMI_anova, ~ TREATMENT:MAT_STAT, at = list(GEN = 'F2', MAT_STAT=c('DOM', 'SUB')), nesting = c("MAT_STAT %in% GEN"))
emm_data <- as.data.frame(treat_mat)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

# no F1
f2_data <- milestones_org %>%
  filter(GEN == "F2")

# 800* 500: SEMI_F2_TM
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = MAT_STAT)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Grandmaternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Maternal status', labels = c('dominant', 'subordinate'))+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Grandmaternal\nstatus", y = "Age (days)\n") +
  geom_point(data = f2_data, aes(x = TREATMENT, y = TRANS_age, color=TREATMENT, shape = MAT_STAT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()

rm(treat_mat, f2_data, emm_data)

### PLOTS: SEMI ####
# coefficient plots #### 
posterior_desc <- describe_posterior(
  SEMI_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(SEMI_anova),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1, 25),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GENF2', 'F2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATSUB', 'F2 Maternal status = SUB')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATDOM', 'F2 Maternal status = DOM')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATF1', 'F1')

custom_order <- c("DT:F2:F2 Maternal status = SUB","SC:F2:F2 Maternal status = SUB","F2:F2 Maternal status = SUB",
                  "DT:F2:F2 Maternal status = DOM","SC:F2:F2 Maternal status = DOM","F2:F2 Maternal status = DOM",
                  "F2:F2 Maternal status = SUB:Male", "F2:F2 Maternal status = DOM:Male",
                  "F2 Maternal status = DOM:Male", "F2 Maternal status = SUB:Male",
                  "DT:F2 Maternal status = SUB", "SC:F2 Maternal status = SUB", "F2 Maternal status = SUB", 
                  "DT:F2 Maternal status = DOM", "SC:F2 Maternal status = DOM", "F2 Maternal status = DOM" ,
                  'F2:Male',
                  'DT:Male','SC:Male',"Male", 
                  'DT:F2','SC:F2','F2',
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_SEMI 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1,3), name = "F0 status", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

################################################################################
## DIG -> CC ####
FULL_data <- readRDS('../data/FULL_data.rds')

# mean(FULL_data$Milestone_Age) # 118
# sd(FULL_data$Milestone_Age)# 18

priors_FULL <- c(
  # General prior for all slopes
  set_prior("normal(0, 15)", class = "b"), # wide range for Rx effects

  # Prior for the intercept
  set_prior("normal(120, 20)", class = "Intercept")
)

FULL_anova<- brm(formula = Milestone_Age| weights(gen_weights) ~
                   (TREATMENT * GEN * MAT_STAT) + (MAT_STAT * GEN * SEX) + (1|ID),
                 data = FULL_data,
                 family = student(link='identity'), 
                 chains = 4, iter = 3000, warmup = 750, seed = 42234223, control = list(max_treedepth = 20),
                 cores=4, backend = 'cmdstanr', init= 'random',
                 prior = priors_FULL,
                 threads = threading(4),
                 file="FULL_anova_ACROSS")


rm(priors_FULL)

### Results: FULL ####
FULL_anova <- readRDS('../models/FULL_anova_ACROSS.rds')

plot(FULL_anova)
pp_check(FULL_anova, ndraws=100)

loo_R2(FULL_anova, moment_match=T)

bayes_R2(FULL_anova)

summary(FULL_anova)


describe_posterior(
  FULL_anova,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = rope_range(FULL_anova),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)


### EMMs ####
### TG ####
(treat_gen <- emmeans(FULL_anova, pairwise ~ TREATMENT:GEN, nesting = c("MAT_STAT %in% GEN")))

pd(treat_gen)

p_significance(treat_gen, threshold = rope_range(FULL_anova))

# Plot TG
treat_gen <- emmeans(FULL_anova, ~ TREATMENT:GEN, nesting = c("MAT_STAT %in% GEN"))
emm_data <- as.data.frame(treat_gen)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))
emm_data$GEN <- as.factor(emm_data$GEN)
levels(emm_data$GEN) <- c('F1', 'F2')
FULL_data$GEN <- as.factor(FULL_data$GEN)
levels(FULL_data$GEN) <- c('F1', 'F2')

milestones_org <- FULL_data %>%
  group_by(ID, TREATMENT, SEX, GEN, MAT_STAT) %>%
  summarise(TRANS_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")

# 700* 500: FULL_TG
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = GEN)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "F0 status", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Generation')+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "F0 status", y = "Age (days)\n") +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = TRANS_age, shape = GEN, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()

rm(treat_gen)

### F2 TM ### 
(F2_only <- emmeans(FULL_anova, pairwise ~ MAT_STAT:TREATMENT, at = list(GEN = 'F2', MAT_STAT=c('DOM', 'SUB')), nesting = c("MAT_STAT %in% GEN")))

pd(F2_only)

p_significance(F2_only, threshold = rope_range(FULL_anova))

# Plot TM
treat_mat <- emmeans(FULL_anova, ~ MAT_STAT:TREATMENT, at = list(GEN = 'F2', MAT_STAT=c('DOM', 'SUB')), nesting = c("MAT_STAT %in% GEN"))
emm_data <- as.data.frame(treat_mat)
emm_data$TRANS_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

# no F1
f2_data <- milestones_org %>%
  filter(GEN == 'F2')

# 700* 500: FULL_TM
ggplot(emm_data, aes(x = TREATMENT, y = TRANS_age, color = TREATMENT, shape = MAT_STAT)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "F0 status", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Maternal status', labels=c('dominant', 'subordinate'))+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "F0 status", y = "Age (days)\n") +
  geom_point(data = f2_data, aes(x = TREATMENT, y = TRANS_age, shape = MAT_STAT, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()


rm(treat_mat, emm_data, f2_data)

#### PLOTS: FULL ####
# coefficient plots #### 
posterior_desc <- describe_posterior(
  FULL_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(FULL_anova),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1, 25),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GENF2', 'F2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATSUB', 'F2 Maternal status = SUB')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATDOM', 'F2 Maternal status = DOM')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATF1', 'F1')

custom_order <- c("DT:F2:F2 Maternal status = SUB","SC:F2:F2 Maternal status = SUB","F2:F2 Maternal status = SUB",
                  "DT:F2:F2 Maternal status = DOM","SC:F2:F2 Maternal status = DOM","F2:F2 Maternal status = DOM",
                  "F2:F2 Maternal status = SUB:Male", "F2:F2 Maternal status = DOM:Male",
                  "F2 Maternal status = DOM:Male", "F2 Maternal status = SUB:Male",
                  "DT:F2 Maternal status = SUB", "SC:F2 Maternal status = SUB", "F2 Maternal status = SUB", 
                  "DT:F2 Maternal status = DOM", "SC:F2 Maternal status = DOM", "F2 Maternal status = DOM" ,
                  'F2:Male',
                  'DT:Male','SC:Male',"Male", 
                  'DT:F2','SC:F2','F2',
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_FULL 800*900
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "F0 status", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

rm(SEMI_data, FULL_data, SEMI_anova, FULL_anova)

### Peak DIG
DIG_PEAK_data <- readRDS('../data/PEAK_data_PRIOR.rds')


# mean(DIG_PEAK_data$Milestone_Age) #100
# sd(DIG_PEAK_data$Milestone_Age)# 19.92

priors_PEAK <- c(
  # General prior for all slopes
  set_prior("normal(0, 15)", class = "b"), # wide range for Rx effects

  # Prior for the intercept
  set_prior("normal(100, 20)", class = "Intercept")
)

DIG_PEAK_anova<- brm(formula = Milestone_Age | weights(gen_weights) ~
                       (TREATMENT * GEN * MAT_STAT) + (MAT_STAT * GEN * SEX) + (1|ID),
                     data = DIG_PEAK_data,
                     family = student(link = 'identity'),
                     chains = 4, iter = 3000, warmup = 750, seed = 23425235, control = list(max_treedepth = 20),
                     cores=4, backend = 'cmdstanr', init='random',
                     prior = priors_PEAK,
                     threads = threading(4),
                     file="DIG_PEAK_anova_ACROSS")
rm(priors_PEAK)

# Results_ PEAK ###
DIG_PEAK_anova <- readRDS("../models/DIG_PEAK_anova_ACROSS_PRIOR.rds")

plot(DIG_PEAK_anova)
pp_check(DIG_PEAK_anova, ndraws=100)

loo_R2(DIG_PEAK_anova, moment_match=T)

bayes_R2(DIG_PEAK_anova)

summary(DIG_PEAK_anova)

describe_posterior(
  DIG_PEAK_anova,
  effects = "all", 
  component = "all",
  rope_range = rope_range(DIG_PEAK_anova),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

#### EMMs ####
### TG ###
(treat_gen <- emmeans(DIG_PEAK_anova, pairwise ~ TREATMENT:GEN, nesting = c("MAT_STAT %in% GEN")))

pd(treat_gen)

p_significance(treat_gen, threshold = rope_range(DIG_PEAK_anova))

# Plot TG
treat_gen <- emmeans(DIG_PEAK_anova, ~ TREATMENT:GEN, nesting = c("MAT_STAT %in% GEN"))
emm_data <- as.data.frame(treat_gen)
emm_data$PEAK_age <- emm_data$emmean
emm_data$HPD_low <- emm_data$lower.HPD
emm_data$HPD_high <- emm_data$upper.HPD
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))
emm_data$GEN <- as.factor(emm_data$GEN)
levels(emm_data$GEN) <- c('F1', 'F2')
levels(DIG_PEAK_data$GEN) <- c('F1', 'F2')

milestones_org <- DIG_PEAK_data %>%
  group_by(ID, TREATMENT, SEX, GEN, MAT_STAT) %>%
  summarise(TRANS_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")

# 800* 500: PEAK_DIG_TG
ggplot(emm_data, aes(x = TREATMENT, y = PEAK_age, color = TREATMENT, shape = GEN)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "F0 status", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Generation', labels = c('F1', 'F2'))+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "F0 status", y = "Age (days)\n") +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = TRANS_age, shape = GEN, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()

rm(treat_gen)

### F2 only: TM ####
(F2_only <- emmeans(DIG_PEAK_anova, pairwise ~ MAT_STAT:TREATMENT, at = list(GEN = 'F2', MAT_STAT=c('DOM', 'SUB')), nesting = c("MAT_STAT %in% GEN")))

pd(F2_only)

p_significance(F2_only, threshold = rope_range(DIG_PEAK_anova))

# Plot TM
treat_mat <- emmeans(DIG_PEAK_anova, ~ MAT_STAT:TREATMENT, at = list(GEN = 'F2', MAT_STAT=c('DOM', 'SUB')), nesting = c("MAT_STAT %in% GEN"))
emm_data <- as.data.frame(treat_mat)
emm_data$PEAK_age <- (emm_data$emmean)
emm_data$HPD_low <- (emm_data$lower.HPD)
emm_data$HPD_high <- (emm_data$upper.HPD)
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

# no F1
f2_data <- milestones_org %>%
  filter(GEN == 'F2')

# 800* 500: PEAK_DIG_TM
ggplot(emm_data, aes(x = TREATMENT, y = PEAK_age, color = TREATMENT, shape = MAT_STAT)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "F0 status", labels = c('DC', 'SC', 'DT')) +
  scale_shape_discrete(name='Maternal status', labels = c('dominant', 'subordinate'))+
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "F0 status", y = "Age (days)\n") +
  geom_point(data = f2_data, aes(x = TREATMENT, y = TRANS_age, shape = MAT_STAT, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3)  +  
  theme_clean()

rm(treat_mat_sex, treat_mat, emm_data, f2_data)

# coefficient plot ####
posterior_desc <- describe_posterior(
  DIG_PEAK_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(DIG_PEAK_anova),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1, 25),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'GENF2', 'F2')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATSUB', 'F2 Maternal status = SUB')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATDOM', 'F2 Maternal status = DOM')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'MAT_STATF1', 'F1')

custom_order <- c("DT:F2:F2 Maternal status = SUB","SC:F2:F2 Maternal status = SUB","F2:F2 Maternal status = SUB",
                  "DT:F2:F2 Maternal status = DOM","SC:F2:F2 Maternal status = DOM","F2:F2 Maternal status = DOM",
                  "F2:F2 Maternal status = SUB:Male", "F2:F2 Maternal status = DOM:Male",
                  "F2 Maternal status = DOM:Male", "F2 Maternal status = SUB:Male",
                  "DT:F2 Maternal status = SUB", "SC:F2 Maternal status = SUB", "F2 Maternal status = SUB", 
                  "DT:F2 Maternal status = DOM", "SC:F2 Maternal status = DOM", "F2 Maternal status = DOM" ,
                  'F2:Male',
                  'DT:Male','SC:Male',"Male", 
                  'DT:F2','SC:F2','F2',
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_PEAK_DIG 700*800
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1,3), name = "F0 status", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

rm(GEN_data, full_data, posterior_desc, GEN_data, PROP_priors,
   DIG_data, DIG_PEAK_anova, emms_peak_FULL,  peak_data, MIN_AGE, MAX_AGE, theme_clean, custom_order,
   priors_PEAK, DIG_PEAK_data, environment_vars, F2_only, orthogonalised_vars)

### cleanup ###
rm(GEN_data, B_prop, priors, transition_data, GEN_data, environment_vars,
   F2_only, SEMI_anova, FULL_anova, priors_FULL, priors_SEMI, SEMI_data,
   FULL_data, posterior_desc, orthogonalised_vars, custom_order, theme_clean)
