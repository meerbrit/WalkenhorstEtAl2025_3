###### Model weights of 2nd Gen FLUT litters and determine offset to expected body mass#####
#### BWalkenhorst 2024 ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() 

# Load required libraries
library(dplyr)
library(zoo)
library(lubridate)
library(ggplot2)
library(readxl)
library(writexl)
library(brms)
library(rstan)
library(tidyverse)
library(tidybayes)
library(ggokabeito) # colour palette
library(extrafont)# use font_import() on first use
library(bayestestR)

# set seed to duplicate results
set.seed(42)

# half normal prior and weak normal for intercept
priors_halfnormal <- c(set_prior('normal(0,0.5)', class = 'b', lb = 0), set_prior("normal(0,1)", class = "Intercept"))

# reset working directory
final_data <- readRDS('../data/WEIGHT_final_data.rds')

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

#### MODEL ####
WEIGHT_AGE <- brms::brm(formula = scale(Weight) ~ scale(AGE_D) + SEX + scale(Rainfall30D) + (1|ID), 
                                     data = final_data, family = gaussian(link='identity'),
                                     chains = 4, iter = 10000, warmup = 2500, cores = 4, backend = "cmdstanr", 
                                     prior = priors_halfnormal,
                                     control = list(max_treedepth = 15, adapt_delta=0.999), init=0, 
                                     threads = threading(4),
                                     file ="WEIGHT_AGE")

#### Model details ####
WEIGHT_AGE <- readRDS("../models/WEIGHT_AGE.rds")
summary(WEIGHT_AGE)

plot(WEIGHT_AGE)
pp_check(WEIGHT_AGE, ndraws=100)

# get the rope range  = -0.1 * SDy, 0.1 * SDy
# as its scaled, sd = 1 so -0.1 and 0.1 it is!
ropeRange <- c(-0.1* sd(scale(final_data$Weight)), 0.1 * sd(scale(final_data$Weight)))

describe_posterior(
  WEIGHT_AGE,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = ropeRange,  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)


loo_R2(WEIGHT_AGE)

bayes_R2(WEIGHT_AGE)


performance::variance_decomposition(WEIGHT_AGE)

### predict weight for FLUT pups ####
# FLUT data with average weight and rainfall
flut_data <- read_excel(
  "../data/2ND_GEN_Weight_Rain.xlsx", 
)
flut_data <- flut_data %>%
  mutate(TREATMENT = as.factor(TREATMENT), SEX = as.factor(SEX), ID = as.factor(ID))

# get the predictions from the model
#weight_pred <- predict(WEIGHT_AGE, newdata = copy_flut, seed=23, allow_new_levels=T) 

#load predictions if not calculated
weight_pred <- readRDS("FLUT_RAIN_weight.rds")

flut_data$WEIGHT_PRED <- as.numeric(weight_pred[,1])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)
flut_data$WEIGHT_CI.lo <- as.numeric(weight_pred[,3])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)
flut_data$WEIGHT_CI.hi <- as.numeric(weight_pred[,4])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)

############## determine BODY MASS OFFSET ####
flut_data$WEIGHT_DIFF <- flut_data$AvgWeight - flut_data$WEIGHT_PRED 

#use percentage
flut_data$WEIGHT_DIFF_PER <-as.numeric(((flut_data$AvgWeight - flut_data$WEIGHT_PRED)/flut_data$WEIGHT_PRED) *100)

# save the weight offset
write_xlsx(flut_data, "../data/DEV_DATA_2ndGEN.xlsx")

