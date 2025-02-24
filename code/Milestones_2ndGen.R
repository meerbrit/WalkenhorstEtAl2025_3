### Calculate milestone ages for FLUT study subjects (by ID!) and include all predictors ####
### BWalkenhorst, 2024 #######################################################################
#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(ggplot2) 
library(ggpubr) 
library(tidyverse)
library(tidybayes) 
library(brms)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # 

set.seed(23)

PROP_data <- read_rds('../data/PROP_data.rds')
B_prop <- readRDS('../models/B_prop.rds')

MIN_AGE = 1
MAX_AGE = 180
PROP_DIFF = 0.01
NUM_SAMPLES_PER_ID = 100

# ==============================================================================
#  FUNCTION TO EXTRACT POSTERIOR CALL PROPORTIONS
# ==============================================================================
extract_posterior_call_props <- function(row) {
  rec_age_c <- seq(MIN_AGE, MAX_AGE, by = 1)
  age_z_vals <- (rec_age_c - mean(PROP_data$REC_AGE)) / sd(PROP_data$REC_AGE)
  
  # Expand AGE_z while keeping all other predictors constant
  new_data <- expand_grid(
    AGE_z = age_z_vals
  ) %>%
    mutate(
      REC_AGE = rec_age_c,  # Assign corresponding REC_AGE
      TREATMENT = row$TREATMENT,
      SEX = row$SEX,
      MaternalRank = row$MAT_STAT,
      WEIGHT_z = row$WEIGHT_z,
      COMP_NORM_z = row$COMP_z,
      GS_z = row$GS_z,
      RAIN_z = row$RAIN_z,
      MOTHER_ID = row$MOTHER_ID,
      LITTER_CODE = row$LITTER_CODE,
      ID = row$ID,
      Total_calls = 1,
      .row = row_number()  # Assign `.row` index for merging
    )
  
  # Generate posterior predictions
  PROP_pred <- B_prop %>%
    epred_draws(newdata = new_data, re_formula = NULL, allow_new_levels = FALSE, robust = TRUE, ndraws = NULL)  # Keep ndraws draws
  
  # Merge posterior predictions with `new_data`
  PROP_pred_tidy <- PROP_pred %>%
    left_join(new_data, by = ".row") %>%  # Merge with full dataset
    mutate(
      Call_type = .category,  # Preserve call type
      Call_prop = .epred  # Rename posterior prediction column
    ) %>%
    select(-.epred, -.category, -.row)  # Remove redundant columns
  
  # **Remove duplicated columns** (Keep only one version)
  PROP_pred_tidy <- PROP_pred_tidy %>%
    select(!ends_with(".y")) %>%
    rename_with(~ str_remove(., "\\.x$"))
  
  return(PROP_pred_tidy)
}

# ==============================================================================
#  FUNCTION TO FIND SEMI, PEAK DIG, AND FULL AGES FROM POSTERIOR SAMPLES
# ==============================================================================
find_transition_ages <- function(call_data, ID_DATA, num_samples = NUM_SAMPLES_PER_ID) {
  call_data <- call_data %>%
    mutate(Call_type = case_when(
      .category %in% c("SumDIG") ~ "DIG",
      .category %in% c("SumBEG") ~ "REP",
      .category %in% c("SumCC") ~ "CC",
      TRUE ~ as.character(.category)  # Keep other types unchanged
    )) %>%
    select(-.category)  # Remove redundant column after renaming
  
  # Group & reshape to wide format
  call_data <- call_data %>%
    group_by(ID, REC_AGE, Call_type, .draw) %>%
    summarise(Call_prop = mean(Call_prop), .groups = "drop") %>% 
    pivot_wider(names_from = Call_type, values_from = Call_prop) %>%
    na.omit()  # Remove missing values before further calculations

  
  # Normalize proportions per age so they sum to 1
  call_data <- call_data %>%
    mutate(TOTAL_PROP = DIG + REP + CC) %>%
    mutate(
      DIG = DIG / TOTAL_PROP,
      REP = REP / TOTAL_PROP,
      CC = CC / TOTAL_PROP
    ) %>%
    select(-TOTAL_PROP)
  
  # Ensure DIG > REP with a minimum threshold and enforce age constraints
  semi_age <- call_data %>%
    filter(DIG > REP + PROP_DIFF) %>%
    group_by(.draw) %>%
    summarise(SEMI_AGE = min(REC_AGE, na.rm = TRUE), .groups = "drop") %>%
    mutate(SEMI_AGE = ifelse(SEMI_AGE < 30, NA, SEMI_AGE))
  
  # Find the peak DIG age (max DIG proportion) after SEMI_AGE
  peak_dig_age <- call_data %>%
    left_join(semi_age, by = ".draw") %>%  # Join early to access SEMI_AGE
    filter(REC_AGE > SEMI_AGE) %>%  # Ensure PEAK_DIG_AGE is after SEMI_AGE
    group_by(.draw) %>%
    filter(DIG == max(DIG, na.rm = TRUE)) %>%  # Find max DIG proportion
    summarise(PEAK_DIG_AGE = first(REC_AGE), .groups = "drop")
  
  full_age <- call_data %>%
    left_join(peak_dig_age, by = ".draw") %>%  # Join early to access PEAK_DIG_AGE
    filter(REC_AGE > PEAK_DIG_AGE) %>%  # Ensure FULL_AGE is after PEAK_DIG_AGE
    filter(CC > DIG + PROP_DIFF) %>%  # DIG must be lower than CC
    group_by(.draw) %>%
    summarise(FULL_AGE = min(REC_AGE, na.rm = TRUE), .groups = "drop")
  
  age_data <- semi_age %>%
    left_join(peak_dig_age, by = ".draw") %>%  # Merge peak DIG age
    left_join(full_age, by = ".draw")          # Merge full age
  
  # **Ensure we have exactly 100 samples per ID**
  transition_ages <- age_data %>%
    drop_na()   # Remove NAs
    
  min_samples <- nrow(transition_ages)
    
  transition_ages <- transition_ages %>%
    slice_sample(n = min(min_samples, num_samples), replace = F)  # Adjust if fewer valid rows exist
  
  # Add metadata
  transition_ages <- transition_ages %>%
    mutate(
      ID = unique(ID_DATA$ID),
      TREATMENT = unique(ID_DATA$TREATMENT),
      SEX = unique(ID_DATA$SEX),
      MAT_STAT = unique(ID_DATA$MAT_STAT),
      WEIGHT = unique(ID_DATA$WEIGHT),
      COMP = unique(ID_DATA$COMP),
      GS = unique(ID_DATA$GS),
      RAIN = unique(ID_DATA$RAIN),
      MOTHER_ID = unique(ID_DATA$MOTHER_ID),
      LITTER_CODE = unique(ID_DATA$LITTER_CODE)
    )
  
  return(transition_ages)
}

# ==============================================================================
#  APPLY TO ALL INDIVIDUALS: Extract Posterior Call Proportions & Transitions
# ==============================================================================
# Create dataset with unique individuals but sample means for standardisation
FLUT_ID_data <- PROP_data %>%
  group_by(ID) %>%
  summarise(
    WEIGHT = mean(PROP_data$WEIGHT_DIFF_PER),
    COMP = mean(PROP_data$COMP_NORM),
    GS = mean(PROP_data$GROUPSIZE),
    RAIN = mean(PROP_data$MonthlyRainfall),
    WEIGHT_z = 0, COMP_z = 0, GS_z = 0, RAIN_z = 0,
    MAT_STAT = first(MaternalRank),
    MOTHER_ID = first(MOTHER_ID),
    LITTER_CODE = first(LITTER_CODE),
    ID = first(ID),
    TREATMENT = as.factor(first(TREATMENT)),
    SEX = first(SEX)
  )

# Initialize an empty tibble to store results
full_transition_data <- tibble()

# Loop through all subjects
for (i in 1:nrow(FLUT_ID_data)) {
  row <- FLUT_ID_data[i, ]
  posterior_call_props <- extract_posterior_call_props(row)
  transition_ages <- find_transition_ages(posterior_call_props, row, num_samples = NUM_SAMPLES_PER_ID)

  full_transition_data <- bind_rows(full_transition_data, transition_ages)

  # Print progress
  percent_done <- round((i / nrow(FLUT_ID_data)) * 100, 2)
  print(paste("Finished processing ID:", row$ID, "-", percent_done, "% complete"))
}


# Save posterior transition ages
saveRDS(full_transition_data, "../data/MILESTONES_POSTERIOR.rds")

### cleanup ###
rm(full_transition_data, B_prop, PROP_data)