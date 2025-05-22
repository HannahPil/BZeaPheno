#------------------------------setup--------------------------------------------

install.packages("lme4")
install.packages("sommer")
library(lme4)
library(sommer)
library(dplyr)

bzeadf <- read.csv("./UPDATED_CLY23_D4_Fieldbook.csv")
head(bzeadf)

#------------------------shared cleaning and setup------------------------------

# Exclude non-trait columns
drop_cols <- c("CLY23_D4", "Rep", "old_genotype", "genotype_id",
               "Species", "Notes", "SPAD1_Date", "SPAD2_Date")
trait_cols <- setdiff(colnames(bzeadf), drop_cols)

# Clean and prepare data
clean_data <- bzeadf %>%
  filter(!(genotype_id %in% c("B73", "Purple Check"))) %>%
  mutate(
    genotype_id = factor(genotype_id),
    Rep = factor(Rep),
    across(c(ST, StPu, StPi), ~ ifelse(. == "Y", 1, ifelse(. == "N", 0, NA)), .names = "converted_{.col}"),
    Kinki = as.numeric(Kinki)
  )

# Replace original binary columns with converted versions
clean_data$ST   <- clean_data$converted_ST
clean_data$StPu <- clean_data$converted_StPu
clean_data$StPi <- clean_data$converted_StPi

# Drop converted_ columns
clean_data <- clean_data %>% select(-starts_with("converted_"))

#------------------------test lme4 model with DTS-------------------------------

# Remove checks and missing values
bzea_clean <- bzeadf %>%
  filter(!(genotype_id %in% c("B73", "Purple Check"))) %>%
  mutate(
    genotype_id = factor(genotype_id),
    Rep = factor(Rep),
    DTS = as.numeric(DTS)
  ) %>%
  filter(!is.na(DTS))  # remove missing DTS

# Fit the linear mixed model
dts_model <- lmer(DTS ~ (1|genotype_id) + (1|Rep), data = bzea_clean)

# Extract variance components
varcomp <- as.data.frame(VarCorr(dts_model))
vg <- varcomp$vcov[varcomp$grp == "genotype_id"]  # genetic variance
ve <- varcomp$vcov[varcomp$grp == "Residual"]     # residual variance

# Estimate H² with replication accounted for automatically by lmer
H2 <- vg / (vg + ve)
cat("Broad-sense heritability (H²) for DTS:", round(H2, 3), "\n")

#------------------------run for all traits (lme4)------------------------------

results_lme4 <- data.frame(
  trait = character(),
  H2 = numeric(),
  n_obs = numeric(),
  n_genotypes = numeric(),
  stringsAsFactors = FALSE
)

for (trait in trait_cols) {
  this_data <- clean_data %>%
    dplyr::select(genotype_id, Rep, all_of(trait)) %>%
    filter(!is.na(.data[[trait]]))
  
  if (nrow(this_data) < 20) next
  
  this_data[[trait]] <- as.numeric(this_data[[trait]])
  
  model_formula <- as.formula(paste(trait, "~ (1|genotype_id) + (1|Rep)"))
  tryCatch({
    model <- lmer(model_formula, data = this_data, REML = TRUE)
    vc <- as.data.frame(VarCorr(model))
    Vg <- vc$vcov[vc$grp == "genotype_id"]
    Ve <- vc$vcov[vc$grp == "Residual"]
    H2 <- Vg / (Vg + Ve)
    
    results_lme4 <- rbind(results_lme4, data.frame(
      trait = trait,
      H2 = round(H2, 3),
      n_obs = nrow(this_data),
      n_genotypes = n_distinct(this_data$genotype_id)
    ))
  }, error = function(e) {
    message(paste("lme4 model failed for trait:", trait))
  })
}

#look at it
results_lme4

#plot it
library(ggplot2)

ggplot(results_lme4, aes(x = reorder(trait, -H2), y = H2)) +
  geom_col(fill = "steelblue") +
  labs(x = "Trait", y = "H²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#---------------------------also trying sommer----------------------------------

results_sommer <- data.frame(
  trait = character(),
  H2_sommer = numeric(),
  n_obs = numeric(),
  n_genotypes = numeric(),
  stringsAsFactors = FALSE
)

for (trait in trait_cols) {
  this_data <- clean_data %>%
    dplyr::select(genotype_id, Rep, all_of(trait)) %>%
    filter(!is.na(.data[[trait]]))
  
  if (nrow(this_data) < 20) next
  
  this_data[[trait]] <- as.numeric(this_data[[trait]])
  this_data$genotype_id <- as.factor(this_data$genotype_id)
  this_data$Rep <- as.factor(this_data$Rep)
  
  tryCatch({
    model_sommer <- mmes(
      fixed = as.formula(paste(trait, "~ 1")),
      random = list(
        genotype_id = usm("genotype_id"),
        Rep = usm("Rep")
      ),
      data = this_data
    )
    
    Vg <- model_sommer$sigma["genotype_id", "Variance"]
    Ve <- model_sommer$sigma["units", "Variance"]
    H2 <- Vg / (Vg + Ve)
    
    results_sommer <- rbind(results_sommer, data.frame(
      trait = trait,
      H2_sommer = round(H2, 3),
      n_obs = nrow(this_data),
      n_genotypes = n_distinct(this_data$genotype_id)
    ))
  }, error = function(e) {
    message(paste("mmes model failed for trait:", trait))
  })
}

str(this_data)
summary(this_data)
#look at it
results_sommer

#ok that doesn't work because it can't handle 0 values in entire reps
#trying just with DTS/DTA:

# Subset and clean for DTS and DTA
traits_to_try <- c("DTS", "DTA")

# Create an empty results table
results_sommer_subset <- data.frame(
  trait = character(),
  H2_sommer = numeric(),
  n_genotypes = numeric(),
  stringsAsFactors = FALSE
)

for (trait in traits_to_try) {
  
  # Filter only relevant columns and rows with non-missing values
  this_data <- clean_data %>%
    dplyr::select(genotype_id, Rep, all_of(trait)) %>%
    filter(!is.na(.data[[trait]])) %>%
    mutate(
      genotype_id = factor(genotype_id),
      Rep = factor(Rep)
    )
  
  # Keep only genotypes with data in all 3 reps
  reps_per_geno <- this_data %>%
    group_by(genotype_id) %>%
    summarize(n_reps = n_distinct(Rep), .groups = "drop") %>%
    filter(n_reps == 3)
  
  this_data <- this_data %>%
    filter(genotype_id %in% reps_per_geno$genotype_id) %>%
    droplevels()
  
  # Skip if not enough data
  if (nrow(this_data) < 30) next
  
  tryCatch({
    model <- mmes(
      fixed = as.formula(paste(trait, "~ 1")),
      random = list(
        genotype_id = usm(~ genotype_id),
        Rep = usm(~ Rep)
      ),
      data = this_data
    )
    
    Vg <- model$sigma["genotype_id", "Variance"]
    Ve <- model$sigma["units", "Variance"]
    H2 <- Vg / (Vg + Ve)
    
    results_sommer_subset <- rbind(results_sommer_subset, data.frame(
      trait = trait,
      H2_sommer = round(H2, 3),
      n_genotypes = n_distinct(this_data$genotype_id)
    ))
  }, error = function(e) {
    message(paste("sommer failed for", trait, ":", e$message))
  })
}

# View the result
results_sommer_subset
head(clean_data)
head(this_data)
#-----------------------------compare results-----------------------------------

herit_both <- left_join(results_lme4, results_sommer, by = c("trait", "n_obs", "n_genotypes"))
herit_both
