#------------------------------setup--------------------------------------------

library(lme4)
library(dplyr)
library(ggplot2)

#-------------------------------read in dfs-------------------------------------

B5 <- read.csv("corr_B5.csv",
               colClasses = c(
                 "CLY25_B5"   = "factor",
                 "Rep"        = "factor",
                 "Genotype"   = "factor",
                 "accession"  = "factor",
                 "species"    = "factor",
                 "species_mex"= "factor",
                 "origin"     = "factor",
                 "notes1"     = "factor",
                 "notes2"     = "factor",
                 "GDDTS_corr" = "numeric",
                 "GDDTA_corr" = "numeric",
                 "PH_corr"    = "numeric",
                 "EH_corr"    = "numeric",
                 "EN_corr"    = "numeric",
                 "NBR_corr"   = "numeric",
                 "SPAD_corr"  = "numeric"
               )
)

#------------------------select traits------------------------------------------

# traits are the numeric columns
trait_cols_B5 <- names(B5)[sapply(B5, is.numeric)]

# remove checks
clean_B5 <- B5 %>%
  filter(!(Genotype %in% c("B73","Purple Check")))

#------------------------run for all traits (lme4)------------------------------

results_B5_lme4 <- data.frame(
  trait = character(),
  H2 = numeric(),
  n_obs = numeric(),
  n_genotypes = numeric(),
  stringsAsFactors = FALSE
)

for (trait in trait_cols_B5) {
  this_data <- clean_B5 %>%
    dplyr::select(Genotype, Rep, all_of(trait)) %>%
    filter(!is.na(.data[[trait]]))
  
  if (nrow(this_data) < 20) next
  
  model_formula <- as.formula(paste(trait, "~ (1|Genotype) + (1|Rep)"))
  tryCatch({
    model <- lmer(model_formula, data = this_data, REML = TRUE)
    vc <- as.data.frame(VarCorr(model))
    Vg <- vc$vcov[vc$grp == "Genotype"]
    Ve <- vc$vcov[vc$grp == "Residual"]
    H2 <- Vg / (Vg + Ve)
    
    results_B5_lme4 <- rbind(results_B5_lme4, data.frame(
      trait = trait,
      H2 = round(H2, 3),
      n_obs = nrow(this_data),
      n_genotypes = n_distinct(this_data$Genotype)
    ))
  }, error = function(e) {
    message(paste("lme4 model failed for trait:", trait))
  })
}

# look at it
results_B5_lme4

# plot it
ggplot(results_B5_lme4, aes(x = reorder(trait, -H2), y = H2)) +
  geom_col(fill = "steelblue") +
  labs(x = "Trait", y = "H² (2025)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))

results_B5_lme4$trait

#--------------------------- D4 --------------------------------------------------------------------------------------------

# read in dfs
D4 <- read.csv("corr_D4.csv",
               colClasses = c(
                 "CLY23_D4"   = "factor",
                 "Rep"        = "factor",
                 "Genotype"   = "factor",
                 "Species"    = "factor",
                 "origin"     = "factor",
                 "notes1"     = "factor",
                 "notes2"     = "factor",
                 "GDDTS_corr" = "numeric",
                 "GDDTA_corr" = "numeric",
                 "PH_corr"    = "numeric",
                 "EH_corr"    = "numeric",
                 "NBR_corr"   = "numeric",
                 "SPAD1_corr" = "numeric",
                 "SPAD2_corr" = "numeric"
               )
)
trait_cols_D4 <- names(D4)[sapply(D4, is.numeric)]

# remove checks
clean_D4 <- D4 %>%
  filter(!(Genotype %in% c("B73","Purple Check")))

#------------------------run for all traits (lme4)------------------------------

results_D4_lme4 <- data.frame(
  trait = character(),
  H2 = numeric(),
  n_obs = numeric(),
  n_genotypes = numeric(),
  stringsAsFactors = FALSE
)

for (trait in trait_cols_D4) {
  this_data <- clean_D4 %>%
    dplyr::select(Genotype, Rep, all_of(trait)) %>%
    filter(!is.na(.data[[trait]]))
  
  if (nrow(this_data) < 20) next
  
  model_formula <- as.formula(paste(trait, "~ (1|Genotype) + (1|Rep)"))
  tryCatch({
    model <- lmer(model_formula, data = this_data, REML = TRUE)
    vc <- as.data.frame(VarCorr(model))
    Vg <- vc$vcov[vc$grp == "Genotype"]
    Ve <- vc$vcov[vc$grp == "Residual"]
    H2 <- Vg / (Vg + Ve)
    
    results_D4_lme4 <- rbind(results_D4_lme4, data.frame(
      trait = trait,
      H2 = round(H2, 3),
      n_obs = nrow(this_data),
      n_genotypes = n_distinct(this_data$Genotype)
    ))
  }, error = function(e) {
    message(paste("lme4 model failed for trait:", trait))
  })
}

# look at it
results_D4_lme4

# plot it
ggplot(results_D4_lme4, aes(x = reorder(trait, -H2), y = H2)) +
  geom_col(fill = "forestgreen") +
  labs(x = "Trait", y = "H² (2023)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))

results_D4_lme4$trait
