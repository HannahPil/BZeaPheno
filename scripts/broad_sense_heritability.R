#------------------------------setup--------------------------------------------

library(lme4)
library(dplyr)
library(ggplot2)

#------------------------------- B5 --------------------------------------------

# read in df
B5 <- read.csv("corr_B5.csv")

# keep types sensible
B5 <- B5 %>%
  mutate(
    Plot        = factor(Plot),
    Rep         = factor(Rep),
    Genotype    = factor(Genotype),
    accession   = factor(accession),
    species     = factor(species),
    species_mex = factor(species_mex),
    origin      = factor(origin),
    notes1      = factor(notes1),
    notes2      = factor(notes2)
  )

# traits are ONLY columns ending in _corr
trait_cols_B5 <- names(B5)[grepl("_corr$", names(B5))]

# remove checks
clean_B5 <- B5 %>%
  filter(!(Genotype %in% c("B73", "Purple Check")))

# run for all traits (lme4)
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
    
    results_B5_lme4 <- rbind(
      results_B5_lme4,
      data.frame(
        trait = trait,
        H2 = round(H2, 3),
        n_obs = nrow(this_data),
        n_genotypes = n_distinct(this_data$Genotype),
        stringsAsFactors = FALSE
      )
    )
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

# list trait names
results_B5_lme4$trait

#------------------------------- D4 --------------------------------------------

# read in df
D4 <- read.csv("corr_D4.csv")

# keep types sensible
D4 <- D4 %>%
  mutate(
    Plot     = factor(Plot),
    Rep      = factor(Rep),
    Genotype = factor(Genotype),
    Species  = factor(Species),
    origin   = factor(origin),
    notes1   = factor(notes1),
    notes2   = factor(notes2)
  )

# traits are ONLY columns ending in _corr
trait_cols_D4 <- names(D4)[grepl("_corr$", names(D4))]

# remove checks
clean_D4 <- D4 %>%
  filter(!(Genotype %in% c("B73", "Purple Check")))

# run for all traits (lme4)
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
    
    results_D4_lme4 <- rbind(
      results_D4_lme4,
      data.frame(
        trait = trait,
        H2 = round(H2, 3),
        n_obs = nrow(this_data),
        n_genotypes = n_distinct(this_data$Genotype),
        stringsAsFactors = FALSE
      )
    )
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

# list trait names
results_D4_lme4$trait

#------------------------ combine plots ----------------------------------------

library(patchwork)

p_B5 <- ggplot(results_B5_lme4, aes(x = reorder(trait, -H2), y = H2)) +
  geom_col(fill = "steelblue") +
  labs(x = "Trait", y = "H² (2025)", title = "CLY25-B5") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))

p_D4 <- ggplot(results_D4_lme4, aes(x = reorder(trait, -H2), y = H2)) +
  geom_col(fill = "forestgreen") +
  labs(x = "Trait", y = "H² (2023)", title = "CLY23-D4") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))

# combine into two panels
p_combined <- p_B5 + p_D4 + plot_layout(ncol = 2)

# PRINT FIRST
p_combined

#save
ggsave(
  filename = "output/heritability_B5_D4_panels.png",
  plot = p_combined,
  width = 7,
  height = 5,
  dpi = 300
)
