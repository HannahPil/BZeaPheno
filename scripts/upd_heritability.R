library(dplyr)
library(ggplot2)
library(patchwork)
library(tidytext)

#------------------------------ helper -----------------------------------------

calc_h2_anova <- function(df, trait_cols, genotype_col = "Genotype", rep_col = "Rep",
                          drop_genotypes = c("B73", "Purple Check"),
                          min_obs = 20) {
  
  df <- df %>%
    filter(!(.data[[genotype_col]] %in% drop_genotypes))
  
  results <- data.frame(
    trait = character(),
    H2 = numeric(),
    n_obs = numeric(),
    n_genotypes = numeric(),
    n_reps = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (trait in trait_cols) {
    
    dat <- df %>%
      select(all_of(c(genotype_col, rep_col, trait))) %>%
      filter(!is.na(.data[[trait]]))
    
    if (nrow(dat) < min_obs) next
    
    n_reps <- length(unique(dat[[rep_col]]))
    if (n_reps < 2) next
    
    complete_genos <- dat %>%
      count(.data[[genotype_col]]) %>%
      filter(n == n_reps) %>%
      pull(.data[[genotype_col]])
    
    dat <- dat %>%
      filter(.data[[genotype_col]] %in% complete_genos)
    
    if (nrow(dat) < min_obs) next
    
    mod <- lm(as.formula(paste(trait, "~", genotype_col)), data = dat)
    aov_tab <- anova(mod)
    
    ss_g <- aov_tab$`Sum Sq`[1]
    ss_e <- aov_tab$`Sum Sq`[2]
    if (ss_e == 0) next
    
    results <- rbind(
      results,
      data.frame(
        trait = trait,
        H2 = round(ss_g / (ss_g + ss_e), 3),
        n_obs = nrow(dat),
        n_genotypes = length(unique(dat[[genotype_col]])),
        n_reps = n_reps,
        stringsAsFactors = FALSE
      )
    )
  }
  
  results
}

#------------------------------ trait map --------------------------------------

trait_map <- read.csv("data/trait_rename.csv", stringsAsFactors = FALSE)
names(trait_map) <- c("original", "upd", "category")

category_order <- c("FT", "Morphology", "Photosynthesis", "Stem")
category_colors <- c(
  "FT" = "#6a5acd",
  "Morphology" = "#2e8b57",
  "Photosynthesis" = "#ff8c00",
  "Stem" = "#4682b4"
)

#------------------------------- B5 (2025) -------------------------------------

B5 <- read.csv("data/corr_B5.csv")

B5 <- B5 %>%
  mutate(
    Plot = factor(Plot),
    Rep = factor(Rep),
    Genotype = factor(Genotype)
  )

results_B5_aov <- calc_h2_anova(
  B5,
  names(B5)[grepl("_corr$", names(B5))]
) %>%
  left_join(trait_map, by = c("trait" = "original")) %>%
  mutate(
    trait_label = upd,
    trait_category = factor(category, levels = category_order),
    year = "2025"
  )

#------------------------------- D4 (2023) -------------------------------------

D4 <- read.csv("data/corr_D4.csv")

D4 <- D4 %>%
  mutate(
    Plot = factor(Plot),
    Rep = factor(Rep),
    Genotype = factor(Genotype)
  )

results_D4_aov <- calc_h2_anova(
  D4,
  names(D4)[grepl("_corr$", names(D4))]
) %>%
  left_join(trait_map, by = c("trait" = "original")) %>%
  mutate(
    trait_label = upd,
    trait_category = factor(category, levels = category_order),
    year = "2023"
  )

#------------------------------- combined plot ---------------------------------

results_all <- bind_rows(results_D4_aov, results_B5_aov)
results_all$year <- factor(results_all$year, levels = c("2023", "2025"))


results_all <- results_all %>%
  mutate(trait_label_ord = reorder_within(trait_label, -H2, year))

p_combined <- ggplot(
  results_all,
  aes(x = trait_label_ord, y = H2, fill = trait_category)
) +
  geom_col() +
  facet_wrap(~ year, ncol = 2, scales = "free_x") +
  labs(
    x = "Trait",
    y = "H² (SS_G / SS_Total)",
    fill = "Category"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_fill_manual(values = category_colors, drop = FALSE) +
  scale_x_reordered() +
  theme_bw(base_size = 18) +
  theme(
    strip.background = element_rect(fill = "grey85", color = "grey40"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_combined

ggsave(
  filename = "output/heritability_2023_2025_panels_anovaSS_grouped.png",
  plot = p_combined,
  width = 9,
  height = 5,
  dpi = 300
)
