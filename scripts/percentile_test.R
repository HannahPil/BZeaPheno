library(tidyverse)

#----------------------------- read + filter -----------------------------------

B5 <- read.csv("corr_B5.csv")

B5_f <- B5 %>%
  filter((is.na(notes1) | notes1 == "") & (is.na(notes2) | notes2 == "") &
           !is.na(species) & species != "Check") %>%
  filter(!is.na(SPAD_corr)) %>%
  mutate(
    species = factor(species, levels = c("B73", sort(setdiff(unique(species), "B73"))))
  )

#----------------------------- summaries: mean/sd + 95th -----------------------

group_means <- B5_f %>%
  group_by(species) %>%
  summarise(
    mean_spad = mean(SPAD_corr),
    sd_spad   = sd(SPAD_corr),
    q95_spad  = quantile(SPAD_corr, 0.95),
    .groups   = "drop"
  )

#----------------------------- sample counts -----------------------------------

sample_counts <- B5_f %>%
  count(species, name = "n_samples")

print(sample_counts)
cat("Total samples used:", nrow(B5_f), "\n")

#----------------------------- plot: jitter + mean±sd + mean point + 95th point -

p <- ggplot(B5_f, aes(x = species, y = SPAD_corr, fill = species)) +
  geom_jitter(
    shape = 21, size = 4, alpha = 0.8,
    position = position_jitter(width = 0.1)
  ) +
  # mean +/- sd (error bar)
  geom_errorbar(
    data = group_means,
    aes(x = species, ymin = mean_spad - sd_spad, ymax = mean_spad + sd_spad),
    width = 0, size = 1,
    inherit.aes = FALSE
  ) +
  # mean point (black)
  geom_point(
    data = group_means,
    aes(x = species, y = mean_spad),
    shape = 16, size = 3,
    inherit.aes = FALSE
  ) +
  # 95th percentile point (blue)
  geom_point(
    data = group_means,
    aes(x = species, y = q95_spad),
    shape = 17, size = 3,
    inherit.aes = FALSE,
    color = "blue"
  ) +
  labs(
    x = "Species",
    y = "SPAD (corrected)",
    title = "SPAD_corr by species (B5)",
    subtitle = "black point = mean; error bar = mean ± SD; blue triangle = 95th percentile"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

# print plot
p

# save plot
ggsave(
  filename = "output/B5_SPAD_corr_mean_sd_and_q95_by_species.png",
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)


#----------------species_mex-----------------------
library(tidyverse)
library(quantreg)

#----------------------------- read + filter -----------------------------------

B5 <- read.csv("corr_B5.csv")

B5_f <- B5 %>%
  filter((is.na(notes1) | notes1 == "") & (is.na(notes2) | notes2 == "") &
           !is.na(species_mex) & species_mex != "Check") %>%
  filter(!is.na(SPAD_corr)) %>%
  mutate(
    species_mex = factor(
      species_mex,
      levels = c("B73", sort(setdiff(unique(species_mex), "B73")))
    )
  )

#----------------------------- summaries for plotting --------------------------

group_means <- B5_f %>%
  group_by(species_mex) %>%
  summarise(
    mean_spad = mean(SPAD_corr),
    sd_spad   = sd(SPAD_corr),
    q95_spad  = quantile(SPAD_corr, 0.95),
    .groups   = "drop"
  )

# sample counts used
print(B5_f %>% count(species_mex, name = "n_samples"))
cat("Total samples used:", nrow(B5_f), "\n")

#----------------------------- tail test (quantile regression) -----------------
# compare upper tail (top 5%) vs B73

B5_f$species_mex <- relevel(B5_f$species_mex, ref = "B73")

rq_model <- rq(SPAD_corr ~ species_mex, data = B5_f, tau = 0.95)
rq_sum   <- summary(rq_model, se = "boot")

cat("\n--- Quantile regression (tau = 0.95; tail test vs B73) ---\n")
print(rq_sum)

# extract p-values for species vs B73 (robust to column naming differences)
coef_df <- as.data.frame(rq_sum$coefficients) %>%
  tibble::rownames_to_column("term")

# first numeric column = estimate; last column = p-value (works across quantreg outputs)
est_col <- names(coef_df)[2]
p_col   <- names(coef_df)[ncol(coef_df)]

tail_results <- coef_df %>%
  filter(term != "(Intercept)") %>%
  transmute(
    species_mex = gsub("^species_mex", "", term),
    estimate_q95_diff = .data[[est_col]],
    p_value = .data[[p_col]]
  )

print(tail_results)


# significance labels for annotation
tail_results <- tail_results %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  filter(sig != "")

# set annotation y just above the 95th percentile point for each species
y_offset <- 0.03 * (max(B5_f$SPAD_corr, na.rm = TRUE) - min(B5_f$SPAD_corr, na.rm = TRUE))

annot_df <- group_means %>%
  mutate(species_mex = as.character(species_mex)) %>%
  inner_join(tail_results, by = "species_mex") %>%
  mutate(
    species_mex = factor(species_mex, levels = levels(B5_f$species_mex)),
    y_star = q95_spad + y_offset
  )

#----------------------------- plot: jitter + mean±sd + mean point + 95th + stars

p <- ggplot(B5_f, aes(x = species_mex, y = SPAD_corr, fill = species_mex)) +
  geom_jitter(
    shape = 21, size = 4, alpha = 0.8,
    position = position_jitter(width = 0.1),
    color = "grey30"
  ) +
  geom_errorbar(
    data = group_means,
    aes(x = species_mex, ymin = mean_spad - sd_spad, ymax = mean_spad + sd_spad),
    width = 0, size = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_means,
    aes(x = species_mex, y = mean_spad, shape = "Mean"),
    size = 3,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_point(
    data = group_means,
    aes(x = species_mex, y = q95_spad, shape = "95th percentile"),
    size = 4,
    inherit.aes = FALSE,
    color = "blue"
  ) +
  geom_text(
    data = annot_df,
    aes(x = species_mex, y = y_star, label = sig),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_shape_manual(
    name = "Summary statistic",
    values = c("Mean" = 16, "95th percentile" = 17)
  ) +
  labs(
    x = "Species (mex)",
    y = "SPAD (corrected)",
    title = "B5 (2025): SPAD_corr by species_mex",
    subtitle = "Blue triangle = 95th percentile; stars = quantile regression tail test vs B73 (tau = 0.95)",
    shape = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  guides(fill = "none")

# print
p

# save
ggsave(
  filename = "output/B5_SPAD_corr_tailtest_by_species_mex.png",
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)

#-------------------predicted N------------------------------------
#----------------------------- read + filter -----------------------------------

B5 <- read.csv("corr_B5.csv")

B5_f <- B5 %>%
  filter((is.na(notes1) | notes1 == "") & (is.na(notes2) | notes2 == "") &
           !is.na(species_mex) & species_mex != "Check") %>%
  filter(!is.na(Predicted_N_corr)) %>%
  mutate(
    species_mex = factor(
      species_mex,
      levels = c("B73", sort(setdiff(unique(species_mex), "B73")))
    )
  )

#----------------------------- summaries for plotting --------------------------

group_means <- B5_f %>%
  group_by(species_mex) %>%
  summarise(
    mean_predN = mean(Predicted_N_corr),
    sd_predN   = sd(Predicted_N_corr),
    q95_predN  = quantile(Predicted_N_corr, 0.95),
    .groups    = "drop"
  )

# sample counts used
print(B5_f %>% count(species_mex, name = "n_samples"))
cat("Total samples used:", nrow(B5_f), "\n")

#----------------------------- tail test (quantile regression) -----------------
# compare upper tail (top 5%) vs B73

B5_f$species_mex <- relevel(B5_f$species_mex, ref = "B73")

rq_model <- rq(Predicted_N_corr ~ species_mex, data = B5_f, tau = 0.95)
rq_sum   <- summary(rq_model, se = "boot")

cat("\n--- Quantile regression (tau = 0.95; tail test vs B73) ---\n")
print(rq_sum)

# extract estimates + p-values (robust to column naming differences)
coef_df <- as.data.frame(rq_sum$coefficients) %>%
  tibble::rownames_to_column("term")

est_col <- names(coef_df)[2]
p_col   <- names(coef_df)[ncol(coef_df)]

tail_results <- coef_df %>%
  filter(term != "(Intercept)") %>%
  transmute(
    species_mex = gsub("^species_mex", "", term),
    estimate_q95_diff = .data[[est_col]],
    p_value = .data[[p_col]]
  )

print(tail_results)

# significance labels for annotation
tail_results <- tail_results %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  ) %>%
  filter(sig != "")

# set annotation y just above the 95th percentile point for each species
y_offset <- 0.03 * (max(B5_f$Predicted_N_corr, na.rm = TRUE) - min(B5_f$Predicted_N_corr, na.rm = TRUE))

annot_df <- group_means %>%
  mutate(species_mex = as.character(species_mex)) %>%
  inner_join(tail_results, by = "species_mex") %>%
  mutate(
    species_mex = factor(species_mex, levels = levels(B5_f$species_mex)),
    y_star = q95_predN + y_offset
  )

#----------------------------- plot: jitter + mean±sd + mean point + 95th + stars

p <- ggplot(B5_f, aes(x = species_mex, y = Predicted_N_corr, fill = species_mex)) +
  geom_jitter(
    shape = 21, size = 4, alpha = 0.8,
    position = position_jitter(width = 0.1),
    color = "grey30"
  ) +
  geom_errorbar(
    data = group_means,
    aes(x = species_mex, ymin = mean_predN - sd_predN, ymax = mean_predN + sd_predN),
    width = 0, size = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_means,
    aes(x = species_mex, y = mean_predN, shape = "Mean"),
    size = 3,
    inherit.aes = FALSE,
    color = "black"
  ) +
  geom_point(
    data = group_means,
    aes(x = species_mex, y = q95_predN, shape = "95th percentile"),
    size = 4,
    inherit.aes = FALSE,
    color = "blue"
  ) +
  geom_text(
    data = annot_df,
    aes(x = species_mex, y = y_star, label = sig),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_shape_manual(
    name = "Summary statistic",
    values = c("Mean" = 16, "95th percentile" = 17)
  ) +
  labs(
    x = "Species (mex)",
    y = "Predicted N (corrected)",
    title = "B5 (2025): Predicted_N_corr by species_mex",
    subtitle = "Blue triangle = 95th percentile; stars = quantile regression tail test vs B73 (tau = 0.95)",
    shape = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  guides(fill = "none")

# print
p

# save
ggsave(
  filename = "output/B5_Predicted_N_corr_tailtest_by_species_mex.png",
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)