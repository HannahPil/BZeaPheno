#--------------------RIDGE
#--------------------RIDGE
#--------------------RIDGE
#--------------------RIDGE

library(tidyverse)
library(ggridges)

dir.create("output", showWarnings = FALSE)

#----------------------------- read in corrected dfs ---------------------------

B5 <- read.csv("data/corr_B5.csv")
D4 <- read.csv("data/corr_D4.csv")

#----------------------------- filter ------------------------------------------

B5_f <- B5 %>%
  filter(species != "Check" & notes1 == "" & notes2 == "")

D4_f <- D4 %>%
  filter(Species != "Check" & notes1 == "" & notes2 == "")

#----------------------------- species order + colors --------------------------

species_order <- c("B73", "Bals", "Zdip", "Hueh", "Zlux", "Mex")

species_colors <- c(
  "B73"  = "#00bfc4",
  "Bals" = "#f564e3",
  "Zdip" = "#f8766d",
  "Hueh" = "#b79f00",
  "Zlux" = "#0ebe43",
  "Mex"  = "#619cff"
)

shape_vals <- setNames(rep(21, length(species_order)), species_order)

#----------------------------- helper: avg non-B73, keep B73 reps --------------

make_plotdat <- function(df, value_col, experiment_label) {
  nonB73 <- df %>%
    filter(Genotype != "B73") %>%
    transmute(
      Genotype,
      species_plot = species_mex,
      value = .data[[value_col]]
    ) %>%
    group_by(Genotype, species_plot) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  
  B73 <- df %>%
    filter(Genotype == "B73") %>%
    transmute(
      Genotype,
      species_plot = species_mex,
      value = .data[[value_col]]
    )
  
  bind_rows(nonB73, B73) %>%
    mutate(
      Experiment = experiment_label,
      Experiment = factor(Experiment, levels = c("2023", "2025")),
      species_plot = factor(species_plot, levels = species_order)
    )
}

#==============================================================================
# SPAD RIDGES
#   2025 (B5): SPAD_corr
#   2023 (D4): SPAD2_corr
#==============================================================================

spad_2025 <- make_plotdat(B5_f, "SPAD_corr",  "2025")
spad_2023 <- make_plotdat(D4_f, "SPAD2_corr", "2023")
spad_all  <- bind_rows(spad_2023, spad_2025)

# reverse species order for ridges
spad_all$species_plot <- factor(
  spad_all$species_plot,
  levels = rev(species_order)
)

p_spad_ridges <- ggplot(
  spad_all,
  aes(x = value, y = species_plot, fill = species_plot)
) +
  # layer 1: fill + jittered points, but no outline (so points aren't "on top of the line")
  geom_density_ridges(
    aes(
      point_color = species_plot,
      point_fill  = species_plot,
      point_shape = species_plot
    ),
    alpha = 0.2,
    color = NA,
    point_alpha = 0.5,
    jittered_points = TRUE
  ) +
  # layer 2: outline only, drawn on top
  geom_density_ridges(
    alpha = 0,
    fill = NA,
    color = "black",
    linewidth = 0.5,
    jittered_points = FALSE
  ) +
  facet_wrap(~ Experiment, ncol = 2) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_color", values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_fill",  values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_shape", values = shape_vals) +
  labs(
    x = "SPAD leaf greenness index",
    y = "Species",
    title = "SPAD by species_mex"
  ) +
  scale_x_continuous(limits = c(-10, 80)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none",
    text = element_text(size = 17)
  )

p_spad_ridges


ggsave(
  filename = "output/SPAD_species_mex_2023_2025_ridges.png",
  plot = p_spad_ridges,
  width = 14,
  height = 6,
  dpi = 300
)

#==============================================================================
# GDDTA RIDGES
#   assumes GDDTA_corr exists in BOTH B5 and D4
#==============================================================================

gddta_2025 <- make_plotdat(B5_f, "GDDTA_corr", "2025")
gddta_2023 <- make_plotdat(D4_f, "GDDTA_corr", "2023")
gddta_all  <- bind_rows(gddta_2023, gddta_2025)

gddta_all$species_plot <- factor(
  gddta_all$species_plot,
  levels = rev(species_order)
)

p_gddta_ridges <- ggplot(
  gddta_all,
  aes(x = value, y = species_plot, fill = species_plot)
) +
  geom_density_ridges(
    aes(
      point_color = species_plot,
      point_fill  = species_plot,
      point_shape = species_plot
    ),
    alpha = 0.2,
    color = NA,
    point_alpha = 0.5,
    jittered_points = TRUE
  ) +
  geom_density_ridges(
    alpha = 0,
    fill = NA,
    color = "black",
    linewidth = 0.5,
    jittered_points = FALSE
  ) +
  facet_wrap(~ Experiment, ncol = 2) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_color", values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_fill",  values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_shape", values = shape_vals) +
  labs(
    x = "GDDTA",
    y = "Species",
    title = "GDDTA by species_mex"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none",
    text = element_text(size = 17)
  )

p_gddta_ridges


ggsave(
  filename = "output/GDDTA_species_mex_2023_2025_ridges.png",
  plot = p_gddta_ridges,
  width = 14,
  height = 6,
  dpi = 300
)

#--------------------------------------------------------------
library(tidyverse)
library(ggridges)

dir.create("output", showWarnings = FALSE)

#----------------------------- read B5 corrected df ----------------------------

B5 <- read.csv("data/corr_B5.csv")

B5_f <- B5 %>%
  filter(species != "Check" & notes1 == "" & notes2 == "")

#----------------------------- species order + colors --------------------------

species_order <- c("B73", "Bals", "Zdip", "Hueh", "Zlux", "Mex")

species_colors <- c(
  "B73"  = "#00bfc4",
  "Bals" = "#f564e3",
  "Zdip" = "#f8766d",
  "Hueh" = "#b79f00",
  "Zlux" = "#0ebe43",
  "Mex"  = "#619cff"
)

shape_vals <- setNames(rep(21, length(species_order)), species_order)

#----------------------------- avg non-B73, keep B73 reps ----------------------

pred_nonB73 <- B5_f %>%
  filter(Genotype != "B73") %>%
  group_by(Genotype, species_mex) %>%
  summarise(
    value = mean(Predicted_N_corr, na.rm = TRUE),
    .groups = "drop"
  )

pred_B73 <- B5_f %>%
  filter(Genotype == "B73") %>%
  transmute(
    Genotype,
    species_mex,
    value = Predicted_N_corr
  )

pred_all <- bind_rows(pred_nonB73, pred_B73) %>%
  mutate(
    species_plot = factor(species_mex, levels = rev(species_order))
  )

#----------------------------- ridge plot --------------------------------------

p_pred_ridges <- ggplot(
  pred_all,
  aes(x = value, y = species_plot, fill = species_plot)
) +
  # layer 1: fill + jittered points (no outline)
  geom_density_ridges(
    aes(
      point_color = species_plot,
      point_fill  = species_plot,
      point_shape = species_plot
    ),
    alpha = 0.2,
    color = NA,
    point_alpha = 0.6,
    jittered_points = TRUE
  ) +
  # layer 2: outline only (drawn on top)
  geom_density_ridges(
    alpha = 0,
    fill = NA,
    color = "black",
    linewidth = 0.5,
    jittered_points = FALSE
  ) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_color", values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_fill",  values = species_colors, drop = FALSE) +
  scale_discrete_manual(aesthetics = "point_shape", values = shape_vals) +
  labs(
    x = "Predicted N",
    y = "",
    title = ""
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none",
    text = element_text(size = 17)
  )

p_pred_ridges

ggsave(
  filename = "output/Predicted_N_corr_species_mex_B5_ridges.png",
  plot = p_pred_ridges,
  width = 8,
  height = 6,
  dpi = 300
)
