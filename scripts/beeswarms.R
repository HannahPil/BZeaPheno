library(tidyverse)
library(ggbeeswarm)

dir.create("output", showWarnings = FALSE)

#----------------------------- read in corrected dfs ---------------------------

B5 <- read.csv("data/corr_B5.csv")
D4 <- read.csv("data/corr_D4.csv")

#----------------------------- filter ------------------------------------------

B5_f <- B5 %>%
  filter(species != "Check" & notes1 == "" & notes2 == "")

D4_f <- D4 %>%
  filter(Species != "Check" & notes1 == "" & notes2 == "")

#----------------------------- plot data: avg non-B73, keep B73 reps -----------

B5_nonB73 <- B5_f %>%
  filter(Genotype != "B73") %>%
  mutate(species_plot = species_mex) %>%
  group_by(Genotype, species_plot) %>%
  summarise(
    SPAD_value = mean(SPAD_corr, na.rm = TRUE),
    .groups = "drop"
  )

B5_B73 <- B5_f %>%
  filter(Genotype == "B73") %>%
  transmute(
    Genotype,
    species_plot = species_mex,
    SPAD_value = SPAD_corr
  )

B5_plotdat <- bind_rows(B5_nonB73, B5_B73)

D4_nonB73 <- D4_f %>%
  filter(Genotype != "B73") %>%
  mutate(species_plot = species_mex) %>%
  group_by(Genotype, species_plot) %>%
  summarise(
    SPAD_value = mean(SPAD2_corr, na.rm = TRUE),
    .groups = "drop"
  )

D4_B73 <- D4_f %>%
  filter(Genotype == "B73") %>%
  transmute(
    Genotype,
    species_plot = species_mex,
    SPAD_value = SPAD2_corr
  )

D4_plotdat <- bind_rows(D4_nonB73, D4_B73)

#----------------------------- summarize mean ± sd -----------------------------

group_means_B5 <- B5_plotdat %>%
  group_by(species_plot) %>%
  summarise(
    mean_spad = mean(SPAD_value, na.rm = TRUE),
    sd_spad   = sd(SPAD_value, na.rm = TRUE),
    .groups = "drop"
  )

group_means_D4 <- D4_plotdat %>%
  group_by(species_plot) %>%
  summarise(
    mean_spad = mean(SPAD_value, na.rm = TRUE),
    sd_spad   = sd(SPAD_value, na.rm = TRUE),
    .groups = "drop"
  )

#----------------------------- combine + relabel --------------------------------

B5_plotdat$Experiment <- "2025"
D4_plotdat$Experiment <- "2023"

combined_plotdat <- bind_rows(
  B5_plotdat %>% rename(value = SPAD_value),
  D4_plotdat %>% rename(value = SPAD_value)
)

group_means_B5$Experiment <- "2025"
group_means_D4$Experiment <- "2023"

combined_means <- bind_rows(
  group_means_B5,
  group_means_D4
)

# enforce panel order (2023 left, 2025 right)
combined_plotdat$Experiment <- factor(combined_plotdat$Experiment,
                                      levels = c("2023", "2025"))

combined_means$Experiment <- factor(combined_means$Experiment,
                                    levels = c("2023", "2025"))

# enforce species order
species_order <- c("B73", "Bals", "Zdip", "Hueh", "Zlux", "Mex")

combined_plotdat$species_plot <- factor(combined_plotdat$species_plot,
                                        levels = species_order)

combined_means$species_plot <- factor(combined_means$species_plot,
                                      levels = species_order)

# custom species colors
species_colors <- c(
  "B73"  = "#00bfc4",
  "Bals" = "#f564e3",
  "Zdip" = "#f8766d",
  "Hueh" = "#b79f00",
  "Zlux" = "#0ebe43",
  "Mex"  = "#619cff"
)

#----------------------------- faceted plot ------------------------------------

p_combined <- ggplot(combined_plotdat,
                     aes(x = species_plot, y = value)) +
  geom_quasirandom(
    aes(fill = species_plot, color = species_plot),
    shape = 21,
    size = 3,
    stroke = 1,
    alpha = 0.5,
    width = 0.3
  ) +
  geom_errorbar(
    data = combined_means,
    aes(
      x = species_plot,
      ymin = mean_spad - sd_spad,
      ymax = mean_spad + sd_spad
    ),
    width = 0,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = combined_means,
    aes(x = species_plot, y = mean_spad),
    shape = 16,
    size = 4,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Experiment, ncol = 2) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  scale_color_manual(values = species_colors, drop = FALSE) +
  labs(
    x = "",
    y = "SPAD leaf greenness index",
    title = ""
  ) +
  theme_minimal(base_size = 17) +
  theme(
    axis.text.x = element_text(),
    
    # facet strip styling
    strip.background = element_rect(fill = "grey85", color = "grey40"),
    strip.text = element_text(size = 16, face = "bold"),
    
    # add panel border
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank()
  ) +
  guides(fill = "none", color = "none")

p_combined

#----------------------------- save --------------------------------------------

ggsave(
  filename = "output/SPAD_species_mex_2023_2025.png",
  plot = p_combined,
  width = 14,
  height = 6,
  dpi = 300
)


#-----------------------now for gddta
#----------------------------- plot data: avg non-B73, keep B73 reps (GDDTA) ----

# B5 (GDDTA_corr, genotype column is "Genotype")
B5_nonB73_gddta <- B5_f %>%
  filter(Genotype != "B73") %>%
  mutate(species_plot = species_mex) %>%
  group_by(Genotype, species_plot) %>%
  summarise(
    GDDTA_value = mean(GDDTA_corr, na.rm = TRUE),
    .groups = "drop"
  )

B5_B73_gddta <- B5_f %>%
  filter(Genotype == "B73") %>%
  transmute(
    Genotype,
    species_plot = species_mex,
    GDDTA_value = GDDTA_corr
  )

B5_plotdat_gddta <- bind_rows(B5_nonB73_gddta, B5_B73_gddta)

# D4 (GDDTA_corr, genotype column is "Genotype")
D4_nonB73_gddta <- D4_f %>%
  filter(Genotype != "B73") %>%
  mutate(species_plot = species_mex) %>%
  group_by(Genotype, species_plot) %>%
  summarise(
    GDDTA_value = mean(GDDTA_corr, na.rm = TRUE),
    .groups = "drop"
  )

D4_B73_gddta <- D4_f %>%
  filter(Genotype == "B73") %>%
  transmute(
    Genotype,
    species_plot = species_mex,
    GDDTA_value = GDDTA_corr
  )

D4_plotdat_gddta <- bind_rows(D4_nonB73_gddta, D4_B73_gddta)

#----------------------------- combine + relabel --------------------------------

B5_plotdat_gddta$Experiment <- "2025"
D4_plotdat_gddta$Experiment <- "2023"

combined_plotdat_gddta <- bind_rows(B5_plotdat_gddta, D4_plotdat_gddta)

# summarize mean ± sd
group_means_B5_gddta <- B5_plotdat_gddta %>%
  group_by(species_plot) %>%
  summarise(
    mean_val = mean(GDDTA_value, na.rm = TRUE),
    sd_val   = sd(GDDTA_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Experiment = "2025")

group_means_D4_gddta <- D4_plotdat_gddta %>%
  group_by(species_plot) %>%
  summarise(
    mean_val = mean(GDDTA_value, na.rm = TRUE),
    sd_val   = sd(GDDTA_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Experiment = "2023")

combined_means_gddta <- bind_rows(group_means_B5_gddta, group_means_D4_gddta)

# enforce panel order (2023 left, 2025 right)
combined_plotdat_gddta$Experiment <- factor(combined_plotdat_gddta$Experiment,
                                            levels = c("2023", "2025"))
combined_means_gddta$Experiment <- factor(combined_means_gddta$Experiment,
                                          levels = c("2023", "2025"))

# enforce species order + colors (assumes you already defined these)
combined_plotdat_gddta$species_plot <- factor(combined_plotdat_gddta$species_plot,
                                              levels = species_order)
combined_means_gddta$species_plot <- factor(combined_means_gddta$species_plot,
                                            levels = species_order)

#----------------------------- faceted plot ------------------------------------

p_gddta <- ggplot(combined_plotdat_gddta,
                  aes(x = species_plot, y = GDDTA_value)) +
  geom_quasirandom(
    aes(fill = species_plot, color = species_plot),
    shape = 21,
    size = 3,
    stroke = 1,
    alpha = 0.5,
    width = 0.3
  ) +
  geom_errorbar(
    data = combined_means_gddta,
    aes(
      x = species_plot,
      ymin = mean_val - sd_val,
      ymax = mean_val + sd_val
    ),
    width = 0,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = combined_means_gddta,
    aes(x = species_plot, y = mean_val),
    shape = 16,
    size = 4,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Experiment, ncol = 2) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  scale_color_manual(values = species_colors, drop = FALSE) +
  labs(
    x = "",
    y = "GDDTA",
    title = "GDDTA by species_mex"
  ) +
  theme_minimal(base_size = 17) +
  theme(
    axis.text.x = element_text(),
    strip.background = element_rect(fill = "grey85", color = "grey40"),
    strip.text = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.background = element_blank()
  ) +
  guides(fill = "none", color = "none")

p_gddta

#----------------------------- save --------------------------------------------

ggsave(
  filename = "output/GDDTA_species_mex_2023_2025.png",
  plot = p_gddta,
  width = 14,
  height = 6,
  dpi = 300
)
