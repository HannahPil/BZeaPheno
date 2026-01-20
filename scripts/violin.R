library(tidyverse)
library(ggpattern)

#------------------------------inputs-------------------------------------------

species_order <- c("Dura", "Nobo", "Mesa", "Chal", "B73", "Bals", "Zdip", "Hueh", "Zlux")

species_colors <- c(
  "Dura" = "#609bfe",
  "Nobo" = "#609bfe",
  "Mesa" = "#609bfe",
  "Chal" = "#609bfe",
  "B73"  = "#03bec4",
  "Bals" = "#f364e2",
  "Zdip" = "#f8756d",
  "Hueh" = "#b69d00",
  "Zlux" = "#00b837"
)

#------------------------------read + filter------------------------------------

B5 <- read.csv("corr_B5.csv")
D4 <- read.csv("corr_D4.csv")

# match your filtering logic: notes empty + species not check + trait not NA
B5f <- B5 %>%
  filter((is.na(notes1) | notes1 == "") & (is.na(notes2) | notes2 == "") &
           !is.na(species) & species != "Check") %>%
  filter(!is.na(GDDTS_corr) & !is.na(GDDTA_corr)) %>%
  mutate(
    year    = "2025",
    species = factor(species, levels = species_order)
  )

D4f <- D4 %>%
  filter((is.na(notes1) | notes1 == "") & (is.na(notes2) | notes2 == "") &
           !is.na(Species) & Species != "Check") %>%
  filter(!is.na(GDDTS_corr) & !is.na(GDDTA_corr)) %>%
  mutate(
    year    = "2023",
    species = factor(Species, levels = species_order)
  ) %>%
  select(-Species)

#------------------------------combine long-------------------------------------

combined_long <- bind_rows(
  B5f %>% select(species, year, GDDTS_corr, GDDTA_corr),
  D4f %>% select(species, year, GDDTS_corr, GDDTA_corr)
) %>%
  pivot_longer(
    cols = c(GDDTS_corr, GDDTA_corr),
    names_to = "trait",
    values_to = "value"
  ) %>%
  mutate(
    trait = recode(trait,
                   GDDTS_corr = "GDDTS (corrected)",
                   GDDTA_corr = "GDDTA (corrected)")
  )

#------------------------------plot 1: grid------------------------------------

p_grid <- ggplot(combined_long, aes(x = species, y = value, fill = species)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = mean_sdl, mult = 1,
               geom = "pointrange", color = "black") +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  facet_grid(trait ~ year, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Corrected GDDTS and GDDTA Distributions by Species (2023 vs 2025)",
    x = "Species",
    y = "GDD"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 20, hjust = 0.5),
    strip.text   = element_text(size = 16)
  ) +
  scale_y_continuous(breaks = seq(0, 1300, 100))

# print
p_grid

# save
ggsave(
  filename = "output/corr_GDDTS_GDDTA_violin_grid_2023_2025.png",
  plot = p_grid,
  width = 16,
  height = 9,
  dpi = 300
)

#------------------------------plot 2: dodged----------------------------------

pd <- position_dodge(width = 0.8)

p_pattern <- ggplot(
  combined_long,
  aes(
    x = species,
    y = value,
    fill = species,
    pattern = year
  )
) +
  geom_violin_pattern(
    trim = FALSE,
    position = pd,
    aes(group = interaction(species, year)),
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.03,
    pattern_alpha = 0.6
  ) +
  stat_summary(
    fun.data = mean_sdl,
    mult = 1,
    geom = "pointrange",
    color = "black",
    position = pd,
    aes(group = interaction(species, year))
  ) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  scale_pattern_manual(
    values = c(
      "2023" = "stripe",
      "2025" = "none"
    )
  ) +
  facet_grid(trait ~ ., scales = "free_y", switch = "y") +
  theme_minimal() +
  guides(
    fill = "none",
    pattern = guide_legend(override.aes = list(fill = "white", pattern = c("stripe", "none")))
  ) +
  labs(
    title = "Corrected GDDTS and GDDTA by Species (striped = 2023, solid = 2025)",
    x = "Species",
    y = "GDD",
    pattern = "Year"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.placement = "outside",
    strip.background = element_blank()
  )

# print
p_pattern

# save
ggsave(
  filename = "output/corr_GDDTS_GDDTA_violin_stacked_patterned_by_year.png",
  plot = p_pattern,
  width = 12,
  height = 9,
  dpi = 300
)

#----------------- PH + EH (corrected) -----------------------------------------
# base filtered dfs (no trait-specific NA filtering yet)
B5_base <- B5 %>%
  filter((is.na(notes1) | notes1 == "") & (is.na(notes2) | notes2 == "") &
           !is.na(species) & species != "Check") %>%
  mutate(
    year    = "2025",
    species = factor(species, levels = species_order)
  )

D4_base <- D4 %>%
  filter((is.na(notes1) | notes1 == "") & (is.na(notes2) | notes2 == "") &
           !is.na(Species) & Species != "Check") %>%
  mutate(
    year    = "2023",
    species = factor(Species, levels = species_order)
  ) %>%
  select(-Species)

# build long df for PH/EH (corrected)
combined_long_h <- bind_rows(
  B5_base %>% select(species, year, PH_corr, EH_corr),
  D4_base %>% select(species, year, PH_corr, EH_corr)
) %>%
  pivot_longer(
    cols = c(PH_corr, EH_corr),
    names_to = "trait",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    trait = recode(trait,
                   PH_corr = "Plant height (corrected)",
                   EH_corr = "Ear height (corrected)"),
    year = factor(year, levels = c("2023", "2025"))
  )

# plot (same styling as your GDD plot)
p_pattern_h <- ggplot(
  combined_long_h,
  aes(x = species, y = value, fill = species, pattern = year)
) +
  geom_violin_pattern(
    trim = FALSE,
    position = pd,
    aes(group = interaction(species, year)),
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.03,
    pattern_alpha = 0.6
  ) +
  stat_summary(
    fun.data = mean_sdl,
    mult = 1,
    geom = "pointrange",
    color = "black",
    position = pd,
    aes(group = interaction(species, year))
  ) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  scale_pattern_manual(values = c("2023" = "stripe", "2025" = "none")) +
  facet_grid(trait ~ ., scales = "free_y", switch = "y") +
  theme_minimal() +
  guides(
    fill = "none",
    pattern = guide_legend(
      override.aes = list(
        fill = "white",
        pattern = c("stripe", "none"),
        pattern_fill = "black"
      )
    )
  ) +
  labs(
    title = "Corrected plant height and ear height by species (striped = 2023, solid = 2025)",
    x = "Species",
    y = "Height (cm)",
    pattern = "Year"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.placement = "outside",
    strip.background = element_blank()
  )

# print
p_pattern_h

# save (to output/)
ggsave(
  filename = "output/corr_PH_EH_violin_stacked_patterned_by_year.png",
  plot = p_pattern_h,
  width = 12,
  height = 9,
  dpi = 300
)

