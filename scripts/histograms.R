library(tidyverse)

#----------------------------- read in corrected dfs ---------------------------

B5 <- read.csv("corr_B5.csv")
D4 <- read.csv("corr_D4.csv")

#----------------------------- filter (match your example) ---------------------
# example filter: Species != "Check" & Notes == ""
# in these files, we’ll use species/Species and notes1/notes2

B5_f <- B5 %>%
  filter(species != "Check" & notes1 == "" & notes2 == "")

D4_f <- D4 %>%
  filter(Species != "Check" & notes1 == "" & notes2 == "")

#----------------------------- traits: ONLY _corr ------------------------------

traits_B5 <- names(B5_f)[grepl("_corr$", names(B5_f))]
traits_D4 <- names(D4_f)[grepl("_corr$", names(D4_f))]

#----------------------------- pivot to long for faceting ----------------------

B5_long <- B5_f %>%
  pivot_longer(cols = all_of(traits_B5), names_to = "trait", values_to = "value") %>%
  filter(!is.na(value))

D4_long <- D4_f %>%
  pivot_longer(cols = all_of(traits_D4), names_to = "trait", values_to = "value") %>%
  filter(!is.na(value))

#----------------------------- plots (one per field) ---------------------------

p_B5 <- ggplot(B5_long, aes(x = value)) +
  geom_histogram(color = "black", alpha = 0.7, bins = 30) +
  facet_wrap(~ trait, scales = "free") +
  labs(
    title = "B5: Distribution of corrected traits",
    x = "Trait value",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_D4 <- ggplot(D4_long, aes(x = value)) +
  geom_histogram(color = "black", alpha = 0.7, bins = 30) +
  facet_wrap(~ trait, scales = "free") +
  labs(
    title = "D4: Distribution of corrected traits",
    x = "Trait value",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# print plots
p_B5
p_D4

# save plots (to output/)
ggsave(
  filename = "output/B5_corr_trait_hist_panels.png",
  plot = p_B5,
  width = 14,
  height = 8,
  dpi = 300
)

ggsave(
  filename = "output/D4_corr_trait_hist_panels.png",
  plot = p_D4,
  width = 14,
  height = 8,
  dpi = 300
)
