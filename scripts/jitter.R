library(tidyverse)
library(dplyr)
library(ggplot2)

bzeadf <- read.csv("./UPDATED_CLY23_D4_FieldBook.csv")

bzeadf_f <- bzeadf %>% filter(Species != "Check" & Notes == "")
bzeadf_ff <- bzeadf_f %>%
  mutate(LA = BW * BL * 0.75)

group_means <- bzeadf_ff %>%
  group_by(Species) %>%
  summarise(
    mean_PH = mean(PH, na.rm = TRUE),
    sd_PH = sd(PH, na.rm = TRUE),
    .groups = "drop"
  )

# Jitter plot: PH by Species
ggplot(bzeadf_ff, aes(x = Species, y = PH, fill = Species)) +
  geom_jitter(
    shape = 21, size = 1, alpha = 0.8,
    position = position_jitter(width = 0.1)
  ) +
  geom_errorbar(
    data = group_means,
    aes(x = Species, ymin = mean_PH - sd_PH, ymax = mean_PH + sd_PH),
    width = 0.1, size = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_means,
    aes(x = Species, y = mean_PH),
    shape = 16, size = 1,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Species",
    y = "Plant Height (PH)",
    title = "Plant Height by Species"
  ) +
  theme_minimal()
