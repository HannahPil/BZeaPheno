library(ggplot2)
library(ggbeeswarm)

B5 <- read.csv("data/corr_B5.csv")

B5 <- B5 |>
  subset(species_mex != "Check")

B5$species_mex <- factor(
  B5$species_mex,
  levels = c("B73", "Bals", "Zdip", "Hueh", "Zlux", "Mex")
)

species_cols <- c(
  "B73"  = "#03bec4",
  "Bals" = "#f364e2",
  "Zdip" = "#f8756d",
  "Hueh" = "#b69d00",
  "Zlux" = "#00b837",
  "Mex"  = "#609bfe"
)

p <- ggplot(
  B5,
  aes(x = species_mex, y = Predicted_N_corr, color = species_mex)
) +
  geom_quasirandom(
    width = 0.25,
    alpha = 0.6,
    size = 2
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    color = "black",
    linewidth = 1,
    size = 0.8
  ) +
  scale_color_manual(values = species_cols) +
  labs(
    x = "Species",
    y = "Predicted N"
  ) +
  theme_classic(base_size = 17) +
  theme(
    legend.position = "none"
  )

p

ggsave(
  filename = "output/B5_predictedN_bees_meanSD.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)
