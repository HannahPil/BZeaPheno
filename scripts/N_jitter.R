library(ggplot2)
library(dplyr)

df <- read.csv("BZea_Block2_N-predictions.csv")

# keep only inbreds (this excludes B73xLH287)
df_inbred <- df %>%
  filter(inbred_hybrid == "inbred")

# reorder taxa: B73 first, then all remaining taxa alphabetically
df_inbred <- df_inbred %>%
  mutate(
    Taxa = factor(
      Taxa,
      levels = c(
        "B73",
        sort(setdiff(unique(Taxa), "B73"))
      )
    )
  )

# taxa-level stats (mean ± sd) for Predicted_N
group_means <- df_inbred %>%
  group_by(Taxa) %>%
  summarise(
    mean_N = mean(Predicted_N, na.rm = TRUE),
    sd_N   = sd(Predicted_N, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_inbred, aes(x = Taxa, y = Predicted_N, fill = Taxa)) +
  geom_jitter(
    shape = 21, size = 3, alpha = 0.8,
    position = position_jitter(width = 0.1)
  ) +
  geom_errorbar(
    data = group_means,
    aes(
      x = Taxa,
      ymin = mean_N - sd_N,
      ymax = mean_N + sd_N
    ),
    width = 0,
    linewidth = 1,      # replaces size for lines in newer ggplot2
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_means,
    aes(x = Taxa, y = mean_N),
    shape = 16, size = 2,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Taxa",
    y = "Predicted N",
    title = "Predicted N by Taxa (inbreds only)"
  ) +
  theme_minimal()


#==============================================
#averaged scans
#==============================================

# load packages
library(ggplot2)
library(dplyr)

# read data
df <- read.csv("BZea_Block2_N-predictions.csv")

# keep only inbreds and set taxa order (b73 first)
df_inbred <- df %>%
  filter(inbred_hybrid == "inbred") %>%
  mutate(
    Taxa = factor(
      Taxa,
      levels = c(
        "B73",
        sort(setdiff(unique(Taxa), "B73"))
      )
    )
  )

# average predicted n by plot (plot-level means)
df_plotmeans <- df_inbred %>%
  group_by(Plot, Taxa) %>%           # keep taxa so each plot stays linked
  summarise(
    Predicted_N = mean(Predicted_N, na.rm = TRUE),
    .groups = "drop"
  )

# taxa-level stats (mean ± sd) using plot means
group_means <- df_plotmeans %>%
  group_by(Taxa) %>%
  summarise(
    mean_N = mean(Predicted_N, na.rm = TRUE),
    sd_N   = sd(Predicted_N, na.rm = TRUE),
    .groups = "drop"
  )

# jitter plot using plot means
ggplot(df_plotmeans, aes(x = Taxa, y = Predicted_N, fill = Taxa)) +
  geom_jitter(
    shape = 21, size = 3, alpha = 0.8,
    position = position_jitter(width = 0.1)
  ) +
  geom_errorbar(
    data = group_means,
    aes(
      x = Taxa,
      ymin = mean_N - sd_N,
      ymax = mean_N + sd_N
    ),
    width = 0,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_means,
    aes(x = Taxa, y = mean_N),
    shape = 16, size = 2,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Taxa",
    y = "Predicted N",
    title = "Predicted N by Taxa (inbreds only, plot means)"
  ) +
  theme_minimal()


#--------------spatially corrected scans-----------------------------
df_spat <- read.csv("B5_block2_PredictedN_corrected.csv")


# jitter plot using plot means
ggplot(df_plotmeans, aes(x = Taxa, y = Predicted_N, fill = Taxa)) +
  geom_jitter(
    shape = 21, size = 3, alpha = 0.8,
    position = position_jitter(width = 0.1)
  ) +
  geom_errorbar(
    data = group_means,
    aes(
      x = Taxa,
      ymin = mean_N - sd_N,
      ymax = mean_N + sd_N
    ),
    width = 0,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_means,
    aes(x = Taxa, y = mean_N),
    shape = 16, size = 2,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Taxa",
    y = "Predicted N",
    title = "Predicted N by Taxa (inbreds only, plot means)"
  ) +
  theme_minimal()