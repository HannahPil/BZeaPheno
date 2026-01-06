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


#--------------spatially corrected scans (calibration curve set)-----------------------------
library(ggplot2)
library(dplyr)

df_spat <- read.csv("B5_block2_PredictedN_corrected.csv")

df_spat <- df_spat %>%
  filter(!is.na(Taxa)) %>%              # <-- drop NA taxa
  mutate(
    Taxa = factor(
      Taxa,
      levels = c(
        "B73",
        sort(setdiff(unique(Taxa), "B73"))
      )
    )
  )

# taxa-level stats
group_means <- df_spat %>%
  group_by(Taxa) %>%
  summarise(
    mean_N = mean(Predicted_N_corr, na.rm = TRUE),
    sd_N   = sd(Predicted_N_corr, na.rm = TRUE),
    .groups = "drop"
  )

# jitter plot
ggplot(df_spat, aes(x = Taxa, y = Predicted_N_corr, fill = Taxa)) +
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
    y = "Predicted N (corrected)",
    title = "Predicted N by Taxa (corrected plot-level values)"
  ) +
  theme_minimal()

#---------------averaged block 2 scans------------------------
library(ggplot2)
library(dplyr)

pred_raw <- read.csv("B5_block2_predictedN.csv")

pred_clean <- pred_raw %>%
  filter(
    Taxa != "#N/A",
    Genotype != "#N/A",
    Predicted_N != "#N/A"
  ) %>%
  mutate(
    # convert Predicted_N to numeric if it was read as character
    Predicted_N = as.numeric(Predicted_N),
    Taxa = factor(
      Taxa,
      levels = c("B73", sort(setdiff(unique(Taxa), "B73")))
    )
  )


# average Predicted_N by plot (keeping genotype + taxa)
pred_plot_means <- pred_clean %>%
  group_by(Plot, Genotype, Taxa) %>%
  summarise(
    mean_PredN = mean(Predicted_N, na.rm = TRUE),
    .groups = "drop"
  )

# taxa-level stats from plot means
pred_taxa_stats <- pred_plot_means %>%
  group_by(Taxa) %>%
  summarise(
    taxa_mean = mean(mean_PredN, na.rm = TRUE),
    taxa_sd   = sd(mean_PredN, na.rm = TRUE),
    .groups = "drop"
  )

# plot
ggplot(pred_plot_means, aes(x = Taxa, y = mean_PredN, fill = Taxa)) +
  geom_jitter(
    shape = 21, size = 3, alpha = 0.8,
    position = position_jitter(width = 0.1)
  ) +
  geom_point(
    data = pred_taxa_stats,
    aes(x = Taxa, y = taxa_mean),
    shape = 16, size = 2,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(2, 5)) +
  labs(
    x = "Taxa",
    y = "Predicted N",
    title = "Predicted N by Taxa (averaged by plot) — B5 Block 2"
  ) +
  theme_minimal()


#-------spatially corrected and averaged scans----------------
library(ggplot2)
library(dplyr)

df_spat <- read.csv("B5_block2_PredictedN_corrected.csv")

# drop NA taxa, reorder
df_spat <- df_spat %>%
  filter(!is.na(Taxa)) %>%
  mutate(
    Taxa = factor(
      Taxa,
      levels = c(
        "B73",
        sort(setdiff(unique(Taxa), "B73"))
      )
    )
  )

# average by plot
plot_means <- df_spat %>%
  group_by(Plot, Genotype, Taxa) %>%
  summarise(
    mean_N = mean(Predicted_N_corr, na.rm = TRUE),
    .groups = "drop"
  )

# taxa-level stats from plot means
group_means <- plot_means %>%
  group_by(Taxa) %>%
  summarise(
    mean_N = mean(mean_N, na.rm = TRUE),
    sd_N   = sd(mean_N, na.rm = TRUE),
    .groups = "drop"
  )

# plot plot-means
ggplot(plot_means, aes(x = Taxa, y = mean_N, fill = Taxa)) +
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
  scale_y_continuous(limits = c(2, 5)) +
  labs(
    x = "Taxa",
    y = "Predicted N (corrected)",
    title = "Predicted N by Taxa"
  ) +
  theme_minimal()


#----correlation between SPAD and predicted_N----------------
spad <- read.csv("corr_B5.csv")
pred_N <- read.csv("B5_block2_PredictedN_corrected.csv")

df <- spad %>%
  select(Plot, SPAD_corr) %>%
  inner_join(
    pred_N %>% select(Plot, Predicted_N_corr),
    by = "Plot"
  )

cor(df$SPAD, df$Predicted_N_corr, use = "complete.obs")

# compute correlation
r_val <- cor(df$SPAD_corr, df$Predicted_N_corr, use = "complete.obs")
r2_val <- r_val^2
label_text <- paste0("r² = ", round(r2_val, 3))

# plot with correlation annotation
ggplot(df, aes(x = SPAD_corr, y = Predicted_N_corr)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text", x = Inf, y = Inf, label = label_text,
           hjust = 1.7, vjust = 3.4, size = 7) +
  theme_minimal() +
  labs(
    x = "SPAD (corrected)",
    y = "Predicted N (corrected)",
    title = "Correlation Between SPAD and Predicted_N_corr"
  )
