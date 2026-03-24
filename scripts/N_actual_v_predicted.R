# ============================ setup ===========================================

library(tidyverse)

dir.create("output", showWarnings = FALSE)

# ============================ load data =======================================

df <- read.csv("data/2025_curve_actual_v_predicted.csv")

# ============================ calculate R2 ====================================

fit <- lm(Predicted ~ Actual, data = df)
r2  <- summary(fit)$r.squared

# ============================ plot ============================================

p <- ggplot(df, aes(x = Actual, y = Predicted)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate(
    "text",
    x = max(df$Actual),
    y = min(df$Predicted),
    label = paste0("R² = ", round(r2, 3)),
    hjust = 1,
    vjust = 0,
    size = 7
  ) +
  labs(
    x = "Actual N",
    y = "Predicted N",
    title = ""
  ) +
  theme_minimal(base_size = 18)

# print plot
p

# save plot
ggsave(
  filename = file.path("output/octexp_analysis", "2025_actual_vs_predicted.png"),
  plot = p,
  width = 5,
  height = 5,
  dpi = 300
)
