# ==============================================================================
# DRONE vs FIELDBOOK HEIGHT CORRELATION
# Quick scatter plot: drone height on 8/24/2023 vs hand-measured PH from D4 fieldbook
# ==============================================================================
library(tidyverse)

# ------------------------------- output dir -----------------------------------
out_dir <- file.path("output", "drone_vs_fieldbook")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- load data ------------------------------------
data_dir <- "data"

timeseries <- read.csv(file.path(data_dir, "maize_height_timeseries.csv"),
                       stringsAsFactors = FALSE)

fieldbook  <- read.csv(file.path(data_dir, "UPDATED_CLY23_D4_FieldBook.csv"),
                       stringsAsFactors = FALSE)

# ------------------------------- filter & join --------------------------------

# drone height on 8/24/2023
drone_824 <- timeseries |>
  filter(flight_date == "8/24/2023") |>
  transmute(
    plot       = as.integer(plot),
    drone_cm   = height_cm
  )

# fieldbook hand-measured PH (cm)
fb_ph <- fieldbook |>
  transmute(
    plot     = as.integer(CLY23_D4),
    fb_PH    = as.numeric(PH),
    Genotype = genotype_id,
    species  = species_mex
  ) |>
  filter(!is.na(fb_PH))

joined <- inner_join(drone_824, fb_ph, by = "plot")

cat("Matched plots:", nrow(joined), "\n")

# ------------------------------- correlation ----------------------------------
cor_val <- cor(joined$drone_cm, joined$fb_PH, use = "complete.obs")
cat("Pearson r:", round(cor_val, 3), "\n")

# ------------------------------- plot -----------------------------------------
p <- ggplot(joined, aes(x = fb_PH, y = drone_cm)) +
  geom_point(
    aes(fill = species),
    shape = 21,
    size  = 2.5,
    alpha = 0.6,
    color = "grey30"
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  annotate(
    "text",
    x    = min(joined$fb_PH, na.rm = TRUE) + 5,
    y    = max(joined$drone_cm, na.rm = TRUE) - 5,
    label = paste0("r = ", round(cor_val, 3)),
    hjust = 0,
    size  = 5,
    fontface = "bold"
  ) +
  labs(
    x     = "Fieldbook plant height (cm)",
    y     = "Drone height 8/24/2023 (cm)",
    title = "Drone vs fieldbook plant height",
    fill  = "taxa"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

print(p)

ggsave(
  file.path(out_dir, "drone_824_vs_fieldbook_PH.png"),
  p,
  width  = 6,
  height = 5
)
