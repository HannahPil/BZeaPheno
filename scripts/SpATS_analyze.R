# ==============================================================================
# SpATS SPATIAL ANALYSIS
# Fit spatial models (via statgenHTP / SpATS) to field trial phenotypes for the
# B5 (2025) and D4 (2023) experiments, write per-trait corrected value tables
# and spatial plots, and build combined wide tables of corrected values.
# ==============================================================================

library(SpATS)
library(statgenHTP)
library(dplyr)


# ==============================================================================
# B5 EXPERIMENT (2025)
# ==============================================================================

# ------------------------------- load manifest --------------------------------

manifest_B5 <- read.csv("data/B5_spats_manifest.csv")

manifest_B5 <- manifest_B5 |>
  transform(
    Plot_id  = factor(Plot_id),
    Genotype = factor(Genotype),
    Rep      = factor(Rep),
    Row      = as.numeric(Row),
    Range    = as.numeric(Range),
    DTS      = as.numeric(DTS),
    DTA      = as.numeric(DTA),
    GDDTS    = as.numeric(GDDTS),
    GDDTA    = as.numeric(GDDTA),
    PH       = as.numeric(PH),
    EN       = as.numeric(EN),
    EH       = as.numeric(EH),
    NBR      = as.numeric(NBR),
    SPAD     = as.numeric(SPAD),
    LAE      = as.numeric(LAE),
    StPi     = as.numeric(StPi),
    StPu     = as.numeric(StPu)
  )


# ------------------------------- build timepoints -----------------------------

manifest_B5$TP <- as.Date("2024-07-01")  # arbitrary date (statgenHTP requires one)

phenoTP_B5 <- createTimePoints(
  dat            = manifest_B5,
  experimentName = "B5_experiment",
  genotype       = "Genotype",
  timePoint      = "TP",
  repId          = "Rep",
  plotId         = "Plot_id",
  rowNum         = "Range",
  colNum         = "Row",
  addCheck       = TRUE,
  checkGenotypes = "B73",
  timeFormat     = "%Y-%m-%d"
)


# ------------------------------- fit all traits -------------------------------
# for each trait: fit SpATS, save per-trait corrected CSV, save spatial plot,
# and append *_corr column to a combined wide table.

traits_B5 <- c("DTS","DTA","GDDTS","GDDTA","PH","EN","EH","NBR","SPAD","LAE",
               "StPi","StPu")

B5_corr_wide <- manifest_B5 |>
  dplyr::select(Plot_id) |>
  dplyr::distinct() |>
  dplyr::arrange(Plot_id)

for (tr in traits_B5) {

  fit <- fitModels(
    TP     = phenoTP_B5,
    trait  = tr,
    engine = "SpATS",
    what   = "fixed"
  )

  # per-trait corrected table
  corr <- getCorrected(fit) |>
    dplyr::rename(Plot_id = plotId) |>
    dplyr::left_join(
      manifest_B5[, c("Plot_id","Genotype","Rep","Row","Range")],
      by = "Plot_id"
    ) |>
    dplyr::arrange(Plot_id)

  write.csv(corr, paste0("data/B5_", tr, "_2025.csv"), row.names = FALSE)

  # spatial plot
  png(
    filename = file.path("output", "spatial",
                         paste0("B5_", tr, "_spatial_raw_2025.png")),
    width    = 1800,
    height   = 1400,
    res      = 200
  )
  plot(fit, timePoints = "2024-07-01", plotType = "spatial", spaTrend = "raw")
  dev.off()

  # append to wide table
  corr_col <- paste0(tr, "_corr")
  B5_corr_wide <- B5_corr_wide |>
    dplyr::left_join(corr[, c("Plot_id", corr_col)], by = "Plot_id") |>
    dplyr::arrange(Plot_id)
}

write.csv(B5_corr_wide, "data/B5_all_traits_corr_2025.csv", row.names = FALSE)


# ==============================================================================
# B5 PREDICTED N (block 2)
# ==============================================================================

b5_N <- read.csv("data/B5_block2_predictedN.csv")

b5_N <- b5_N |>
  transform(
    Plot        = factor(Plot),
    Genotype    = factor(Genotype),
    Predicted_N = as.numeric(Predicted_N),
    Row         = as.numeric(Row),
    Range       = as.numeric(Range)
  )

# plot-level averages; only keep Predicted_N rows with M_Distance <= 4
b5_plot_avg <- b5_N |>
  group_by(Plot) |>
  summarise(
    Row         = first(Row),
    Range       = first(Range),
    Genotype    = first(Genotype),
    Accession   = first(Accession),
    Taxa        = first(Taxa),
    M_Distance  = mean(M_Distance),
    Predicted_N = mean(Predicted_N[M_Distance <= 4]),
    M_Distance  = mean(M_Distance[M_Distance <= 4]),
    .groups     = "drop"
  )

b5_plot_avg$TP <- as.Date("2024-07-01")

pheno_b5_N <- createTimePoints(
  dat            = b5_plot_avg,
  experimentName = "B5_experiment",
  genotype       = "Genotype",
  timePoint      = "TP",
  plotId         = "Plot",
  rowNum         = "Range",
  colNum         = "Row",
  addCheck       = TRUE,
  checkGenotypes = "B73",
  timeFormat     = "%Y-%m-%d"
)

fit_B5_N <- fitModels(
  TP     = pheno_b5_N,
  trait  = "Predicted_N",
  engine = "SpATS",
  what   = "fixed"
)

# spatial plot
png(
  filename = file.path("output", "spatial",
                       "B5_block2_PredictedN_spatial_raw_BLUPs_2025.png"),
  width    = 1800,
  height   = 1400,
  res      = 200
)
plot(fit_B5_N, timePoints = "2024-07-01", plotType = "spatial", spaTrend = "raw")
dev.off()

# corrected values
corr_B5_N <- getCorrected(fit_B5_N) |>
  rename(Plot = plotId) |>
  left_join(
    b5_plot_avg |> dplyr::select(Plot, Genotype, Accession, Taxa, Row, Range),
    by = "Plot"
  )

write.csv(corr_B5_N, "data/B5_block2_PredictedN_corrected.csv", row.names = FALSE)


# ==============================================================================
# D4 EXPERIMENT (2023)
# ==============================================================================

# ------------------------------- load manifest --------------------------------

manifest_D4 <- read.csv("data/D4_spats_manifest.csv")

manifest_D4 <- manifest_D4 |>
  transform(
    Plot_id  = factor(Plot_id),
    Genotype = factor(Genotype),
    Rep      = factor(Rep),
    Row      = as.numeric(Row),
    Range    = as.numeric(Range),
    GDDTS    = as.numeric(GDDTS),
    GDDTA    = as.numeric(GDDTA),
    PH       = as.numeric(PH),
    EN       = as.numeric(EN),
    EH       = as.numeric(EH),
    NBR      = as.numeric(NBR),
    SPAD1    = as.numeric(SPAD1),
    SPAD2    = as.numeric(SPAD2),
    StPu     = as.numeric(StPu),
    StPi     = as.numeric(StPi),
    BW       = as.numeric(BW),
    BL       = as.numeric(BL),
    SL       = as.numeric(SL)
  )


# ------------------------------- build timepoints -----------------------------

manifest_D4$TP <- as.Date("2024-07-01")

phenoTP_D4 <- createTimePoints(
  dat            = manifest_D4,
  experimentName = "D4_experiment",
  genotype       = "Genotype",
  timePoint      = "TP",
  repId          = "Rep",
  plotId         = "Plot_id",
  rowNum         = "Range",
  colNum         = "Row",
  addCheck       = TRUE,
  checkGenotypes = c("B73"),
  timeFormat     = "%Y-%m-%d"
)


# ------------------------------- fit all traits -------------------------------
# for each trait: fit SpATS, save per-trait corrected CSV, save spatial plot,
# and append *_corr column to a combined wide table.

traits_D4 <- c("GDDTS","GDDTA","PH","EH","NBR","SPAD1","SPAD2",
               "StPu","StPi","BW","BL","SL")

D4_corr_wide <- manifest_D4 |>
  dplyr::select(Plot_id) |>
  dplyr::distinct() |>
  dplyr::arrange(Plot_id)

for (tr in traits_D4) {

  fit <- fitModels(
    TP     = phenoTP_D4,
    trait  = tr,
    engine = "SpATS",
    what   = "fixed"
  )

  # per-trait corrected table
  corr <- getCorrected(fit) |>
    dplyr::rename(Plot_id = plotId) |>
    dplyr::left_join(
      manifest_D4[, c("Plot_id","Genotype","Rep","Row","Range")],
      by = "Plot_id"
    ) |>
    dplyr::arrange(Plot_id)

  write.csv(corr, paste0("data/D4_", tr, "_2023.csv"), row.names = FALSE)

  # spatial plot
  png(
    filename = file.path("output", "spatial",
                         paste0("D4_", tr, "_spatial_raw_2023.png")),
    width    = 1800,
    height   = 1400,
    res      = 200
  )
  plot(fit, timePoints = "2024-07-01", plotType = "spatial", spaTrend = "raw")
  dev.off()

  # append to wide table
  corr_col <- paste0(tr, "_corr")
  D4_corr_wide <- D4_corr_wide |>
    dplyr::left_join(corr[, c("Plot_id", corr_col)], by = "Plot_id") |>
    dplyr::arrange(Plot_id)
}

write.csv(D4_corr_wide, "data/D4_all_traits_corr_2023.csv", row.names = FALSE)
