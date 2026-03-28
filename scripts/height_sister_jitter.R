# ==============================================================================
# SISTER-LINE HEIGHT COMPARISON AT FOCUS GENE
# Compare plant height (from drone timeseries) of teosinte-introgression carriers
# vs their BC1 sister lines at a specific gene locus.
# Required files in data/:
#   maize_height_timeseries.csv
#   D4_spats_manifest.csv
#   corr_D4.csv
#   results_list_new_name.rds
#   RNAmetadata.csv
#   edgeR_log2cpm_TMM_filtered.csv
# ==============================================================================

library(tidyverse)

# ##############################################################################
# ##                                                                          ##
# ##   >>> CHANGE THIS GENE ID / COORDINATES TO ANALYZE A DIFFERENT GENE <<<  ##
# ##                                                                          ##
# ##############################################################################

FOCUS_GENE <- "Zm00001eb012750"
FOCUS_CHR  <- "chr1"
FOCUS_START <- 42029359
FOCUS_END   <- 42034183

# flight date to use for height comparison
HEIGHT_DATE <- "5/31/2023"

# ------------------------------- output dir -----------------------------------

out_dir <- file.path("output", "height_sister_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- load data ------------------------------------

data_dir <- "data"

heights  <- read.csv(file.path(data_dir, "maize_height_timeseries.csv"))
manifest <- read.csv(file.path(data_dir, "D4_spats_manifest.csv"))
corr_D4  <- read.csv(file.path(data_dir, "corr_D4.csv"))
teogeno  <- readRDS(file.path(data_dir, "results_list_new_name.rds"))
teogeno  <- teogeno[!duplicated(names(teogeno))]

# ------------------------------- height at target date ------------------------

height_at_date <- heights |>
  filter(flight_date == HEIGHT_DATE) |>
  select(plot, height_cm)

# ------------------------------- map plots to genotypes -----------------------

plot_geno <- manifest |>
  select(Plot_id, Genotype) |>
  distinct()

# add species from corr_D4
species_lookup <- corr_D4 |>
  select(Genotype, Species, species_mex) |>
  distinct()

plot_info <- height_at_date |>
  inner_join(plot_geno, by = c("plot" = "Plot_id")) |>
  left_join(species_lookup, by = "Genotype")

# ------------------------------- introgression status at focus gene ------------

# for each genotype, check if teogeno has an introgression overlapping the gene
get_intro_status <- function(genotype, teogeno, chr, start, end) {
  key <- paste0(genotype, ".B")
  if (!key %in% names(teogeno)) return(FALSE)

  segs <- teogeno[[key]]
  intro <- segs |>
    filter(V4 == "Introgression", V1 == chr)

  if (nrow(intro) == 0) return(FALSE)

  any(intro$V2 <= end & intro$V3 >= start)
}

unique_genos <- unique(plot_info$Genotype)

intro_status <- tibble(
  Genotype = unique_genos,
  has_teo = map_lgl(unique_genos, ~get_intro_status(
    .x, teogeno, FOCUS_CHR, FOCUS_START, FOCUS_END
  ))
)

plot_info <- plot_info |>
  left_join(intro_status, by = "Genotype")

# ------------------------------- BC1 families ---------------------------------

# exclude checks
plot_info <- plot_info |>
  filter(!Genotype %in% c("B73", "Purple Check"))

# BC1 family = first 10 characters of genotype name
plot_info <- plot_info |>
  mutate(
    bc1_family = substr(Genotype, 1, 10),
    allele = ifelse(has_teo, "Teosinte", "B73")
  )

# keep only families that have BOTH teo and B73 alleles
families_with_both <- plot_info |>
  group_by(bc1_family) |>
  filter(n_distinct(allele) == 2) |>
  ungroup()

cat("\n=== Lines carrying teosinte introgression at", FOCUS_GENE, "===\n")
teo_carriers <- plot_info |> filter(has_teo) |> select(Genotype, Species, bc1_family) |> distinct()
print(as.data.frame(teo_carriers), row.names = FALSE)

cat("\nFamilies with teo + sister pairs:", length(unique(families_with_both$bc1_family)), "\n")

# ------------------------------- colors ---------------------------------------

taxa_colors <- c(
  "B73"  = "#03bec4",
  "Bals" = "#f364e2",
  "Zdip" = "#f8756d",
  "Hueh" = "#b69d00",
  "Zlux" = "#00b837",
  "Dura" = "#609bfe",
  "Nobo" = "#609bfe",
  "Mesa" = "#609bfe",
  "Chal" = "#609bfe",
  "Mex"  = "#609bfe",
  "Check" = "grey60"
)

# ------------------------------- effect sizes ---------------------------------

family_means <- families_with_both |>
  group_by(bc1_family, allele) |>
  summarise(
    mean_ht = mean(height_cm, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = allele,
    values_from = c(mean_ht, n)
  ) |>
  filter(!is.na(mean_ht_Teosinte), !is.na(mean_ht_B73)) |>
  mutate(effect = mean_ht_Teosinte - mean_ht_B73) |>
  left_join(
    families_with_both |>
      filter(has_teo) |>
      select(bc1_family, species_mex) |>
      distinct(),
    by = "bc1_family"
  )

cat("\n=== Effect sizes (teo - B73) within BC1 families ===\n")
print(as.data.frame(family_means |> select(bc1_family, species_mex, effect, n_Teosinte, n_B73)), row.names = FALSE)

write.csv(
  family_means,
  file.path(out_dir, paste0(FOCUS_GENE, "_height_sister_effect_sizes.csv")),
  row.names = FALSE
)

# ------------------------------- waterfall bar chart ---------------------------

family_means_sorted <- family_means |>
  filter(!is.na(species_mex)) |>
  arrange(effect) |>
  mutate(bc1_family = factor(bc1_family, levels = bc1_family))

p_waterfall <- ggplot(family_means_sorted, aes(x = bc1_family, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(aes(fill = species_mex), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = taxa_colors, name = "Taxa") +
  labs(
    x = "BC1 family",
    y = "Effect size (Teo - B73, cm)",
    title = paste("Introgression effect on early height by family:", FOCUS_GENE),
    subtitle = paste0("Height on ", HEIGHT_DATE, "; positive = teosinte allele increases height")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title  = element_text(face = "bold")
  )

print(p_waterfall)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_height_waterfall.png")),
  p_waterfall,
  width = max(6, nrow(family_means_sorted) * 0.4 + 2),
  height = 5
)

# ------------------------------- paired waterfall (height + expression) --------

# load RNA-seq data
rna_meta <- read.csv(file.path(data_dir, "RNAmetadata.csv"), stringsAsFactors = FALSE)
tmm      <- read.csv(file.path(data_dir, "edgeR_log2cpm_TMM_filtered.csv"),
                      check.names = FALSE, row.names = 1)

# extract focus gene expression per sample
focus_expr <- as.numeric(tmm[FOCUS_GENE, ])
names(focus_expr) <- colnames(tmm)

# build sample-level table: sample → genotype → bc1_family → introgression status
rna_samples <- rna_meta |>
  filter(
    sample_id %in% colnames(tmm),
    !genotype %in% c("B73", "Purple Check", "NA")
  ) |>
  mutate(
    bc1_family = substr(genotype, 1, 10),
    expression = focus_expr[sample_id]
  ) |>
  filter(!is.na(expression))

# determine introgression status at focus gene for each sequenced genotype
rna_genos <- unique(rna_samples$genotype)

rna_intro <- tibble(
  genotype = rna_genos,
  has_teo_rna = map_lgl(rna_genos, ~get_intro_status(
    .x, teogeno, FOCUS_CHR, FOCUS_START, FOCUS_END
  ))
)

rna_samples <- rna_samples |>
  left_join(rna_intro, by = "genotype") |>
  mutate(allele = ifelse(has_teo_rna, "Teosinte", "B73"))

# compute expression effect per BC1 family (teo mean - B73 mean)
# keep only families with both alleles represented in RNA data
expr_family_means <- rna_samples |>
  group_by(bc1_family, allele) |>
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    n_expr = n(),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from  = allele,
    values_from = c(mean_expr, n_expr)
  ) |>
  filter(!is.na(mean_expr_Teosinte), !is.na(mean_expr_B73)) |>
  mutate(expr_effect = mean_expr_Teosinte - mean_expr_B73)

cat("\n=== Expression effect sizes (families with RNA data for both alleles) ===\n")
print(as.data.frame(expr_family_means), row.names = FALSE)

# merge height and expression effects for families that have both
paired_effects <- family_means_sorted |>
  select(bc1_family, species_mex, effect) |>
  rename(ht_effect = effect) |>
  left_join(
    expr_family_means |> select(bc1_family, expr_effect),
    by = "bc1_family"
  )

# pivot long for faceted panel chart
paired_long <- paired_effects |>
  filter(!is.na(expr_effect)) |>
  pivot_longer(
    cols = c(ht_effect, expr_effect),
    names_to = "measure",
    values_to = "effect"
  ) |>
  mutate(
    measure = recode(
      measure,
      ht_effect   = "Height effect (cm)",
      expr_effect = "Expression effect (log2 CPM)"
    ),
    measure = factor(measure, levels = c("Height effect (cm)", "Expression effect (log2 CPM)"))
  )

# reorder families by height effect
fam_order <- paired_effects |>
  filter(!is.na(expr_effect)) |>
  arrange(ht_effect) |>
  pull(bc1_family) |>
  as.character()

paired_long <- paired_long |>
  mutate(bc1_family = factor(bc1_family, levels = fam_order))

p_paired <- ggplot(paired_long, aes(x = bc1_family, y = effect, fill = species_mex)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(width = 0.7, alpha = 0.85) +
  facet_wrap(~ measure, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = taxa_colors, name = "Taxa") +
  labs(
    x = "BC1 family",
    y = "Effect size (Teo - B73)",
    title = paste("Early height and expression effects by family:", FOCUS_GENE),
    subtitle = paste0("Height on ", HEIGHT_DATE, "; only families with RNA-seq data for both alleles")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title  = element_text(face = "bold"),
    strip.text  = element_text(size = 13, face = "bold")
  )

print(p_paired)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_paired_waterfall.png")),
  p_paired,
  width = max(6, length(fam_order) * 0.5 + 2),
  height = 8
)

# also save the height-only waterfall with a marker for sequenced families
# (families with RNA data get a bold outline)
sequenced_families <- unique(expr_family_means$bc1_family)

family_means_sorted <- family_means_sorted |>
  mutate(sequenced = bc1_family %in% sequenced_families)

p_waterfall_marked <- ggplot(family_means_sorted, aes(x = bc1_family, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(
    aes(fill = species_mex, linewidth = sequenced, color = sequenced),
    width = 0.7, alpha = 0.85
  ) +
  scale_fill_manual(values = taxa_colors, name = "Taxa") +
  scale_linewidth_manual(values = c("FALSE" = 0.2, "TRUE" = 1.2), guide = "none") +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "black"),
    labels = c("FALSE" = "No RNA data", "TRUE" = "RNA-seq confirmed"),
    name = "Expression"
  ) +
  labs(
    x = "BC1 family",
    y = "Effect size (Teo - B73, cm)",
    title = paste("Introgression effect on early height by family:", FOCUS_GENE),
    subtitle = paste0("Height on ", HEIGHT_DATE, "; bold outline = RNA-seq confirmed higher expression")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    plot.title  = element_text(face = "bold")
  )

print(p_waterfall_marked)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_height_waterfall_rna_marked.png")),
  p_waterfall_marked,
  width = max(6, nrow(family_means_sorted) * 0.4 + 2),
  height = 5
)

# ------------------------------- faceted sister jitter plot --------------------

build_facet_jitter <- function(df_fam, fam_levels, show_legend = TRUE,
                               title_suffix = "", nrow = NULL) {
  df_fam <- df_fam |>
    mutate(
      bc1_family = factor(bc1_family, levels = fam_levels),
      allele     = factor(allele, levels = c("B73", "Teosinte"))
    )

  fstats <- df_fam |>
    group_by(bc1_family, allele) |>
    summarise(
      mean_ht = mean(height_cm, na.rm = TRUE),
      sd_ht   = sd(height_cm, na.rm = TRUE),
      n       = n(),
      .groups = "drop"
    ) |>
    mutate(sd_ht = ifelse(is.na(sd_ht), 0, sd_ht))

  p <- ggplot(df_fam, aes(x = allele, y = height_cm)) +
    geom_jitter(
      aes(fill = species_mex),
      shape    = 21,
      size     = 3,
      alpha    = 0.75,
      color    = "grey30",
      position = position_jitter(width = 0.2)
    ) +
    geom_errorbar(
      data = fstats,
      aes(
        x    = allele,
        ymin = mean_ht - sd_ht,
        ymax = mean_ht + sd_ht
      ),
      width     = 0,
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = fstats,
      aes(x = allele, y = mean_ht),
      shape = 16,
      size  = 4,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ bc1_family, nrow = nrow) +
    scale_fill_manual(values = taxa_colors, name = "teosinte taxa") +
    scale_x_discrete(labels = c("B73" = "B73", "Teosinte" = "Teo")) +
    labs(
      x = "genotype",
      y = paste0("plant height on ", HEIGHT_DATE, " (cm)"),
      title = paste0("Per-family early height: ", FOCUS_GENE, title_suffix),
      subtitle = "Each panel = one BC1 family; black = mean ± SD"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold"),
      strip.text  = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 11, face = "bold")
    )

  if (!show_legend) p <- p + theme(legend.position = "none")
  return(p)
}

# sort families by effect size
family_order <- levels(family_means_sorted$bc1_family)

paired_df <- families_with_both |>
  filter(!is.na(species_mex)) |>
  filter(bc1_family %in% family_means_sorted$bc1_family)

# -- all families plot --
p_facet <- build_facet_jitter(paired_df, family_order)
print(p_facet)

n_fam <- length(unique(paired_df$bc1_family))

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_height_faceted_families.png")),
  p_facet,
  width  = min(16, max(6, ceiling(sqrt(n_fam)) * 3)),
  height = min(14, max(4, ceiling(n_fam / ceiling(sqrt(n_fam))) * 3))
)

# -- top 5 families by absolute effect size (no legend) --
top5_families <- family_means_sorted |>
  arrange(desc(abs(effect))) |>
  slice_head(n = 5) |>
  arrange(effect) |>
  pull(bc1_family) |>
  as.character()

paired_df_top5 <- paired_df |>
  filter(bc1_family %in% top5_families)

p_top5 <- build_facet_jitter(
  paired_df_top5,
  top5_families,
  show_legend  = FALSE,
  title_suffix = " (top 5)",
  nrow         = 1
)
print(p_top5)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_height_faceted_top5.png")),
  p_top5,
  width  = 16,
  height = 4
)

# ------------------------------- pooled teo vs B73 jitter ---------------------

# all lines with introgression at focus gene vs all without (not just sisters)
pooled_df <- plot_info |>
  filter(!is.na(species_mex)) |>
  mutate(
    Genotype_group = factor(
      ifelse(has_teo, "Teosinte Introgression", "B73 Background"),
      levels = c("B73 Background", "Teosinte Introgression")
    ),
    point_color = ifelse(has_teo, as.character(species_mex), "B73")
  )

group_stats <- pooled_df |>
  group_by(Genotype_group) |>
  summarise(
    mean_ht = mean(height_cm, na.rm = TRUE),
    sd_ht   = sd(height_cm, na.rm = TRUE),
    .groups = "drop"
  )

p_pooled <- ggplot(pooled_df, aes(x = Genotype_group, y = height_cm)) +
  geom_jitter(
    aes(fill = point_color),
    shape = 21,
    size = 3,
    alpha = 0.75,
    color = "grey30",
    position = position_jitter(width = 0.2)
  ) +
  geom_errorbar(
    data = group_stats,
    aes(
      x = Genotype_group,
      ymin = mean_ht - sd_ht,
      ymax = mean_ht + sd_ht
    ),
    width = 0,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_stats,
    aes(x = Genotype_group, y = mean_ht),
    shape = 16,
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_x_discrete(labels = c(
    "B73 Background" = "B73",
    "Teosinte Introgression" = "Teo"
  )) +
  scale_fill_manual(
    values = taxa_colors,
    name = "teosinte taxa"
  ) +
  labs(
    x = "genotype",
    y = paste0("plant height on ", HEIGHT_DATE, " (cm)"),
    title = paste("Effect of introgression on early height:", FOCUS_GENE),
    subtitle = "Teo points colored by taxa; black = mean ± SD"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    plot.title  = element_text(face = "bold")
  )

print(p_pooled)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_height_pooled.png")),
  p_pooled,
  width = 5,
  height = 5
)

# -- same plot but filtered to sequenced lines only ----------------------------

sequenced_genos <- unique(rna_meta$genotype[rna_meta$sample_id %in% colnames(tmm)])

pooled_seq_df <- plot_info |>
  filter(!is.na(species_mex), Genotype %in% sequenced_genos) |>
  mutate(
    Genotype_group = factor(
      ifelse(has_teo, "Teosinte Introgression", "B73 Background"),
      levels = c("B73 Background", "Teosinte Introgression")
    ),
    point_color = ifelse(has_teo, as.character(species_mex), "B73")
  )

group_stats_seq <- pooled_seq_df |>
  group_by(Genotype_group) |>
  summarise(
    mean_ht = mean(height_cm, na.rm = TRUE),
    sd_ht   = sd(height_cm, na.rm = TRUE),
    .groups = "drop"
  )

p_pooled_seq <- ggplot(pooled_seq_df, aes(x = Genotype_group, y = height_cm)) +
  geom_jitter(
    aes(fill = point_color),
    shape = 21,
    size = 3,
    alpha = 0.75,
    color = "grey30",
    position = position_jitter(width = 0.2)
  ) +
  geom_errorbar(
    data = group_stats_seq,
    aes(
      x = Genotype_group,
      ymin = mean_ht - sd_ht,
      ymax = mean_ht + sd_ht
    ),
    width = 0,
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = group_stats_seq,
    aes(x = Genotype_group, y = mean_ht),
    shape = 16,
    size = 4,
    color = "black",
    inherit.aes = FALSE
  ) +
  scale_x_discrete(labels = c(
    "B73 Background" = "B73",
    "Teosinte Introgression" = "Teo"
  )) +
  scale_fill_manual(
    values = taxa_colors,
    name = "teosinte taxa"
  ) +
  labs(
    x = "genotype",
    y = paste0("plant height on ", HEIGHT_DATE, " (cm)"),
    title = paste0("Introgression effect on height at ", FOCUS_GENE),
    subtitle = paste0("Sequenced lines only; height on ", HEIGHT_DATE, "; black = mean ± SD")
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    plot.title  = element_text(face = "bold")
  )

print(p_pooled_seq)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_height_pooled_sequenced_only.png")),
  p_pooled_seq,
  width = 5,
  height = 5
)

print(p_pooled)

ggsave(
  file.path(out_dir, paste0(FOCUS_GENE, "_height_pooled.png")),
  p_pooled,
  width = 5,
  height = 5
)

# ##############################################################################
# ##  HEIGHT vs EXPRESSION CORRELATION PLOTS                                  ##
# ##  For all sequenced lines: expression at FOCUS_GENE (x) vs height (y)     ##
# ##############################################################################

corr_dir <- file.path("output", "height_expression_correlations")
dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)

# genotype-level mean expression at FOCUS_GENE (average across RNA-seq replicates)
geno_expr <- rna_meta |>
  filter(
    sample_id %in% colnames(tmm),
    !genotype %in% c("Purple Check", "NA"),
    !is.na(genotype)
  ) |>
  mutate(expression = focus_expr[sample_id]) |>
  filter(!is.na(expression)) |>
  group_by(genotype, taxa) |>
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

# genotype-level mean height (average across field reps/plots)
geno_height <- heights |>
  filter(flight_date == HEIGHT_DATE) |>
  inner_join(plot_geno, by = c("plot" = "Plot_id")) |>
  filter(!Genotype %in% c("Purple Check")) |>
  group_by(Genotype) |>
  summarise(
    mean_ht = mean(height_cm, na.rm = TRUE),
    n_plots = n(),
    .groups = "drop"
  )

# merge: genotypes that have BOTH expression and height data
corr_df <- geno_expr |>
  inner_join(geno_height, by = c("genotype" = "Genotype")) |>
  left_join(intro_status, by = c("genotype" = "Genotype")) |>
  mutate(
    has_teo = replace_na(has_teo, FALSE),
    allele = ifelse(has_teo, "Teosinte", "B73"),
    allele = factor(allele, levels = c("B73", "Teosinte"))
  )

cat("\n=== Height vs expression correlation data ===\n")
cat("Genotypes with both height and expression:", nrow(corr_df), "\n")
cat("  Teo carriers:", sum(corr_df$has_teo), "\n")
cat("  B73 background:", sum(!corr_df$has_teo), "\n")

# correlation stats
r_all <- cor(corr_df$mean_expr, corr_df$mean_ht, use = "complete.obs")
r2_all <- r_all^2

cat(paste0("\nAll lines: r = ", round(r_all, 3), ", r² = ", round(r2_all, 3), "\n"))

# -- scatter: all lines, colored by introgression status -----------------------

p_corr_status <- ggplot(corr_df, aes(x = mean_expr, y = mean_ht)) +
  geom_point(
    aes(fill = allele),
    shape = 21,
    size = 3,
    alpha = 0.75,
    color = "grey30"
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_fill_manual(
    values = c("B73" = "grey70", "Teosinte" = "darkorange"),
    name = "Allele at gene"
  ) +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("r² = ", round(r2_all, 3)),
    hjust = 1.2, vjust = 1.5, size = 5
  ) +
  labs(
    x = paste0(FOCUS_GENE, " expression (log2 CPM)"),
    y = paste0("plant height on ", HEIGHT_DATE, " (cm)"),
    title = paste("Height vs expression:", FOCUS_GENE),
    subtitle = paste0("n = ", nrow(corr_df), " genotypes; height on ", HEIGHT_DATE)
  ) +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold"))

print(p_corr_status)

ggsave(
  file.path(corr_dir, paste0(FOCUS_GENE, "_height_vs_expr_by_allele.png")),
  p_corr_status,
  width = 7,
  height = 6
)

# -- scatter: all lines, colored by taxa ---------------------------------------

p_corr_taxa <- ggplot(corr_df, aes(x = mean_expr, y = mean_ht)) +
  geom_point(
    aes(fill = taxa),
    shape = 21,
    size = 3,
    alpha = 0.75,
    color = "grey30"
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_fill_manual(values = taxa_colors, name = "Taxa") +
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0("r² = ", round(r2_all, 3)),
    hjust = 1.2, vjust = 1.5, size = 5
  ) +
  labs(
    x = paste0(FOCUS_GENE, " expression (log2 CPM)"),
    y = paste0("plant height on ", HEIGHT_DATE, " (cm)"),
    title = paste("Height vs expression:", FOCUS_GENE),
    subtitle = paste0("n = ", nrow(corr_df), " genotypes; colored by taxa")
  ) +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold"))

print(p_corr_taxa)

ggsave(
  file.path(corr_dir, paste0(FOCUS_GENE, "_height_vs_expr_by_taxa.png")),
  p_corr_taxa,
  width = 7,
  height = 6
)

# -- scatter: teo carriers only, colored by taxa -------------------------------

corr_teo <- corr_df |> filter(has_teo)

if (nrow(corr_teo) >= 3) {
  r_teo <- cor(corr_teo$mean_expr, corr_teo$mean_ht, use = "complete.obs")
  r2_teo <- r_teo^2

  p_corr_teo <- ggplot(corr_teo, aes(x = mean_expr, y = mean_ht)) +
    geom_point(
      aes(fill = taxa),
      shape = 21,
      size = 3.5,
      alpha = 0.8,
      color = "grey30"
    ) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    scale_fill_manual(values = taxa_colors, name = "Taxa") +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("r² = ", round(r2_teo, 3)),
      hjust = 1.2, vjust = 1.5, size = 5
    ) +
    labs(
      x = paste0(FOCUS_GENE, " expression (log2 CPM)"),
      y = paste0("plant height on ", HEIGHT_DATE, " (cm)"),
      title = paste("Teo carriers only — height vs expression:", FOCUS_GENE),
      subtitle = paste0("n = ", nrow(corr_teo), " genotypes carrying teosinte at this locus")
    ) +
    theme_minimal(base_size = 16) +
    theme(plot.title = element_text(face = "bold"))

  print(p_corr_teo)

  ggsave(
    file.path(corr_dir, paste0(FOCUS_GENE, "_height_vs_expr_teo_only.png")),
    p_corr_teo,
    width = 7,
    height = 6
  )
}

# -- scatter: B73-background only, colored by taxa -----------------------------

corr_b73 <- corr_df |> filter(!has_teo)

if (nrow(corr_b73) >= 3) {
  r_b73 <- cor(corr_b73$mean_expr, corr_b73$mean_ht, use = "complete.obs")
  r2_b73 <- r_b73^2

  p_corr_b73 <- ggplot(corr_b73, aes(x = mean_expr, y = mean_ht)) +
    geom_point(
      aes(fill = taxa),
      shape = 21,
      size = 3.5,
      alpha = 0.8,
      color = "grey30"
    ) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    scale_fill_manual(values = taxa_colors, name = "Taxa") +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("r² = ", round(r2_b73, 3)),
      hjust = 1.2, vjust = 1.5, size = 5
    ) +
    labs(
      x = paste0(FOCUS_GENE, " expression (log2 CPM)"),
      y = paste0("plant height on ", HEIGHT_DATE, " (cm)"),
      title = paste("B73-background only — height vs expression:", FOCUS_GENE),
      subtitle = paste0("n = ", nrow(corr_b73), " genotypes without teosinte at this locus")
    ) +
    theme_minimal(base_size = 16) +
    theme(plot.title = element_text(face = "bold"))

  print(p_corr_b73)

  ggsave(
    file.path(corr_dir, paste0(FOCUS_GENE, "_height_vs_expr_b73_only.png")),
    p_corr_b73,
    width = 7,
    height = 6
  )
}
