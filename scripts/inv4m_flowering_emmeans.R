# ==============================================================================
# INV4M FLOWERING-TIME EMMEANS ON inv4m_introgression
#
# Uses the nested model from inv4m_flowering_lmm.R
#   trait ~ inv4m + (1 | bc1_family) + (1 | bc2_family) + (1 | bc2_family:Genotype)
# and extracts emmeans + pairwise contrast for the two-level inv4m factor.
#
# Produces per trait (DTS_corr, DTA_corr):
#   - <trait>_inv4m_emmeans.csv           per-level marginal means with CI
#   - <trait>_inv4m_contrast.csv          yes-minus-no contrast with CI, t, p
#   - inv4m_emmeans_combined.png          combined visual of both traits
# ==============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

set.seed(42)

# ------------------------------- output dir -----------------------------------

out_dir <- file.path("output", "inv4m_flowering_emmeans")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- load + prep ----------------------------------

corr <- read.csv("data/corr_B5.csv", stringsAsFactors = FALSE)

OUTLIER_PLOTS <- c(2756, 453, 304)
dropped <- corr |> filter(Plot %in% OUTLIER_PLOTS)
cat(sprintf("Dropping %d outlier plots: %s\n",
            nrow(dropped), paste(OUTLIER_PLOTS, collapse = ", ")))
corr <- corr |> filter(!Plot %in% OUTLIER_PLOTS)

df <- corr |>
  filter(!Genotype %in% c("B73", "Purple Check"), !is.na(Genotype)) |>
  filter(inv4m_introgression %in% c("yes", "no")) |>
  mutate(
    bc1_family = substr(Genotype, 1, 13),
    bc2_family = substr(Genotype, 1, 16),
    inv4m      = factor(inv4m_introgression, levels = c("no", "yes"))
  )

both_arms_bc1 <- df |>
  distinct(bc1_family, inv4m) |>
  count(bc1_family) |>
  filter(n == 2) |>
  pull(bc1_family)
df <- df |> filter(bc1_family %in% both_arms_bc1)

cat(sprintf("Plots: %d  |  BC1: %d  |  BC2: %d  |  Genotypes: %d\n",
            nrow(df),
            length(unique(df$bc1_family)),
            length(unique(df$bc2_family)),
            length(unique(df$Genotype))))

# ==============================================================================
# Per-trait: fit nested model + emmeans + contrast
# ==============================================================================

run_trait <- function(trait) {

  cat(sprintf("\n################ %s ################\n", trait))

  f <- as.formula(paste0(
    trait,
    " ~ inv4m + (1 | bc1_family) + (1 | bc2_family)",
    " + (1 | bc2_family:Genotype)"
  ))
  m <- suppressMessages(lmerTest::lmer(f, data = df))
  cat("Random-effect SDs:\n")
  print(as.data.frame(VarCorr(m))[, c("grp", "sdcor")])

  emm <- emmeans(m, specs = ~ inv4m)
  emm_df <- as.data.frame(emm) |>
    as_tibble() |>
    rename(ci_lo = lower.CL, ci_hi = upper.CL)

  # Pairwise contrast: yes - no (and no - yes by convention, but we reverse to
  # keep the same sign convention as the other scripts: negative = inv4m earlier)
  ctr <- contrast(emm, method = "revpairwise")
  ctr_df <- as.data.frame(ctr) |>
    as_tibble() |>
    rename(df = df, t = t.ratio, p = p.value) |>
    mutate(
      ci_lo = estimate - qt(0.975, df) * SE,
      ci_hi = estimate + qt(0.975, df) * SE
    )

  write.csv(emm_df,
            file.path(out_dir, paste0(trait, "_inv4m_emmeans.csv")),
            row.names = FALSE)
  write.csv(ctr_df,
            file.path(out_dir, paste0(trait, "_inv4m_contrast.csv")),
            row.names = FALSE)

  cat(sprintf("\n---- %s summary ----\n", trait))
  cat("  EMMs:\n")
  print(emm_df)
  cat("  Contrast:\n")
  print(ctr_df)

  list(trait = trait, emm = emm_df, contrast = ctr_df)
}

results <- lapply(c("DTS_corr", "DTA_corr"), run_trait)
names(results) <- c("DTS_corr", "DTA_corr")

# ==============================================================================
# Combined plot across traits
# ==============================================================================

TRAIT_LABELS <- c(DTS_corr = "days to silking",
                  DTA_corr = "days to anthesis")

emm_all <- bind_rows(lapply(results, function(r) {
  r$emm |> mutate(trait = r$trait)
})) |>
  mutate(trait_label = TRAIT_LABELS[trait])

ctr_all <- bind_rows(lapply(results, function(r) {
  r$contrast |> mutate(trait = r$trait)
})) |>
  mutate(trait_label = TRAIT_LABELS[trait])

# panel 1: per-level EMMs with CI
p_emm <- ggplot(emm_all, aes(x = inv4m, y = emmean, colour = inv4m)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.1, linewidth = 0.8) +
  geom_point(size = 4) +
  facet_wrap(~ trait_label, scales = "free_y") +
  scale_colour_manual(values = c(no = "#ffd701", yes = "#551a8b"),
                      labels = c(no = "non-carrier sister", yes = "inv4m"),
                      guide  = "none") +
  labs(x = NULL, y = "Estimated marginal mean (days)",
       title    = "inv4m marginal means (nested LMM)",
       subtitle = "Error bars = 95% CI") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

# panel 2: contrasts with CI and p
p_ctr <- ggplot(ctr_all,
                aes(x = estimate, y = trait_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.8) +
  geom_point(size = 4, colour = "#551a8b") +
  geom_text(aes(label = sprintf("p = %.3g", p)),
            nudge_y = 0.22, size = 4) +
  labs(x = "Contrast (inv4m \u2212 non-carrier sister), days",
       y = NULL,
       title    = "inv4m contrast from emmeans",
       subtitle = "Negative = inv4m flowers earlier; 95% CI from Satterthwaite df") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_emm)
ggsave(file.path(out_dir, "inv4m_emmeans_panel.png"),
       p_emm, width = 8, height = 4)

print(p_ctr)
ggsave(file.path(out_dir, "inv4m_contrast_panel.png"),
       p_ctr, width = 7, height = 3.5)

cat(sprintf("\nAll outputs written to: %s\n", out_dir))
