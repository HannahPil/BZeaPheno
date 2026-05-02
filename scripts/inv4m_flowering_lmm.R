# ==============================================================================
# INV4M FLOWERING-TIME LINEAR MIXED MODEL
#
# Tests whether the inv4m introgression significantly affects flowering time,
# accounting for the nested structure of the experiment:
#   plots (reps) within genotypes within families.
#
# Runs FOUR models per trait:
#   - BC1-only        : trait ~ inv4m + (1|bc1) + (1|bc1:Genotype)
#                       family = first 13 chars of Genotype
#   - BC2-only        : trait ~ inv4m + (1|bc2) + (1|bc2:Genotype)
#                       family = first 16 chars of Genotype
#   - BC1+BC2 nested  : trait ~ inv4m + (1|bc1) + (1|bc2) + (1|bc2:Genotype)
#   - Full nested     : trait ~ inv4m + (1|accession) + (1|f1) + (1|bc1) +
#                                       (1|bc2) + (1|bc2:Genotype)
#                       accession=7, F1=10, BC1=13, BC2=16 chars.
#
# All non-genotype levels are crossed random intercepts; each is a strict
# prefix of the next so the levels are perfectly nested by construction.
# (The full-nested model may produce singular fits if a level has zero
# variance after the others are accounted for -- that's reported but does
# not invalidate the inv4m fixed-effect estimate.)
#
# Sanity check (each model):
#   - Permutation test: permute inv4m labels WITHIN family 1000x, refit,
#     record coefficient.  Multilevel models permute within the finest
#     pedigree level (BC2).
# ==============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)

set.seed(42)

# ------------------------------- output dir -----------------------------------

out_dir <- file.path("output", "inv4m_flowering_lmm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- load + prep ----------------------------------

corr <- read.csv("data/corr_B5.csv", stringsAsFactors = FALSE)

# drop known outlier plots (same list as the effect-size script)
OUTLIER_PLOTS <- c(2756, 453, 304)
dropped <- corr |> filter(Plot %in% OUTLIER_PLOTS)
cat(sprintf("Dropping %d outlier plots: %s\n",
            nrow(dropped), paste(OUTLIER_PLOTS, collapse = ", ")))
corr <- corr |> filter(!Plot %in% OUTLIER_PLOTS)

# ==============================================================================
# Per-level pipeline: fit LMMs + permutation tests for one family grouping
# ==============================================================================

run_pipeline <- function(family_tag, n_chars, family_label) {

  cat(sprintf("\n################ %s (first %d chars) ################\n",
              family_label, n_chars))

  df <- corr |>
    filter(!Genotype %in% c("B73", "Purple Check"), !is.na(Genotype)) |>
    filter(inv4m_introgression %in% c("yes", "no")) |>
    mutate(
      family = substr(Genotype, 1, n_chars),
      inv4m  = factor(inv4m_introgression, levels = c("no", "yes"))
    )

  both_arms <- df |>
    distinct(family, inv4m) |>
    count(family) |>
    filter(n == 2) |>
    pull(family)

  df <- df |> filter(family %in% both_arms)

  cat(sprintf("Usable %s families (both arms present): %d\n",
              family_label, length(both_arms)))
  cat(sprintf("Total plots: %d\n", nrow(df)))

  # ---- LMM fit ----
  fit_lmm <- function(trait) {
    f <- as.formula(paste0(
      trait, " ~ inv4m + (1 | family) + (1 | family:Genotype)"
    ))
    m <- lmerTest::lmer(f, data = df)

    fe <- summary(m)$coefficients
    row_inv4m <- fe["inv4myes", ]
    ci <- confint(m, parm = "inv4myes", method = "Wald")
    vc <- as.data.frame(VarCorr(m))

    list(
      model    = m,
      trait    = trait,
      estimate = unname(row_inv4m["Estimate"]),
      se       = unname(row_inv4m["Std. Error"]),
      df       = unname(row_inv4m["df"]),
      t        = unname(row_inv4m["t value"]),
      p        = unname(row_inv4m["Pr(>|t|)"]),
      ci_lo    = ci[1, 1],
      ci_hi    = ci[1, 2],
      vc       = vc
    )
  }

  cat(sprintf("\n===== Fitting LMMs (%s) =====\n", family_label))
  fits <- lapply(c("DTS_corr", "DTA_corr"), fit_lmm)
  names(fits) <- c("DTS_corr", "DTA_corr")

  # ---- permutation test ----
  permute_within_family <- function(d) {
    d |>
      group_by(family) |>
      mutate(inv4m = sample(inv4m)) |>
      ungroup()
  }

  perm_test <- function(trait, n_perm = 1000) {
    observed <- fits[[trait]]$estimate
    f <- as.formula(paste0(
      trait, " ~ inv4m + (1 | family) + (1 | family:Genotype)"
    ))
    null_coefs <- numeric(n_perm)
    pb_every <- max(1, n_perm %/% 10)
    for (i in seq_len(n_perm)) {
      d_perm <- permute_within_family(df)
      m_perm <- suppressMessages(suppressWarnings(lme4::lmer(f, data = d_perm)))
      null_coefs[i] <- fixef(m_perm)["inv4myes"]
      if (i %% pb_every == 0) cat(sprintf("  %s %s perm %d/%d\n",
                                          family_tag, trait, i, n_perm))
    }
    p_perm <- (sum(abs(null_coefs) >= abs(observed)) + 1) / (n_perm + 1)
    list(observed = observed, null = null_coefs, p_perm = p_perm)
  }

  cat(sprintf("\n===== Permutation tests (%s, 1000 perms each) =====\n", family_label))
  perms <- lapply(c("DTS_corr", "DTA_corr"), perm_test, n_perm = 1000)
  names(perms) <- c("DTS_corr", "DTA_corr")

  # ---- report ----
  report_one <- function(trait) {
    f  <- fits[[trait]]
    pr <- perms[[trait]]
    cat(sprintf("\n---- %s | %s ----\n", family_label, trait))
    cat(sprintf("  LMM fixed-effect estimate (yes - no): %+.3f days\n", f$estimate))
    cat(sprintf("  95%% CI: [%+.3f, %+.3f]\n", f$ci_lo, f$ci_hi))
    cat(sprintf("  t(df=%.1f) = %.2f,  p = %.4g  (Satterthwaite)\n",
                f$df, f$t, f$p))
    cat(sprintf("  permutation p (1000 within-family perms) = %.4g\n",
                pr$p_perm))
    cat("  variance components:\n")
    for (i in seq_len(nrow(f$vc))) {
      cat(sprintf("    %-32s sd = %.3f  (var = %.3f)\n",
                  f$vc$grp[i], f$vc$sdcor[i], f$vc$vcov[i]))
    }
  }

  cat(sprintf("\n===== RESULTS (%s) =====\n", family_label))
  for (tr in c("DTS_corr", "DTA_corr")) report_one(tr)

  # ---- tidy results table ----
  results_tbl <- tibble(
    family_level = family_label,
    trait        = names(fits),
    estimate     = sapply(fits, `[[`, "estimate"),
    se           = sapply(fits, `[[`, "se"),
    ci_lo        = sapply(fits, `[[`, "ci_lo"),
    ci_hi        = sapply(fits, `[[`, "ci_hi"),
    df           = sapply(fits, `[[`, "df"),
    t            = sapply(fits, `[[`, "t"),
    p_lmm        = sapply(fits, `[[`, "p"),
    p_perm       = sapply(perms, `[[`, "p_perm")
  )

  # ---- diagnostic plots ----
  coef_df <- results_tbl |>
    mutate(trait_label = c(
      DTS_corr = "days to silking",
      DTA_corr = "days to anthesis"
    )[trait])

  p_coef <- ggplot(coef_df, aes(x = estimate, y = trait_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.8) +
    geom_point(size = 4, color = "#551a8b") +
    geom_text(aes(label = sprintf("p=%.3g", p_lmm)),
              nudge_y = 0.25, size = 4) +
    labs(
      x = "Effect of inv4m on flowering time (days, yes \u2212 no)",
      y = NULL,
      title    = paste0("inv4m effect on flowering time (LMM, ", family_label, ")"),
      subtitle = "Fixed effect + 95% Wald CI; family & genotype as random effects"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p_coef)
  ggsave(file.path(out_dir, paste0(family_tag, "_lmm_coefficient_plot.png")),
         p_coef, width = 7, height = 3.5)

  null_df <- bind_rows(lapply(names(perms), function(tr) {
    tibble(trait = tr, null = perms[[tr]]$null, observed = perms[[tr]]$observed)
  }))

  p_null <- ggplot(null_df, aes(x = null)) +
    geom_histogram(bins = 40, fill = "grey70", color = "white") +
    geom_vline(aes(xintercept = observed), color = "#551a8b", linewidth = 1) +
    facet_wrap(~ trait, scales = "free") +
    labs(
      x = "coefficient under within-family permutation",
      y = "count",
      title    = paste0("Permutation null distribution (", family_label, ")"),
      subtitle = "Purple line = observed coefficient from real data"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p_null)
  ggsave(file.path(out_dir, paste0(family_tag, "_permutation_null.png")),
         p_null, width = 9, height = 4)

  results_tbl
}

# ==============================================================================
# Nested pipeline: BOTH BC1 and BC2 as random effects in a single model.
# Since BC2 is strictly nested within BC1 (a BC2 family belongs to one BC1),
# we add both levels plus genotype-in-BC2 as random intercepts.
# ==============================================================================

run_nested_pipeline <- function() {

  family_tag   <- "bc1_bc2_nested"
  family_label <- "BC1 + BC2 nested"

  cat(sprintf("\n################ %s ################\n", family_label))

  df <- corr |>
    filter(!Genotype %in% c("B73", "Purple Check"), !is.na(Genotype)) |>
    filter(inv4m_introgression %in% c("yes", "no")) |>
    mutate(
      bc1_family = substr(Genotype, 1, 13),
      bc2_family = substr(Genotype, 1, 16),
      inv4m      = factor(inv4m_introgression, levels = c("no", "yes"))
    )

  # keep BC1 families that have both arms (same filter as BC1-only run)
  both_arms_bc1 <- df |>
    distinct(bc1_family, inv4m) |>
    count(bc1_family) |>
    filter(n == 2) |>
    pull(bc1_family)
  df <- df |> filter(bc1_family %in% both_arms_bc1)

  cat(sprintf("Plots: %d  |  BC1 families: %d  |  BC2 families: %d\n",
              nrow(df),
              length(unique(df$bc1_family)),
              length(unique(df$bc2_family))))

  # ---- LMM fit ----
  fit_lmm <- function(trait) {
    f <- as.formula(paste0(
      trait,
      " ~ inv4m + (1 | bc1_family) + (1 | bc2_family)",
      " + (1 | bc2_family:Genotype)"
    ))
    m <- lmerTest::lmer(f, data = df)

    fe <- summary(m)$coefficients
    row_inv4m <- fe["inv4myes", ]
    ci <- confint(m, parm = "inv4myes", method = "Wald")
    vc <- as.data.frame(VarCorr(m))

    list(
      model    = m,
      trait    = trait,
      estimate = unname(row_inv4m["Estimate"]),
      se       = unname(row_inv4m["Std. Error"]),
      df       = unname(row_inv4m["df"]),
      t        = unname(row_inv4m["t value"]),
      p        = unname(row_inv4m["Pr(>|t|)"]),
      ci_lo    = ci[1, 1],
      ci_hi    = ci[1, 2],
      vc       = vc
    )
  }

  cat(sprintf("\n===== Fitting LMMs (%s) =====\n", family_label))
  fits <- lapply(c("DTS_corr", "DTA_corr"), fit_lmm)
  names(fits) <- c("DTS_corr", "DTA_corr")

  # ---- permutation test: permute inv4m WITHIN BC2 family (finest level) ----
  permute_within_bc2 <- function(d) {
    d |>
      group_by(bc2_family) |>
      mutate(inv4m = sample(inv4m)) |>
      ungroup()
  }

  perm_test <- function(trait, n_perm = 1000) {
    observed <- fits[[trait]]$estimate
    f <- as.formula(paste0(
      trait,
      " ~ inv4m + (1 | bc1_family) + (1 | bc2_family)",
      " + (1 | bc2_family:Genotype)"
    ))
    null_coefs <- numeric(n_perm)
    pb_every <- max(1, n_perm %/% 10)
    for (i in seq_len(n_perm)) {
      d_perm <- permute_within_bc2(df)
      m_perm <- suppressMessages(suppressWarnings(lme4::lmer(f, data = d_perm)))
      null_coefs[i] <- fixef(m_perm)["inv4myes"]
      if (i %% pb_every == 0) cat(sprintf("  %s %s perm %d/%d\n",
                                          family_tag, trait, i, n_perm))
    }
    p_perm <- (sum(abs(null_coefs) >= abs(observed)) + 1) / (n_perm + 1)
    list(observed = observed, null = null_coefs, p_perm = p_perm)
  }

  cat(sprintf("\n===== Permutation tests (%s, 1000 perms each) =====\n",
              family_label))
  perms <- lapply(c("DTS_corr", "DTA_corr"), perm_test, n_perm = 1000)
  names(perms) <- c("DTS_corr", "DTA_corr")

  # ---- report ----
  report_one <- function(trait) {
    f  <- fits[[trait]]
    pr <- perms[[trait]]
    cat(sprintf("\n---- %s | %s ----\n", family_label, trait))
    cat(sprintf("  LMM fixed-effect estimate (yes - no): %+.3f days\n", f$estimate))
    cat(sprintf("  95%% CI: [%+.3f, %+.3f]\n", f$ci_lo, f$ci_hi))
    cat(sprintf("  t(df=%.1f) = %.2f,  p = %.4g  (Satterthwaite)\n",
                f$df, f$t, f$p))
    cat(sprintf("  permutation p (1000 within-BC2 perms) = %.4g\n",
                pr$p_perm))
    cat("  variance components:\n")
    for (i in seq_len(nrow(f$vc))) {
      cat(sprintf("    %-32s sd = %.3f  (var = %.3f)\n",
                  f$vc$grp[i], f$vc$sdcor[i], f$vc$vcov[i]))
    }
  }

  cat(sprintf("\n===== RESULTS (%s) =====\n", family_label))
  for (tr in c("DTS_corr", "DTA_corr")) report_one(tr)

  results_tbl <- tibble(
    family_level = family_label,
    trait        = names(fits),
    estimate     = sapply(fits, `[[`, "estimate"),
    se           = sapply(fits, `[[`, "se"),
    ci_lo        = sapply(fits, `[[`, "ci_lo"),
    ci_hi        = sapply(fits, `[[`, "ci_hi"),
    df           = sapply(fits, `[[`, "df"),
    t            = sapply(fits, `[[`, "t"),
    p_lmm        = sapply(fits, `[[`, "p"),
    p_perm       = sapply(perms, `[[`, "p_perm")
  )

  # ---- diagnostic plots ----
  coef_df <- results_tbl |>
    mutate(trait_label = c(
      DTS_corr = "days to silking",
      DTA_corr = "days to anthesis"
    )[trait])

  p_coef <- ggplot(coef_df, aes(x = estimate, y = trait_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.8) +
    geom_point(size = 4, color = "#551a8b") +
    geom_text(aes(label = sprintf("p=%.3g", p_lmm)),
              nudge_y = 0.25, size = 4) +
    labs(
      x = "Effect of inv4m on flowering time (days, yes \u2212 no)",
      y = NULL,
      title    = paste0("inv4m effect on flowering time (LMM, ", family_label, ")"),
      subtitle = "Fixed effect + 95% Wald CI; BC1, BC2, genotype all random"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p_coef)
  ggsave(file.path(out_dir, paste0(family_tag, "_lmm_coefficient_plot.png")),
         p_coef, width = 7, height = 3.5)

  null_df <- bind_rows(lapply(names(perms), function(tr) {
    tibble(trait = tr, null = perms[[tr]]$null, observed = perms[[tr]]$observed)
  }))

  p_null <- ggplot(null_df, aes(x = null)) +
    geom_histogram(bins = 40, fill = "grey70", color = "white") +
    geom_vline(aes(xintercept = observed), color = "#551a8b", linewidth = 1) +
    facet_wrap(~ trait, scales = "free") +
    labs(
      x = "coefficient under within-BC2 permutation",
      y = "count",
      title    = paste0("Permutation null distribution (", family_label, ")"),
      subtitle = "Purple line = observed coefficient from real data"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p_null)
  ggsave(file.path(out_dir, paste0(family_tag, "_permutation_null.png")),
         p_null, width = 9, height = 4)

  results_tbl
}

# ==============================================================================
# Fully-nested pipeline: accession | F1 | BC1 | BC2 + genotype
# All four pedigree levels (each a strict prefix of the next) plus genotype
# as crossed random intercepts.
#
# Heads-up: this is more random-effect levels than the data may support.
# Likely outcomes:
#   - one or more variance components collapse to 0 (singular fit warning)
#   - the inv4m fixed-effect estimate stays close to the simpler models, but
#     CI may widen and Satterthwaite df may drop sharply
# If singular fit fires, the model is still usable for the inv4m coefficient,
# but the affected variance component is uninformative.
# ==============================================================================

run_full_nested_pipeline <- function() {

  family_tag   <- "full_nested"
  family_label <- "accession | F1 | BC1 | BC2 nested"

  cat(sprintf("\n################ %s ################\n", family_label))

  df <- corr |>
    filter(!Genotype %in% c("B73", "Purple Check"), !is.na(Genotype)) |>
    filter(inv4m_introgression %in% c("yes", "no")) |>
    mutate(
      accession  = substr(Genotype, 1,  7),
      f1         = substr(Genotype, 1, 10),
      bc1_family = substr(Genotype, 1, 13),
      bc2_family = substr(Genotype, 1, 16),
      inv4m      = factor(inv4m_introgression, levels = c("no", "yes"))
    )

  # filter at BC1 level (most permissive) -- same as nested pipeline above
  both_arms_bc1 <- df |>
    distinct(bc1_family, inv4m) |>
    count(bc1_family) |>
    filter(n == 2) |>
    pull(bc1_family)
  df <- df |> filter(bc1_family %in% both_arms_bc1)

  cat(sprintf("Plots: %d  |  accessions: %d  |  F1s: %d  |  BC1s: %d  |  BC2s: %d\n",
              nrow(df),
              length(unique(df$accession)),
              length(unique(df$f1)),
              length(unique(df$bc1_family)),
              length(unique(df$bc2_family))))

  # ---- LMM fit ----
  fit_lmm <- function(trait) {
    f <- as.formula(paste0(
      trait,
      " ~ inv4m + (1 | accession) + (1 | f1) + (1 | bc1_family)",
      " + (1 | bc2_family) + (1 | bc2_family:Genotype)"
    ))
    # use suppressMessages to keep the output tidy; we still see the fit
    m <- suppressMessages(lmerTest::lmer(f, data = df))

    # check for singular fit and report it
    if (lme4::isSingular(m)) {
      cat(sprintf("  [note] %s: singular fit -- at least one variance is at boundary 0\n",
                  trait))
    }

    fe <- summary(m)$coefficients
    row_inv4m <- fe["inv4myes", ]
    ci <- confint(m, parm = "inv4myes", method = "Wald")
    vc <- as.data.frame(VarCorr(m))

    list(
      model    = m,
      trait    = trait,
      estimate = unname(row_inv4m["Estimate"]),
      se       = unname(row_inv4m["Std. Error"]),
      df       = unname(row_inv4m["df"]),
      t        = unname(row_inv4m["t value"]),
      p        = unname(row_inv4m["Pr(>|t|)"]),
      ci_lo    = ci[1, 1],
      ci_hi    = ci[1, 2],
      vc       = vc
    )
  }

  cat(sprintf("\n===== Fitting LMMs (%s) =====\n", family_label))
  fits <- lapply(c("DTS_corr", "DTA_corr"), fit_lmm)
  names(fits) <- c("DTS_corr", "DTA_corr")

  # ---- permutation: still permute within BC2 (finest pedigree level) ----
  permute_within_bc2 <- function(d) {
    d |>
      group_by(bc2_family) |>
      mutate(inv4m = sample(inv4m)) |>
      ungroup()
  }

  perm_test <- function(trait, n_perm = 1000) {
    observed <- fits[[trait]]$estimate
    f <- as.formula(paste0(
      trait,
      " ~ inv4m + (1 | accession) + (1 | f1) + (1 | bc1_family)",
      " + (1 | bc2_family) + (1 | bc2_family:Genotype)"
    ))
    null_coefs <- numeric(n_perm)
    pb_every <- max(1, n_perm %/% 10)
    for (i in seq_len(n_perm)) {
      d_perm <- permute_within_bc2(df)
      m_perm <- suppressMessages(suppressWarnings(lme4::lmer(f, data = d_perm)))
      null_coefs[i] <- fixef(m_perm)["inv4myes"]
      if (i %% pb_every == 0) cat(sprintf("  %s %s perm %d/%d\n",
                                          family_tag, trait, i, n_perm))
    }
    p_perm <- (sum(abs(null_coefs) >= abs(observed)) + 1) / (n_perm + 1)
    list(observed = observed, null = null_coefs, p_perm = p_perm)
  }

  cat(sprintf("\n===== Permutation tests (%s, 1000 perms each) =====\n",
              family_label))
  perms <- lapply(c("DTS_corr", "DTA_corr"), perm_test, n_perm = 1000)
  names(perms) <- c("DTS_corr", "DTA_corr")

  # ---- report ----
  report_one <- function(trait) {
    f  <- fits[[trait]]
    pr <- perms[[trait]]
    cat(sprintf("\n---- %s | %s ----\n", family_label, trait))
    cat(sprintf("  LMM fixed-effect estimate (yes - no): %+.3f days\n", f$estimate))
    cat(sprintf("  95%% CI: [%+.3f, %+.3f]\n", f$ci_lo, f$ci_hi))
    cat(sprintf("  t(df=%.1f) = %.2f,  p = %.4g  (Satterthwaite)\n",
                f$df, f$t, f$p))
    cat(sprintf("  permutation p (1000 within-BC2 perms) = %.4g\n",
                pr$p_perm))
    cat("  variance components:\n")
    for (i in seq_len(nrow(f$vc))) {
      cat(sprintf("    %-32s sd = %.3f  (var = %.3f)\n",
                  f$vc$grp[i], f$vc$sdcor[i], f$vc$vcov[i]))
    }
  }

  cat(sprintf("\n===== RESULTS (%s) =====\n", family_label))
  for (tr in c("DTS_corr", "DTA_corr")) report_one(tr)

  results_tbl <- tibble(
    family_level = family_label,
    trait        = names(fits),
    estimate     = sapply(fits, `[[`, "estimate"),
    se           = sapply(fits, `[[`, "se"),
    ci_lo        = sapply(fits, `[[`, "ci_lo"),
    ci_hi        = sapply(fits, `[[`, "ci_hi"),
    df           = sapply(fits, `[[`, "df"),
    t            = sapply(fits, `[[`, "t"),
    p_lmm        = sapply(fits, `[[`, "p"),
    p_perm       = sapply(perms, `[[`, "p_perm")
  )

  # ---- diagnostic plots ----
  coef_df <- results_tbl |>
    mutate(trait_label = c(
      DTS_corr = "days to silking",
      DTA_corr = "days to anthesis"
    )[trait])

  p_coef <- ggplot(coef_df, aes(x = estimate, y = trait_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.8) +
    geom_point(size = 4, color = "#551a8b") +
    geom_text(aes(label = sprintf("p=%.3g", p_lmm)),
              nudge_y = 0.25, size = 4) +
    labs(
      x = "Effect of inv4m on flowering time (days, yes − no)",
      y = NULL,
      title    = paste0("inv4m effect on flowering time (LMM, ", family_label, ")"),
      subtitle = "Fixed effect + 95% Wald CI; all pedigree levels random"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p_coef)
  ggsave(file.path(out_dir, paste0(family_tag, "_lmm_coefficient_plot.png")),
         p_coef, width = 7, height = 3.5)

  null_df <- bind_rows(lapply(names(perms), function(tr) {
    tibble(trait = tr, null = perms[[tr]]$null, observed = perms[[tr]]$observed)
  }))

  p_null <- ggplot(null_df, aes(x = null)) +
    geom_histogram(bins = 40, fill = "grey70", color = "white") +
    geom_vline(aes(xintercept = observed), color = "#551a8b", linewidth = 1) +
    facet_wrap(~ trait, scales = "free") +
    labs(
      x = "coefficient under within-BC2 permutation",
      y = "count",
      title    = paste0("Permutation null distribution (", family_label, ")"),
      subtitle = "Purple line = observed coefficient from real data"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  print(p_null)
  ggsave(file.path(out_dir, paste0(family_tag, "_permutation_null.png")),
         p_null, width = 9, height = 4)

  results_tbl
}

# ==============================================================================
# Run all four models, combine results, write single CSV
# ==============================================================================

res_bc1         <- run_pipeline("bc1", 13, "BC1 family")
res_bc2         <- run_pipeline("bc2", 16, "BC2 family")
res_nested      <- run_nested_pipeline()
res_full_nested <- run_full_nested_pipeline()

all_results <- bind_rows(res_bc1, res_bc2, res_nested, res_full_nested)
write.csv(all_results,
          file.path(out_dir, "lmm_results.csv"),
          row.names = FALSE)

cat(sprintf("\nCombined results table: %s\n",
            file.path(out_dir, "lmm_results.csv")))
cat("Per-level plots written to:", out_dir, "\n")
