# ==============================================================================
# INV4M FLOWERING-TIME EFFECT SIZE PLOTS (dabestr)
# Compare DTS and DTA between lines carrying the inv4m introgression and their
# non-carrier sisters within the same family.  Runs at TWO grouping levels:
#   - BC1 family: first 13 chars of Genotype
#   - BC2 family: first 16 chars of Genotype
#
# Produces, for each trait × each family level × plot- and genotype-level:
#   (a) group shared-control (all inv4m vs all non-carrier sisters)
#   (b) paired per-family (inv4m mean vs non-carrier sister mean, paired by family)
#   (c) per-family forest-style plot of effect sizes with CIs
#
# Effect sign convention: mean(non-carrier sister) - mean(inv4m)
#   positive value => inv4m carriers flower earlier than their non-carrier sisters.
# ==============================================================================


library(tidyverse)
library(dabestr)
library(patchwork)

set.seed(42)

# factor values & display labels
INV4M_LABEL  <- "inv4m"
SISTER_LABEL <- "non-carrier sister"

# palette: inv4m first (left), non-carrier sister second (right)
INV4M_COLOR  <- "#551a8b"
SISTER_COLOR <- "#ffd701"
PALETTE <- setNames(
  c(INV4M_COLOR, SISTER_COLOR),
  c(INV4M_LABEL, SISTER_LABEL)
)

# force ggplot's discrete defaults so any plot built downstream (including the
# ones dabestr returns) uses our palette unless it explicitly sets its own scale
options(
  ggplot2.discrete.colour = unname(PALETTE),
  ggplot2.discrete.fill   = unname(PALETTE)
)

# short human trait labels
TRAIT_LABELS <- c(
  DTS_corr = "days to silking (corrected)",
  DTA_corr = "days to anthesis (corrected)"
)

# ------------------------------- output dir -----------------------------------

out_dir <- file.path("output", "inv4m_flowering_effect_size")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- load data ------------------------------------

corr <- read.csv("data/corr_B5.csv", stringsAsFactors = FALSE)

# drop known outlier plots before any analysis
OUTLIER_PLOTS <- c(2756, 453, 304)
dropped <- corr |> filter(Plot %in% OUTLIER_PLOTS)
cat(sprintf("Dropping %d outlier plots: %s\n",
            nrow(dropped), paste(OUTLIER_PLOTS, collapse = ", ")))
corr <- corr |> filter(!Plot %in% OUTLIER_PLOTS)

# ------------------------------- helpers --------------------------------------

save_db_plot <- function(db_eff, file, title, subtitle,
                         width = 7, height = 5,
                         palette = PALETTE) {
  # dabestr bakes its internal subplots to grobs (layers of class
  # GeomDrawGrob) before returning, so adding a scale to the returned object
  # can't change the colours inside those grobs.  Instead, replace ggsci's
  # npg scale/palette functions with our manual equivalents for the duration
  # of this call, so dabestr picks up our colours when it builds the grobs.
  vals <- unname(palette)

  orig_scale_color <- ggsci::scale_color_npg
  orig_scale_fill  <- ggsci::scale_fill_npg
  orig_pal         <- ggsci::pal_npg
  on.exit({
    assignInNamespace("scale_color_npg", orig_scale_color, ns = "ggsci")
    assignInNamespace("scale_fill_npg",  orig_scale_fill,  ns = "ggsci")
    assignInNamespace("pal_npg",         orig_pal,         ns = "ggsci")
  }, add = TRUE)

  assignInNamespace(
    "scale_color_npg",
    function(...) ggplot2::scale_colour_manual(values = vals),
    ns = "ggsci"
  )
  assignInNamespace(
    "scale_fill_npg",
    function(...) ggplot2::scale_fill_manual(values = vals),
    ns = "ggsci"
  )
  assignInNamespace(
    "pal_npg",
    function(...) function(n) rep_len(vals, n),
    ns = "ggsci"
  )

  p <- dabestr::dabest_plot(db_eff)
  p <- p + patchwork::plot_annotation(
    title    = title,
    subtitle = subtitle,
    theme    = theme(plot.title = element_text(face = "bold"))
  )
  print(p)
  ggsave(file.path(out_dir, file), p, width = width, height = height)
  invisible(p)
}

drop_na_trait <- function(data, trait) {
  data[!is.na(data[[trait]]), , drop = FALSE]
}

safe_run <- function(expr, tag) {
  captured <- NULL
  tryCatch(
    withCallingHandlers(
      expr,
      error = function(e) {
        captured <<- sys.calls()
      },
      warning = function(w) {
        message(sprintf("[warn] %s: %s", tag, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      message(sprintf("[skip] %s: %s", tag, conditionMessage(e)))
      if (!is.null(captured)) {
        message("  full call stack at error:")
        for (i in seq_along(captured)) {
          line <- paste(deparse(captured[[i]]), collapse = " ")
          if (nchar(line) > 200) line <- paste0(substr(line, 1, 200), "...")
          message(sprintf("    [%2d] %s", i, line))
        }
      }
      NULL
    }
  )
}

bootstrap_diff <- function(sis_vals, teo_vals, n_boot = 5000) {
  # effect size: non-carrier sister - inv4m  (matches dabestr with idx = c("inv4m","non-carrier sister"))
  # NOTE: use sample.int-based indexing to avoid R's single-value sample() trap.
  # base::sample(c(70), replace = TRUE) returns an integer sampled from 1:70,
  # NOT the value 70. This silently blew up CIs for any arm with n = 1.
  boot_sample <- function(x) x[sample.int(length(x), size = length(x), replace = TRUE)]

  sis_vals <- sis_vals[!is.na(sis_vals)]
  teo_vals <- teo_vals[!is.na(teo_vals)]
  if (length(sis_vals) == 0 || length(teo_vals) == 0) {
    return(c(point = NA, lo = NA, hi = NA))
  }
  obs <- mean(sis_vals) - mean(teo_vals)
  boots <- replicate(n_boot, {
    mean(boot_sample(sis_vals)) - mean(boot_sample(teo_vals))
  })
  c(point = obs, lo = unname(quantile(boots, 0.025)),
    hi    = unname(quantile(boots, 0.975)))
}

# ==============================================================================
# Per-level pipeline: runs (a), (b), (c) for one family grouping (BC1 or BC2)
# ==============================================================================

run_pipeline <- function(family_tag, n_chars, family_label) {

  cat(sprintf("\n################ %s family level (first %d chars) ################\n",
              family_tag, n_chars))

  # ---- build data at this grouping ----
  df <- corr |>
    filter(!Genotype %in% c("B73", "Purple Check"), !is.na(Genotype)) |>
    mutate(
      family = substr(Genotype, 1, n_chars),
      allele = ifelse(inv4m_introgression == "yes", INV4M_LABEL, SISTER_LABEL)
    ) |>
    filter(allele %in% c(INV4M_LABEL, SISTER_LABEL))

  fam_summary <- df |>
    distinct(family, Genotype, allele) |>
    count(family, allele) |>
    pivot_wider(names_from = allele, values_from = n, values_fill = 0) |>
    filter(.data[[INV4M_LABEL]] >= 1, .data[[SISTER_LABEL]] >= 1)

  cat(sprintf("Usable %s families (>=1 inv4m and >=1 non-carrier sister line): %d\n",
              family_tag, nrow(fam_summary)))
  print(as.data.frame(fam_summary))

  df <- df |> filter(family %in% fam_summary$family)

  df_plot <- df
  df_geno <- df |>
    group_by(family, Genotype, allele) |>
    summarise(
      DTS_corr = mean(DTS_corr, na.rm = TRUE),
      DTA_corr = mean(DTA_corr, na.rm = TRUE),
      .groups  = "drop"
    )
  df_family <- df_geno |>
    group_by(family, allele) |>
    summarise(
      DTS_corr = mean(DTS_corr, na.rm = TRUE),
      DTA_corr = mean(DTA_corr, na.rm = TRUE),
      .groups  = "drop"
    ) |>
    group_by(family) |>
    filter(all(c(INV4M_LABEL, SISTER_LABEL) %in% allele)) |>
    ungroup()

  # ---- (a) group shared-control ----
  run_group <- function(data, trait, level_tag) {
    d <- drop_na_trait(data, trait) |> rename(value = !!sym(trait))
    db <- dabestr::load(
      data = d, x = allele, y = value,
      idx  = c(INV4M_LABEL, SISTER_LABEL)
    )
    eff <- dabestr::mean_diff(db)
    save_db_plot(
      eff,
      paste0(family_tag, "_", trait, "_", level_tag, "_a_group.png"),
      title    = paste0("Group comparison: ", TRAIT_LABELS[[trait]]),
      subtitle = paste0("All inv4m lines vs all non-carrier sisters (",
                        family_label, ", ", level_tag, "-level)"),
      width = 6, height = 5
    )
    eff
  }

  for (tr in c("DTS_corr", "DTA_corr")) {
    safe_run(run_group(df_plot, tr, "plot"),
             paste0("(a) ", family_tag, " group plot ", tr))
    safe_run(run_group(df_geno, tr, "geno"),
             paste0("(a) ", family_tag, " group geno ", tr))
  }

  # ---- (b) paired per-family ----
  run_paired_family <- function(trait) {
    d <- drop_na_trait(df_family, trait) |> rename(value = !!sym(trait))
    db <- dabestr::load(
      data     = d, x = allele, y = value,
      idx      = c(INV4M_LABEL, SISTER_LABEL),
      paired   = "baseline",
      id_col   = family
    )
    eff <- dabestr::mean_diff(db)
    save_db_plot(
      eff,
      paste0(family_tag, "_", trait, "_b_paired_family.png"),
      title    = paste0("Family-paired: ", TRAIT_LABELS[[trait]]),
      subtitle = paste0(family_label,
                        ": each family = one pair (inv4m mean vs non-carrier sister mean)"),
      width = 6, height = 5
    )
    eff
  }

  for (tr in c("DTS_corr", "DTA_corr")) {
    safe_run(run_paired_family(tr),
             paste0("(b) ", family_tag, " paired_family ", tr))
  }

  # ---- (c) forest plot ----
  forest_plot <- function(data, trait, level_tag) {
    fams <- unique(data$family)

    eff_df <- map_dfr(fams, function(f) {
      sub <- data |> filter(family == f)
      sis_vals <- sub[sub$allele == SISTER_LABEL, trait, drop = TRUE]
      teo_vals <- sub[sub$allele == INV4M_LABEL,  trait, drop = TRUE]
      r <- bootstrap_diff(sis_vals, teo_vals)
      tibble(
        family = f,
        n_sis  = length(sis_vals),
        n_teo  = length(teo_vals),
        point  = r["point"],
        lo     = r["lo"],
        hi     = r["hi"]
      )
    }) |>
      filter(!is.na(point)) |>
      arrange(point) |>
      mutate(
        family    = factor(family, levels = family),
        fam_label = paste0(as.character(family),
                           " (inv4m=", n_teo, ", nc=", n_sis, ")")
      ) |>
      mutate(fam_label = factor(fam_label, levels = fam_label))

    write.csv(
      eff_df,
      file.path(out_dir,
                paste0(family_tag, "_", trait, "_", level_tag, "_c_forest_effects.csv")),
      row.names = FALSE
    )

    p <- ggplot(eff_df, aes(x = point, y = fam_label)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7) +
      geom_point(size = 3, color = INV4M_COLOR) +
      labs(
        x     = paste0("Mean difference (non-carrier sister \u2212 inv4m), ",
                       TRAIT_LABELS[[trait]]),
        y     = paste0(family_label, " (n per arm)"),
        title = paste0("Per-family effect size: ", TRAIT_LABELS[[trait]],
                       " (", family_label, ", ", level_tag, "-level)"),
        subtitle = "Positive = inv4m flowers earlier; 95% bootstrap CI"
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))

    print(p)
    ggsave(
      file.path(out_dir,
                paste0(family_tag, "_", trait, "_", level_tag, "_c_forest.png")),
      p,
      width  = 7,
      height = max(4, nrow(eff_df) * 0.25 + 1.5)
    )

    eff_df
  }

  for (tr in c("DTS_corr", "DTA_corr")) {
    forest_plot(df_plot, tr, "plot")
    forest_plot(df_geno, tr, "geno")
  }
}

# ==============================================================================
# Run pipeline at both grouping levels
# ==============================================================================

run_pipeline("bc1", 13, "BC1 family")
run_pipeline("bc2", 16, "BC2 family")

cat("\nAll outputs written to:", out_dir, "\n")
