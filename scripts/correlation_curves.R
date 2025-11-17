library(dplyr)
library(ggplot2)

# read in dfs
D4 <- read.csv("corr_D4.csv",
               colClasses = c(
                 "CLY23_D4"   = "factor",
                 "Rep"        = "factor",
                 "Genotype"   = "factor",
                 "Species"    = "factor",
                 "origin"     = "factor",
                 "notes1"     = "factor",
                 "notes2"     = "factor",
                 "GDDTS_corr" = "numeric",
                 "GDDTA_corr" = "numeric",
                 "PH_corr"    = "numeric",
                 "EH_corr"    = "numeric",
                 "NBR_corr"   = "numeric",
                 "SPAD1_corr" = "numeric",
                 "SPAD2_corr" = "numeric"
               )
)

B5 <- read.csv("corr_B5.csv",
               colClasses = c(
                 "CLY25_B5"   = "factor",
                 "Rep"        = "factor",
                 "Genotype"   = "factor",
                 "accession"  = "factor",
                 "species"    = "factor",
                 "species_mex"= "factor",
                 "origin"     = "factor",
                 "notes1"     = "factor",
                 "notes2"     = "factor",
                 "GDDTS_corr" = "numeric",
                 "GDDTA_corr" = "numeric",
                 "PH_corr"    = "numeric",
                 "EH_corr"    = "numeric",
                 "EN_corr"    = "numeric",
                 "NBR_corr"   = "numeric",
                 "SPAD_corr"  = "numeric"
               )
)

# function to do correlation and plot
make_corr_plot <- function(D4, B5, trait, ylab, 
                           xlim = c(750, 1250), ylim = c(700, 950)) {
  
  col_D4 <- paste0(trait, "_2023")
  col_B5 <- paste0(trait, "_2025")
  
  # average raw
  B5_avg <- B5 |>
    group_by(Genotype) |>
    summarise("{col_B5}" := mean(.data[[paste0(trait, "_corr")]], na.rm = TRUE),
              .groups = "drop")
  
  D4_avg <- D4 |>
    group_by(Genotype) |>
    summarise("{col_D4}" := mean(.data[[paste0(trait, "_corr")]], na.rm = TRUE),
              .groups = "drop")
  
  TRAIT <- inner_join(D4_avg, B5_avg, by = "Genotype")
  
  # average without sickly
  B5_nosick <- B5 |>
    filter(is.na(notes1) | notes1 == "") |>
    group_by(Genotype) |>
    summarise("{col_B5}" := mean(.data[[paste0(trait, "_corr")]], na.rm = TRUE),
              .groups = "drop")
  
  D4_nosick <- D4 |>
    filter(is.na(notes1) | notes1 == "") |>
    group_by(Genotype) |>
    summarise("{col_D4}" := mean(.data[[paste0(trait, "_corr")]], na.rm = TRUE),
              .groups = "drop")
  
  TRAIT_nosick <- inner_join(D4_nosick, B5_nosick, by = "Genotype")
  
  # average within 3 sds
  mean23 <- mean(TRAIT[[col_D4]], na.rm = TRUE)
  sd23   <- sd(TRAIT[[col_D4]], na.rm = TRUE)
  mean25 <- mean(TRAIT[[col_B5]], na.rm = TRUE)
  sd25   <- sd(TRAIT[[col_B5]], na.rm = TRUE)
  
  TRAIT_3SDs <- TRAIT |>
    filter(abs(.data[[col_D4]] - mean23) <= 3*sd23,
           abs(.data[[col_B5]] - mean25) <= 3*sd25)
  
  # tag
  TRAIT$Filter <- "Raw"
  TRAIT_nosick$Filter <- "Without sickly plots"
  TRAIT_3SDs$Filter <- "Within 3 SDs"
  
  TRAIT_all <- bind_rows(TRAIT, TRAIT_nosick, TRAIT_3SDs)
  
  # get correlations
  corrs <- TRAIT_all %>%
    group_by(Filter) %>%
    summarise(r = cor(.data[[col_D4]], .data[[col_B5]], use = "complete.obs"))
  
  # plot
  ggplot(TRAIT_all, aes_string(x = col_D4, y = col_B5)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    facet_wrap(~Filter) +
    labs(x = paste0("2023 ", ylab),
         y = paste0("2025 ", ylab)) +
    geom_text(data = corrs,
              aes(x = -Inf, y = Inf,
                  label = paste0("r = ", round(r, 4))),
              inherit.aes = FALSE,
              hjust = -0.1, vjust = 1.2, size = 6) +
    coord_cartesian(xlim = xlim, ylim = ylim)
}

# run for each trait
plot_GDDTS <- make_corr_plot(D4, B5, trait = "GDDTS", ylab = "GDDTS",
                             xlim = c(750, 1250), ylim = c(700, 950))

plot_GDDTA <- make_corr_plot(D4, B5, trait = "GDDTA", ylab = "GDDTA",
                             xlim = c(750, 1250), ylim = c(700, 950))

plot_PH <- make_corr_plot(D4, B5, trait = "PH", ylab = "PH (cm)",
                          xlim = c(80, 270), ylim = c(80, 270))

# example: show gddts plot
plot_GDDTS
plot_GDDTA
plot_PH
