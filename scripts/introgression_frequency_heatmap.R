library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

# format rds-derived object that contains teosinte genotype segments
teogeno <- readRDS("results_list_new_name.rds")

length(teogeno)

# -------------------- UPDATED METADATA SOURCE --------------------
# read in B5 corrected phenotype / metadata table
B5_corr <- read_csv("corr_B5.csv", show_col_types = FALSE)

# derive metadata with genotype names matching teogeno (add .B)
metadata <- B5_corr %>%
  distinct(Genotype, species) %>%
  transmute(
    genotype_raw = as.character(Genotype),
    genotype     = paste0(Genotype, ".B"),
    taxa         = as.character(species)
  ) %>%
  filter(!is.na(taxa))

# ---------------------------------------------------------------

# set window size (bp)
window_bp <- 4000000
cat("window size (bp):", window_bp, "\n")

# make long df of segments for all genotypes, attach taxa, omit NA taxa
seg_df <- imap_dfr(teogeno, ~{
  tibble(
    genotype = .y,
    chr      = .x$V1,
    start    = .x$V2,
    end      = .x$V3,
    region   = .x$V4
  )
}) %>%
  left_join(metadata, by = "genotype") %>%
  filter(!is.na(taxa))

# chromosome lengths from b73 spans
chr_lengths <- seg_df %>%
  filter(region == "B73") %>%
  group_by(chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  mutate(chr_num = as.integer(str_remove(chr, "chr"))) %>%
  arrange(chr_num) %>%
  select(chr, chr_len)

# build genome windows
windows <- chr_lengths %>%
  mutate(win_starts = map(chr_len, ~seq(1, .x, by = window_bp))) %>%
  unnest(win_starts) %>%
  rename(win_start = win_starts) %>%
  mutate(
    win_end = pmin(win_start + window_bp - 1, chr_len),
    win_mid = (win_start + win_end) / 2,
    win_id  = row_number()
  )

# introgression segments only
intro_seg <- seg_df %>%
  filter(region == "Introgression") %>%
  select(genotype, taxa, chr, start, end)

# map introgressions to windows (overlap join)
intro_windows <- intro_seg %>%
  inner_join(windows, by = "chr") %>%
  filter(start <= win_end, end >= win_start) %>%
  distinct(taxa, genotype, chr, win_id, win_mid)

# genotype counts per taxa (denominator)
taxa_n <- metadata %>%
  filter(genotype %in% names(teogeno)) %>%
  filter(!is.na(taxa)) %>%
  group_by(taxa) %>%
  summarise(n_genotypes = n_distinct(genotype), .groups = "drop")

# count introgressed genotypes per taxa x window
taxa_win_counts <- intro_windows %>%
  group_by(taxa, chr, win_id, win_mid) %>%
  summarise(n_intro_genotypes = n_distinct(genotype), .groups = "drop") %>%
  left_join(taxa_n, by = "taxa") %>%
  mutate(freq = n_intro_genotypes / n_genotypes)

# complete missing windows (freq = 0)
taxa_levels <- sort(unique(as.character(taxa_win_counts$taxa)))

taxa_win_full <- windows %>%
  select(chr, win_id, win_mid) %>%
  crossing(taxa = taxa_levels) %>%
  left_join(taxa_n, by = "taxa") %>%
  left_join(
    taxa_win_counts %>% select(taxa, chr, win_id, n_intro_genotypes, freq),
    by = c("taxa", "chr", "win_id")
  ) %>%
  mutate(
    n_intro_genotypes = replace_na(n_intro_genotypes, 0),
    freq = replace_na(freq, 0),
    chr  = factor(chr, levels = paste0("chr", 1:10)),
    taxa = factor(taxa, levels = taxa_levels)
  )

# function: plot one taxon (stacked chromosome bars)
plot_one_taxon <- function(taxon_to_plot) {
  df <- taxa_win_full %>%
    filter(taxa == taxon_to_plot) %>%
    mutate(
      x0 = win_mid - window_bp / 2,
      x1 = win_mid + window_bp / 2
    )
  
  ggplot(df) +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = 0, ymax = 1, fill = freq), color = NA) +
    facet_grid(rows = vars(chr), switch = "y") +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_x_continuous("position", labels = label_number(scale = 1e-6, suffix = " Mb")) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      limits = c(0, max(df$freq)),
      oob = squish
    ) +
    labs(
      title = paste0("Introgression frequency heatmap: ", taxon_to_plot),
      subtitle = paste0(
        "window size = ", format(window_bp, scientific = FALSE),
        " bp; fill = fraction of genotypes with introgression"
      ),
      fill = "freq"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      axis.ticks.y = element_blank()
    )
}

# output folder + save all taxa plots
out_dir <- "introgression_heatmaps_by_taxa"
dir.create(out_dir, showWarnings = FALSE)

plots_by_taxa <- vector("list", length(taxa_levels))
names(plots_by_taxa) <- taxa_levels

for (tx in taxa_levels) {
  p <- plot_one_taxon(tx)
  plots_by_taxa[[tx]] <- p
  
  ggsave(
    filename = file.path(out_dir, paste0("BZea_introgression_heatmap_alt_", tx, "_4Mb.png")),
    plot = p,
    width = 18,
    height = 13,
    dpi = 300
  )
}

# print one example
plots_by_taxa[[taxa_levels[2]]]

#-----------------------ALL TAXA COMBINED---------------------------------------------
# total genotypes included (after your taxa != NA filter)
n_total_genotypes <- seg_df %>%
  distinct(genotype) %>%
  nrow()
cat("genotypes included (non-NA taxa):", n_total_genotypes, "\n")

# introgression segments only (drop taxa)
intro_seg_all <- seg_df %>%
  filter(region == "Introgression") %>%
  select(genotype, chr, start, end)

# map introgressions to windows (overlap join) + distinct per genotype per window
intro_windows_all <- intro_seg_all %>%
  inner_join(windows, by = "chr") %>%
  filter(start <= win_end, end >= win_start) %>%
  distinct(genotype, chr, win_id)

# count introgressed genotypes per window (this is "number of lines")
win_counts_all <- intro_windows_all %>%
  group_by(chr, win_id) %>%
  summarise(n_intro_genotypes = n_distinct(genotype), .groups = "drop")

# complete missing windows (count = 0)
win_full_all <- windows %>%
  select(chr, win_id, win_start, win_end, win_mid) %>%
  left_join(win_counts_all, by = c("chr", "win_id")) %>%
  mutate(n_intro_genotypes = replace_na(n_intro_genotypes, 0))

# order chromosomes correctly
win_full_all <- win_full_all %>%
  mutate(chr = factor(chr, levels = paste0("chr", 1:10)))

fill_max <- max(win_full_all$n_intro_genotypes)
fill_min <- min(win_full_all$n_intro_genotypes)

p_all <- ggplot(win_full_all) +
  geom_rect(
    aes(xmin = win_start, xmax = win_end, ymin = 0, ymax = 1, fill = n_intro_genotypes),
    color = NA
  ) +
  facet_grid(rows = vars(chr), switch = "y") +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous("position", labels = label_number(scale = 1e-6, suffix = " Mb")) +
  scale_fill_viridis_c(
    option = "magma",
    limits = c(50, fill_max),
    oob = squish
  ) +
  labs(
    title = "Introgression count heatmap",
    subtitle = paste0(
      "window size = ", format(window_bp, scientific = FALSE),
      " bp; fill = number of genotypes with introgression overlapping window (n = ",
      n_total_genotypes, ")"
    ),
    fill = "n lines"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.ticks.y = element_blank()
  )

p_all

ggsave(
  filename = file.path("./output", paste0("BZea_introgression_heatmap_all_lines_4Mb.png")),
  plot = p_all,
  width = 18,
  height = 10,
  dpi = 300
)
