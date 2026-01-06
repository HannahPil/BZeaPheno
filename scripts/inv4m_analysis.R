library(lme4)
library(dplyr)
library(lmerTest)

inv4m_df <- read.csv("BZea_inv4m.csv")


inv4m_df <- inv4m_df |>
  transform(
    Plot_id = factor(Plot_id),
    Genotype = factor(Genotype),
    F1_ID = factor(F1_ID),
    Rep = factor(Rep),
    inv4m     = as.factor(inv4m),
    taxa     = as.factor(taxa),
    Accession   = as.factor(Accession),
    DTS_corr   = as.numeric(DTS_corr),
    DTA_corr   = as.numeric(DTA_corr),
    PH_corr = as.numeric(PH_corr)
    #GDDTS_BLUE   = as.numeric(GDDTS_BLUE),
    #GDDTA_BLUE   = as.numeric(GDDTA_BLUE)
  )


#filter to 3 SDs
inv4m_df <- inv4m_df |>
  filter(
    abs(DTS_corr - mean(DTS_corr, na.rm = TRUE)) <= 3*sd(DTS_corr, na.rm = TRUE),
    abs(DTA_corr - mean(DTA_corr, na.rm = TRUE)) <= 3*sd(DTA_corr, na.rm = TRUE)
  )

#filter to just mexicana
inv4m_df <- inv4m_df |>
  filter(
    !taxa %in% c("Zdip","Zlux","Check","B73","Bals", "Hueh")
  )

#set B73 as reference taxa
inv4m_df$taxa <- relevel(inv4m_df$taxa, ref = "Chal")

DTS_lmer <- lmer(
  DTS_corr ~ taxa + inv4m +  (1 | F1_ID),
  data = inv4m_df
)
summary(DTS_lmer)
summary(DTA_lmer)
summary(PH_lmer)


#----------filtered by taxa------------------
#nobo
nobo_only <- inv4m_df |> 
  filter(taxa == "Nobo")

nobo_DTA <- lmer(
  DTA_corr ~ inv4m + (1 | F1_ID),
  data = nobo_only
)
summary(nobo_DTA)


#mesa
mesa_only <- inv4m_df |> 
  filter(taxa == "Mesa")

mesa_DTA <- lmer(
  DTA_corr ~ inv4m + (1 | F1_ID),
  data = mesa_only
)
summary(mesa_DTA)


#chal
chal_only <- inv4m_df |> 
  filter(taxa == "Chal")

chal_DTA <- lmer(
  DTA_corr ~ inv4m + (1 | F1_ID),
  data = chal_only
)
summary(chal_DTA)


#dura
dura_only <- inv4m_df |> 
  filter(taxa == "Dura")

dura_DTA <- lmer(
  DTA_corr ~ inv4m + (1 | F1_ID),
  data = dura_only
)
summary(dura_DTA)

#----------------------------hypothesis test------------------------------------------
# =============================
# libraries
# =============================
library(tidyr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggbeeswarm)
library(ggtext)
library(ggpubr)

# =============================
# data setup
# =============================

# averaged reps *and* split by taxa
df <- inv4m_df |>
  filter(taxa %in% c("Nobo", "Dura", "Mesa", "Chal")) |>
  rename(
    donor = taxa,
    value = DTS_corr
  ) |>
  filter(!is.na(value)) |>
  group_by(donor, Genotype, inv4m) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    genotype = recode(inv4m, "no" = "CTRL", "yes" = "INV4M")
  )

# =============================
# hypothesis test
# =============================
htest <- df |>
  group_by(donor) |>
  filter(n_distinct(genotype) == 2) |>
  t_test(value ~ genotype) |>
  adjust_pvalue(method = "fdr") |>
  add_significance()

# add y positions for brackets
htest_plot <- htest |>
  add_y_position(scales = "free")

# =============================
# custom theme
# =============================
pheno_theme2 <- theme_classic(base_size = 22) +
  theme(
    plot.title = element_markdown(hjust = 0.5, face = "bold"),
    axis.title.y = element_markdown(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    strip.background = element_blank(),
    strip.text = element_markdown(size = 18, face = "bold"),
    legend.position = "none"
  )

# =============================
# relabel for plot
# =============================
df2 <- df |>
  mutate(
    genotype = forcats::fct_relevel(genotype, "CTRL", "INV4M")
  )

htest_plot2 <- htest_plot |>
  mutate(
    group1 = dplyr::recode(group1, "no" = "CTRL", "yes" = "INV4M"),
    group2 = dplyr::recode(group2, "no" = "CTRL", "yes" = "INV4M")
  )

# =============================
# plot
# =============================
ggplot(df2, aes(x = genotype, y = value, color = genotype)) +
  geom_boxplot(width = 0.25, linewidth = 1.5, alpha = 0) +
  geom_quasirandom(size = 2.5) +
  scale_color_manual(values = c("CTRL" = "gold", "INV4M" = "purple4")) +
  labs(
    y = "Corrected DTS",
    x = "Genotype"
  ) +
  stat_pvalue_manual(
    htest_plot2,
    tip.length = 0.01,
    step.height = 0.08,
    size = 8,
    bracket.size = 0.8
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  facet_wrap(~ donor, ncol = 4) +
  pheno_theme2

ggsave(
  filename = "output/B5_DTS_corr_inv4m_test_by_donor_avg.png",
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300
)


#-------hypothesis test combined--------------------------------

# ===============================
# libraries
# ===============================
library(tidyr)
library(rstatix)
library(ggplot2)
library(ggbeeswarm)
library(ggtext)
library(ggpubr)
library(dplyr)

# ===============================
# setup (filter + average reps)
# ===============================
df <- inv4m_df |>
  filter(taxa %in% c("Nobo", "Dura", "Mesa", "Chal")) |>
  rename(value = DTS_corr) |>
  filter(!is.na(value)) |>
  group_by(Genotype, inv4m) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    genotype = recode(inv4m, "no" = "CTRL", "yes" = "INV4M")
  )

# ===============================
# hypothesis test
# ===============================
htest <- df |>
  t_test(value ~ genotype) |>
  adjust_pvalue(method = "fdr") |>
  add_significance() |>
  add_y_position()

print(htest)

# ===============================
# theme
# ===============================
pheno_theme2 <- theme_classic(base_size = 22) +
  theme(
    plot.title = element_markdown(hjust = 0.5, face = "bold"),
    axis.title.y = element_markdown(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    strip.background = element_blank(),
    strip.text = element_markdown(size = 18, face = "bold"),
    legend.position = "none"
  )

# ===============================
# plot
# ===============================
ggplot(df, aes(
  x = genotype,
  y = value,
  color = genotype
)) +
  geom_boxplot(width = 0.25, linewidth = 1.5, alpha = 0) +
  geom_quasirandom(size = 2.5) +
  scale_color_manual(values = c("CTRL" = "gold", "INV4M" = "purple4")) +
  labs(
    y = "Corrected DTS",
    x = "Genotype"
  ) +
  stat_pvalue_manual(
    htest,
    tip.length = 0.01,
    step.increase = 0.08,
    size = 8,
    bracket.size = 0.8
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  pheno_theme2

#----------------------------hypothesis test (PH)------------------------------------------
# =============================
# libraries
# =============================
library(tidyr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggbeeswarm)
library(ggtext)
library(ggpubr)

# =============================
# data setup
# =============================

# averaged reps *and* split by taxa
df <- inv4m_df |>
  filter(taxa %in% c("Nobo", "Dura", "Mesa", "Chal")) |>
  rename(
    donor = taxa,
    value = PH_corr
  ) |>
  filter(!is.na(value)) |>
  group_by(donor, Genotype, inv4m) |>
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    genotype = recode(inv4m, "no" = "CTRL", "yes" = "INV4M")
  )

# =============================
# hypothesis test
# =============================
htest <- df |>
  group_by(donor) |>
  filter(n_distinct(genotype) == 2) |>
  t_test(value ~ genotype) |>
  adjust_pvalue(method = "fdr") |>
  add_significance()

# add y positions for brackets
htest_plot <- htest |>
  add_y_position(scales = "free")

# =============================
# custom theme
# =============================
pheno_theme2 <- theme_classic(base_size = 22) +
  theme(
    plot.title = element_markdown(hjust = 0.5, face = "bold"),
    axis.title.y = element_markdown(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    strip.background = element_blank(),
    strip.text = element_markdown(size = 18, face = "bold"),
    legend.position = "none"
  )

# =============================
# relabel for plot
# =============================
df2 <- df |>
  mutate(
    genotype = forcats::fct_relevel(genotype, "CTRL", "INV4M")
  )

htest_plot2 <- htest_plot |>
  mutate(
    group1 = dplyr::recode(group1, "no" = "CTRL", "yes" = "INV4M"),
    group2 = dplyr::recode(group2, "no" = "CTRL", "yes" = "INV4M")
  )

# =============================
# plot
# =============================
ggplot(df2, aes(x = genotype, y = value, color = genotype)) +
  geom_boxplot(width = 0.25, linewidth = 1.5, alpha = 0) +
  geom_quasirandom(size = 2.5) +
  scale_color_manual(values = c("CTRL" = "gold", "INV4M" = "purple4")) +
  labs(
    y = "Corrected PH",
    x = "Genotype"
  ) +
  stat_pvalue_manual(
    htest_plot2,
    tip.length = 0.01,
    step.height = 0.08,
    size = 8,
    bracket.size = 0.8
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  facet_wrap(~ donor, ncol = 4) +
  pheno_theme2

ggsave(
  filename = "output/B5_PH_corr_inv4m_test_by_donor.png",
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300
)



