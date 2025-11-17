library(ggplot2)
library(ggthemes)
library(dplyr)
library(agricolae)
library(multcompView)
library(ggridges)

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/R/phenotyping_figures/data")
BZea <- read.csv("D4_SPAD_edited.csv")

#set order and color
species_order <- c("B73", "Mex", "Bals", "Zdip", "Hueh", "Zlux")
species_colors <- c("#5f9bfe", "#02bec5", "#f463e2", "#f8756d", "#b69c00", "#00b836")

#set order and color
species_3 <- c("B73", "Zdip", "Zlux")
species_3_colors <- c("white", "#e4211c", "#367db8")

# exclude checks
BZea_df <- BZea %>%
  filter(!Species %in% c("Check"))

BZea_df_3 <- BZea_df %>%
  filter(Species %in% species_3)
         
#order
BZea_df$Species <- factor(BZea_df$Species, levels = species_order)

#plot with all SPAD1s overlapping
SPAD1_plotty <- ggplot(BZea_df, aes(x = SPAD1, fill=Species)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = species_colors) +
  scale_x_continuous(limits = c(0, 80)) +
  labs(x="Leaf Greenness Index", y="Density")

print(SPAD1_plotty)


#plot with B73, zdip, zlux SPAD1s overlapping
SPAD1_3 <- ggplot(BZea_df_3, aes(x = SPAD1, fill=Species)) +
  geom_density(alpha = 0.5, size=.5) +
  theme_minimal() +
  scale_fill_manual(values = species_3_colors) +
  scale_x_continuous(limits = c(0, 80)) +
  labs(x="Leaf Greenness Index", y="Density")

print(SPAD1_3)

#plot with B73, zdip, zlux SPAD2s overlapping
SPAD1_32 <- ggplot(BZea_df_3, aes(x = SPAD2, fill=Species)) +
  geom_density(alpha = 0.5, size=.5) +
  theme_minimal() +
  scale_fill_manual(values = species_3_colors) +
  scale_x_continuous(limits = c(0, 80)) +
  labs(x="Leaf Greenness Index", y="Density")

print(SPAD1_32)



#plot SPAD1 and 2 overlapped but species separate
overlap_plotty <- ggplot(BZea_df, aes(x = SPAD1, fill = "SPAD1")) +
  geom_density(alpha = 0.5) +
  geom_density(aes(x = SPAD2, fill = "SPAD2"), alpha = 0.5, linetype = 2) +
  facet_wrap(~Species, scales = "free_x", ncol = 6) +
  theme_minimal() +
  scale_fill_manual(values = c("#00b836", "#b69c00"),
                        name = "",
                        labels = c("SPAD1", "SPAD2")) +
  scale_x_continuous(limits = c(0, 80)) +
  labs(x="Leaf Greenness Index", y="Density")

print(overlap_plotty)


#----------2025 SPAD---------------------------------------------------

library(ggplot2)
library(ggthemes)
library(dplyr)
library(agricolae)
library(multcompView)
library(ggridges)


setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea phenotyping/BZeaPheno")
BZea_df <- read.csv("SPAD-CLY25-Fieldbook - B5_BZea_eval.csv")

BZea_df$SPAD <- as.numeric(BZea_df$SPAD)
BZea_df$species_mex <- factor(BZea_df$species_mex, levels = c("Mex", "Hueh", "Bals", "Zlux", "Zdip", "B73"))

species_order <- c("Zlux", "Mex", "Bals", "Zdip", "Hueh", "B73")
species_colors <- c("#00b836", "#5f9bfe", "#f463e2", "#f8756d", "#b69c00", "#02bec5")


# exclude checks
BZea_df <- BZea_df %>%
  filter(!species_mex %in% c("Check", "", NA))

#order
BZea_df$species_mex <- factor(BZea_df$species_mex, levels = species_order)

#plot with all SPAD1s overlapping
SPAD1_plotty <- ggplot(BZea_df, aes(x = SPAD, fill=species_mex)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = species_colors) +
  scale_x_continuous(limits = c(0, 80)) +
  labs(x="Leaf Greenness Index", y="Density")

print(SPAD1_plotty)


#------2025 ridge---------------------------------

#rainbow overlapping a little
ggplot(BZea_df, aes(x = SPAD, y = species_mex, fill = species_mex, height = stat(density))) + 
  geom_density_ridges(scale = 3, draw_baseline = FALSE, alpha=0.7)

#quantiles
ggplot(BZea_df, aes(x=SPAD, y=species_mex, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) +
  scale_fill_viridis_d(name = "Quartiles")

#points
q <- ggplot(BZea_df, aes(x = SPAD, y = species_mex, fill = species_mex)) +
  geom_density_ridges(
    aes(point_color = species_mex, point_fill = species_mex, point_shape = species_mex),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 21, 21, 21, 21, 21)) +
  labs(x = "Leaf greenness index", y = "Species") + 
  scale_x_continuous(limits = c(-10, 80)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none"
  )

p <- ggplot(BZea_df, aes(x = SPAD, y = species_mex, fill = species_mex)) +
  geom_density_ridges(
    aes(point_color = species_mex, point_fill = species_mex, point_shape = species_mex),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 21, 21, 21, 21, 21)) +
  labs(x = "Leaf greenness index", y = "Species") +
  scale_x_continuous(limits = c(-10, 80)) +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "none"
  ) +
  geom_vline(
    xintercept = c(3.5, 16.4, 25.8, 37.7, 44.5, 54.8),
    linetype = "dashed",
    color = "gray60",
    alpha = 0.7
  )
p
ggsave(
  "C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea phenotyping/BZeaPheno/output/BZea_ridges_black.png",
  p,
  width = 9,
  height = 5.7,
  dpi = 300
)
#----------------2023 ridge-----------------------------------------------
BZea23 <- read.csv("UPDATED_CLY23_D4_FieldBook.csv")

#set order and color
species_order <- c("Zlux", "Mex", "Bals", "Zdip", "Hueh", "B73")
species_colors <- c("#00b836", "#5f9bfe", "#f463e2", "#f8756d", "#b69c00", "#02bec5")

# exclude checks
BZea_df23 <- BZea23 %>%
  filter(!species_mex %in% c("Check"))

#order
BZea_df23$species_mex <- factor(BZea_df23$species_mex, levels = species_order)

#points
ggplot(BZea_df23, aes(x = SPAD2, y = species_mex, fill = species_mex)) +
  geom_density_ridges(
    aes(point_color = species_mex, point_fill = species_mex, point_shape = species_mex),
    alpha = .2, point_alpha = 1, jittered_points = TRUE
  ) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_shape", values = c(21, 21, 21, 21, 21, 21)) +
  labs(x = "Leaf greenness index", y = "Species") +
  scale_x_continuous(limits = c(-10, 80))

