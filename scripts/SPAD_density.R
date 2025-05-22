library(ggplot2)
library(ggthemes)
library(dplyr)
library(agricolae)
library(multcompView)


setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/R/phenotyping_figures/data")
BZea <- read.csv("D4_SPAD_edited.csv")

#set order and color
species_order <- c("Mex", "B73", "Bals", "Zdip", "Hueh", "Zlux")
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
