# Load necessary libraries
library(ggplot2)

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/R/phenotyping_figures/data")
avg_values <- read.csv("CLY23-D4-clean_phenotypes_for_analysis - avg_values.csv")

# Define the order of species
species_order <- c("Dura", "Nabo", "Mesa", "Chal", "B73", "Bals", "Zdip", "Hueh", "Zlux")

# Convert Species to a factor and set the levels to ensure correct order
avg_values$Species <- factor(avg_values$Species, levels = species_order)

species_colors <- c(
  "Dura" = "#609bfe",
  "Nabo" = "#609bfe",
  "Mesa" = "#609bfe",
  "Chal" = "#609bfe",
  "B73" = "#03bec4",
  "Bals" = "#f364e2",
  "Zdip" = "#f8756d",
  "Hueh" = "#b69d00",
  "Zlux" = "#00b837"
)

avg_values_Zlux <- avg_values %>%
  filter(Species == "Zlux")

avg_values_nocheck <- avg_values %>%
  filter(Species != "B73")

str(avg_values)

avg_values$PH <- as.numeric(avg_values$PH)
avg_values$EH <- as.numeric(avg_values$EH)
avg_values$EN <- as.numeric(avg_values$EN)
avg_values$BW <- as.numeric(avg_values$BW)
avg_values$BL <- as.numeric(avg_values$BL)
avg_values$SL <- as.numeric(avg_values$SL)
avg_values$leaf_area <- as.numeric(avg_values$leaf_area)

# Create the violin plot
plotty <- ggplot(avg_values, aes(x = Species, y = leaf_area, fill = Species)) +
  geom_violin(trim = FALSE) + 
  theme_minimal() + 
  labs(title = "Leaf Area Distribution by Species",
       x = "Species",
       y = "Leaf Area (cm²)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotty + stat_summary(fun.data=mean_sdl, mult=1, 
                      geom="pointrange", color="black")