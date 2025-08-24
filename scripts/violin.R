# Load necessary libraries
library(ggplot2)
getwd()

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea phenotyping/BZeaPheno")

bzea25 <- read.csv("./CLY25_data_analysis - B5_BZea_eval.csv")
str(bzea25)

bzea25 %>% filter(is.na(species))
bzea25 %>% filter(species == "")
setdiff(unique(bzea25f$species), species_order)
unique(bzea25$species)


# filter out stuff from notes1
bzea25f <- bzea25 %>%
  filter((is.na(notes1) | notes1 == "") & !is.na(species) & species != "Check" & !is.na(DTS))

unique(bzea25f$species)
# Define the order of species
species_order <- c("Dura", "Nobo", "Mesa", "Chal", "B73", "Bals", "Zdip", "Hueh", "Zlux")
bzea25f$species <- factor(bzea25f$species, levels = species_order)

# convert columns to numeric
numeric_cols <- c("DTS", "DTA", "GDDTS", "GDDTA", "PH", "EH", "EN", "Prolif")
bzea25f[numeric_cols] <- lapply(bzea25f[numeric_cols], function(x) as.numeric(as.character(x)))

bzea25f$Rep <- as.factor(bzea25f$Rep)

species_colors <- c(
  "Dura" = "#609bfe",
  "Nobo" = "#609bfe",
  "Mesa" = "#609bfe",
  "Chal" = "#609bfe",
  "B73" = "#03bec4",
  "Bals" = "#f364e2",
  "Zdip" = "#f8756d",
  "Hueh" = "#b69d00",
  "Zlux" = "#00b837"
)

# Create the violin plot
plotty <- ggplot(bzea25f, aes(x = species, y = GDDTS, fill = species)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = species_colors) +  
  theme_minimal() +
  labs(title = "2025 GDDTS Distribution by Species",
       x = "Species",
       y = "GDD") +
  ylim(600,1200) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotty + stat_summary(fun.data=mean_sdl, mult=1, 
                      geom="pointrange", color="black")

#-------repeating with old data------------------------------------------------------


bzea23 <- read.csv("./UPDATED_CLY23_D4_FieldBook.csv")
str(bzea23)

bzea23 %>% filter(is.na(species))
bzea23 %>% filter(species == "")

# filter out stuff from notes1
bzea23f <- bzea23 %>%
  filter((is.na(notes1) | notes1 == "") & !is.na(species) & species != "Check" & !is.na(DTS))

unique(bzea23f$species)
# Define the order of species
species_order <- c("Dura", "Nobo", "Mesa", "Chal", "B73", "Bals", "Zdip", "Hueh", "Zlux")
bzea23f$species <- factor(bzea23f$species, levels = species_order)

# convert columns to numeric
numeric_cols <- c("GDDTS", "GDDTA", "DTS", "DTA", "PH", "EH", "EN", "Prolif")
bzea23f[numeric_cols] <- lapply(bzea23f[numeric_cols], function(x) as.numeric(as.character(x)))

bzea23f$Rep <- as.factor(bzea23f$Rep)

# Create the violin plot
plotty2 <- ggplot(bzea23f, aes(x = species, y = GDDTS, fill = species)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = species_colors) +  
  theme_minimal() +
  labs(title = "2023 GDDTS Distribution by Species",
       x = "Species",
       y = "GDD") +
  ylim(600,1200) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plotty2 + stat_summary(fun.data=mean_sdl, mult=1, 
                      geom="pointrange", color="black")


#-----------------combined----------------------------

# Combine EH and PH into long format for both years
bzea23_long <- bzea23f %>%
  select(species, GDDTS, GDDTA) %>%
  mutate(year = "2023")

bzea25_long <- bzea25f %>%
  select(species, GDDTS, GDDTA) %>%
  mutate(year = "2025")

# Combine both years
combined_long <- bind_rows(bzea23_long, bzea25_long)

# Pivot longer to combine EH and PH
library(tidyr)
combined_long <- pivot_longer(combined_long,
                              cols = c(GDDTS, GDDTA),
                              names_to = "trait",
                              values_to = "value")


# Re-set species factor levels
combined_long$species <- factor(combined_long$species,
                                levels = c("Dura", "Nobo", "Mesa", "Chal", "B73", "Bals", "Zdip", "Hueh", "Zlux"))

# Plot: PH and EH across years and species
ggplot(combined_long, aes(x = species, y = value, fill = species)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = mean_sdl, mult = 1,
               geom = "pointrange", color = "black") +
  scale_fill_manual(values = species_colors) +
  facet_grid(trait ~ year, scales = "free_y") +
  theme_minimal() +
  labs(title = "GDDTS and GDDTA Distributions by Species (2023 & 2025)",
       x = "Species", y = "GDD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(combined_long, aes(x = species, y = value, fill = species)) +
  geom_violin(trim = FALSE, 
              position = position_dodge(width = 0.8), 
              aes(group = interaction(species, year))) +
  stat_summary(fun.data = mean_sdl, mult = 1,
               geom = "pointrange", color = "black",
               position = position_dodge(width = 0.8),
               aes(group = interaction(species, year))) +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~trait, scales = "free_y") +
  theme_minimal() +
  guides(fill = "none") +
  labs(title = "GDDTS and GDDTA Distributions by Species and Year",
       x = "Species", y = "GDD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#------------DTS DTA---------------------
# Combine EH and PH into long format for both years
bzea23_long2 <- bzea23f %>%
  select(species, DTS, DTA) %>%
  mutate(year = "2023")

bzea25_long2 <- bzea25f %>%
  select(species, DTS, DTA) %>%
  mutate(year = "2025")

# Combine both years
combined_long2 <- bind_rows(bzea23_long2, bzea25_long2)

# Pivot longer to combine EH and PH
library(tidyr)
combined_long2 <- pivot_longer(combined_long2,
                              cols = c(DTS, DTA),
                              names_to = "trait",
                              values_to = "value")


# Re-set species factor levels
combined_long2$species <- factor(combined_long2$species,
                                levels = c("Dura", "Nobo", "Mesa", "Chal", "B73", "Bals", "Zdip", "Hueh", "Zlux"))

# Plot: PH and EH across years and species
ggplot(combined_long2, aes(x = species, y = value, fill = species)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.data = mean_sdl, mult = 1,
               geom = "pointrange", color = "black") +
  scale_fill_manual(values = species_colors) +
  facet_grid(trait ~ year, scales = "free_y") +
  theme_minimal() +
  labs(title = "DTS and DTA Distributions by Species (2023 & 2025)",
       x = "Species", y = "Days") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(combined_long2, aes(x = species, y = value, fill = species)) +
  geom_violin(trim = FALSE, 
              position = position_dodge(width = 0.8), 
              aes(group = interaction(species, year))) +
  stat_summary(fun.data = mean_sdl, mult = 1,
               geom = "pointrange", color = "black",
               position = position_dodge(width = 0.8),
               aes(group = interaction(species, year))) +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~trait, scales = "free_y") +
  theme_minimal() +
  guides(fill = "none") +
  labs(title = "GDDTS and GDDTA Distributions by Species and Year",
       x = "Species", y = "GDD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

