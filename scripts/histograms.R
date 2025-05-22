install.packages("tidyverse")
library(tidyverse)
library(dplyr)

bzeadf <- read.csv("./UPDATED_CLY23_D4_FieldBook.csv")
  
bzeadf_f <- bzeadf %>% filter(Species != "Check" & Notes == "")
bzeadf_ff <- bzeadf_f %>%
  mutate(LA = BW * BL * 0.75)


ggplot(bzeadf_ff, aes(x = NBR, fill = Species)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.7, position = "identity") +
  facet_wrap(~ Species, scales = "free_y") +
  labs(
    title = "NBR Distribution by Species",
    x = "Number of nodes with brace roots (NBR)",
    y = "Count"
  ) +
  theme_minimal()
