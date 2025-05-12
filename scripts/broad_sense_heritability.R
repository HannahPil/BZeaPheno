install.packages("lme4")
library(lme4)
library(dplyr)


bzeadf <- read.csv("./UPDATED_CLY23_D4_Fieldbook.csv")
head(bzeadf)

#------------------------test lme4 model with DTS-------------------------------

# Remove checks and missing values
bzea_clean <- bzeadf %>%
  filter(!(genotype_id %in% c("B73", "Purple Check"))) %>%  # exclude checks
  mutate(
    genotype_id = factor(genotype_id),
    Rep = factor(Rep),
    DTS = as.numeric(DTS)
  ) %>%
  filter(!is.na(DTS))  # remove missing DTS

# Fit the linear mixed model
dts_model <- lmer(DTS ~ (1|genotype_id) + (1|Rep), data = bzea_clean)

# Extract variance components
varcomp <- as.data.frame(VarCorr(dts_model))
vg <- varcomp$vcov[varcomp$grp == "genotype_id"]  # genetic variance
ve <- varcomp$vcov[varcomp$grp == "Residual"]     # residual variance

# Estimate H² with replication accounted for automatically by lmer
H2 <- vg / (vg + ve)
cat("Broad-sense heritability (H²) for DTS:", round(H2, 3), "\n")

#------------------------run for all traits-------------------------------------
