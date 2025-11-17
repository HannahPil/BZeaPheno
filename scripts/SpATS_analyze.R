install.packages("SpATS")
install.packages("statgenHTP")
library(SpATS)
library(statgenHTP)
library(dplyr)

#set up data

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea phenotyping/BZeaPheno")

manifest_B5 <- read.csv("B5_spats_manifest.csv")

manifest_B5 <- manifest_B5 |>
  transform(
    Plot_id = factor(Plot_id),
    Genotype = factor(Genotype),
    Rep = factor(Rep),
    Row     = as.numeric(Row),
    Range   = as.numeric(Range),
    DTS   = as.numeric(DTS),
    DTA   = as.numeric(DTA),
    GDDTS   = as.numeric(GDDTS),
    GDDTA   = as.numeric(GDDTA),
    PH   = as.numeric(PH),
    EN   = as.numeric(EN),
    EH   = as.numeric(EH),
    NBR   = as.numeric(NBR),
    SPAD   = as.numeric(SPAD)
  )

#----example of SpATS model run-----

#spats1 <- SpATS(
  #response = "GDDTS",
  #spatial  = ~ SAP(Row, Range, nseg = c(2,4), degree = 3, pord = 2),
  #genotype = "Genotype",
  #random   = ~ Rep,
  #data     = manifest_B5,
  #genotype.as.random = FALSE  # true for blups and false for blues
#)
#print(spats1)
#plot(spats1)

#-----statgen model run--------------------------------------------

manifest_B5$TP <- as.Date("2024-07-01")  # any arbitrary date

phenoTP_B5 <- createTimePoints(
  dat            = manifest_B5,
  experimentName = "B5_experiment",
  genotype       = "Genotype",
  timePoint      = "TP",
  repId          = "Rep",
  plotId         = "Plot_id",
  rowNum         = "Range",
  colNum         = "Row",
  addCheck       = TRUE,
  checkGenotypes = "B73",
  timeFormat     = "%Y-%m-%d"     # match fake date format
)


## Fit spatial model for one or more traits
fit_B5_GDDTS <- fitModels(
  TP       = phenoTP_B5,
  trait    = "SPAD",      # choose the trait you want to model
  engine   = "SpATS",
  useCheck = TRUE # splits checks as fixed, others random
)

## Plot spatial trend (optional visualization)
plot(fit_B5_GDDTS, timePoints = "2024-07-01", plotType = "spatial", spaTrend = "raw")

#-----------loop for all traits---------------
#BLUEs
traits_B5 <- c("DTS","DTA","GDDTS","GDDTA","PH","EN","EH","NBR","SPAD")

for(tr in traits_B5) {
  
  fit <- fitModels(
    TP = phenoTP_B5,
    trait = tr,
    engine = "SpATS",
    what = "fixed"
  )
  
  corr <- getCorrected(fit) |>
    rename(Plot_id = plotId) |>
    left_join(manifest_B5[,c("Plot_id","Genotype","Rep","Row","Range")], by="Plot_id") |>
    arrange(Plot_id)   # <-- sort by Plot_id
  
  write.csv(corr, paste0("B5_", tr, "_2025.csv"), row.names = FALSE)
}

#--BLUES--------------
for(tr in traits_B5) {
  
  fit <- fitModels(
    TP = phenoTP_B5,
    trait = tr,
    engine = "SpATS",
    what = "fixed"
  )
  
  geno <- getGenoPred(fit)$genoPred
  
  write.csv(geno, paste0("B5_", tr, "_BLUEs_2025.csv"), row.names = FALSE)
}

#--------------------predicted N analysis--------------------------------------
b5_N <- read.csv("B5_block2_predictedN.csv")

b5_N <- b5_N |>
  transform(
    Plot = factor(Plot),
    Genotype = factor(Genotype),
    Predicted_N = as.numeric(Predicted_N),
    Row     = as.numeric(Row),
    Range   = as.numeric(Range)
  )

b5_plot_avg <- b5_N |>
  group_by(Plot) |>
  summarise(
    Row           = first(Row),
    Range         = first(Range),
    Genotype      = first(Genotype),
    Accession     = first(Accession),
    Taxa          = first(Taxa),
    M_Distance    = mean(M_Distance),
    Predicted_N = mean(Predicted_N[M_Distance <= 4]),     # ONLY KEEP VALUES WHOSE M_DISTANCE ≤ 4
    M_Distance  = mean(M_Distance[M_Distance <= 4]),      #same for M-distance
    .groups = "drop"
  )

b5_plot_avg$TP <- as.Date("2024-07-01")  # any arbitrary date

pheno_b5_N <- createTimePoints(
  dat            = b5_plot_avg,
  experimentName = "B5_experiment",
  genotype       = "Genotype",
  timePoint      = "TP",
  plotId         = "Plot",
  rowNum         = "Range",
  colNum         = "Row",
  addCheck       = TRUE,
  checkGenotypes = "B73",
  timeFormat     = "%Y-%m-%d"     # match fake date format
)

fit_B5_N <- fitModels(
  TP       = pheno_b5_N,
  trait    = "Predicted_N",      
  engine   = "SpATS",
  useCheck = TRUE 
)

# Plot spatial trend (optional visualization)
plot(fit_B5_N, timePoints = "2024-07-01", plotType = "spatial", spaTrend = "raw")

corr_B5_N <- getCorrected(fit_B5_N) |>
  rename(Plot = plotId) |>
  left_join(
    b5_plot_avg |>
      dplyr::select(Plot, Genotype, Accession, Taxa, Row, Range),
    by = "Plot"
  )

write.csv(
  corr_B5_N,
  "B5_block2_PredictedN_corrected.csv",
  row.names = FALSE
)

#-----------------D4 ANALYSIS------------D4 ANALYSIS--------------------D4 ANALYSIS--------------------------D4 ANALYSIS------------------------------#

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea phenotyping/BZeaPheno")

manifest_D4 <- read.csv("D4_spats_manifest.csv")
#map <- read.csv("B5_spats_map.csv", header = FALSE)

manifest_D4 <- manifest_D4 |>
  transform(
    Plot_id = factor(Plot_id),
    Genotype = factor(Genotype),
    Rep = factor(Rep),
    Row     = as.numeric(Row),
    Range   = as.numeric(Range),
    GDDTS   = as.numeric(GDDTS),
    GDDTA   = as.numeric(GDDTA),
    PH   = as.numeric(PH),
    EN   = as.numeric(EN),
    EH   = as.numeric(EH),
    NBR   = as.numeric(NBR),
    SPAD1   = as.numeric(SPAD1),
    SPAD2   = as.numeric(SPAD2)
  )
# SpATS model run

spatsD4 <- SpATS(
  response = "PH",
  spatial  = ~ SAP(Row, Range, nseg = c(2,4), degree = 3, pord = 2),
  genotype = "Genotype",
  random   = ~ Rep,
  data     = manifest_D4,
  genotype.as.random = FALSE  # true for blups and false for blues
)
print(spatsD4)


##-----statgen model run--------------------------------------------

manifest_D4$TP <- as.Date("2024-07-01")  # any arbitrary date

phenoTP_D4 <- createTimePoints(
  dat            = manifest_D4,
  experimentName = "D4_experiment",
  genotype       = "Genotype",
  timePoint      = "TP",
  repId          = "Rep",
  plotId         = "Plot_id",
  rowNum         = "Range",
  colNum         = "Row",
  addCheck       = TRUE,
  checkGenotypes = c("B73"),
  timeFormat     = "%Y-%m-%d"     # match fake date format
)


## Fit spatial model for one or more traits
fit_D4_PH <- fitModels(
  TP       = phenoTP_D4,
  trait    = "PH",      # choose the trait you want to model
  engine   = "SpATS",
  useCheck = TRUE # splits checks as fixed, others random
)

## Plot spatial trend (optional visualization)
plot(fit_D4_PH, timePoints = "2024-07-01", plotType = "spatial", spaTrend = "raw")

#-----------loop for all traits---------------
traits_D4 <- c("GDDTS","GDDTA","PH","EH","NBR","SPAD1","SPAD2")

for(tr in traits_D4) {
  fit <- fitModels(
    TP = phenoTP_D4,
    trait = tr,
    engine = "SpATS",
    useCheck = TRUE
  )
  
  corr <- getCorrected(fit) |>
    rename(Plot_id = plotId) |>
    left_join(manifest_D4[,c("Plot_id","Genotype","Rep","Row","Range")], by="Plot_id")
  
  write.csv(corr, paste0("D4_", tr, "_2023.csv"), row.names = FALSE)
}

plot(fit_D4_SPAD2, timePoints = "2024-07-01", plotType = "spatial", spaTrend = "raw")
