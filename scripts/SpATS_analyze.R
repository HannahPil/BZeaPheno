install.packages("SpATS")
library(SpATS)

# Load Data Run 1

setwd("C:/Users/Hannah Pil/Documents/gemmalab/BZea/BZea phenotyping/BZeaPheno")

manifest_FT <- read.csv("B5_spats_manifest.csv")
map <- read.csv("B5_spats_map.csv", header = FALSE)

manifest_FT <- manifest_FT |>
  transform(
    Plot_id = factor(Plot_id),
    Row     = as.numeric(Row),
    Range   = as.numeric(Range),
    GDDTS   = as.numeric(GDDTS),
    GDDTA   = as.numeric(GDDTA),
    PH   = as.numeric(PH),
    EN   = as.numeric(EN),
    EH   = as.numeric(EH),
    NBR   = as.numeric(NBR),
    SPAD   = as.numeric(SPAD)
  )
#check number of columns and rows 

ncol(map)
nrow(map)
# SpATS model run

spats1 <- SpATS(response="EH",
                spatial=~SAP(Row, Range, nseg = c(2,4), degree = 3, pord = 2),
                 genotype = "Plot_id",
                 data = manifest_FT,
                 control = list(tolerance = 1e-3, maxit = 500),
                 genotype.as.random = TRUE
               )
print(spats1)

# Make some plots 

plot(spats1)

plot(SpATS::variogram(spats1))

#Make corrected values dataframe 

spatially_corrected_B5 <- as.data.frame(spatsmodel1$fitted)


#In run1_corrected_pheno add Genotype, Row, Column

run1_corrected_pheno[c("Genotype", "Row", "Column")] <- run1_pheno[c("Genotype", "Row", "Column")]

#save updated spatial run data

write.csv(run1_corrected_pheno, "C:/Users/zehta/Github/PTxB73-Cold-Experiment/Results/run1_corrected_pheno.csv", row.names = FALSE)

#######################################################################################################


### start a new BLUP‐accumulator df with one row per LineID
blup_df <- data.frame(spatsmodel1 = run1_pheno$Genotype, stringsAsFactors = FALSE)

# extract the BLUP vector
geno_names  <- spatsmodel1$terms$Genotype
geno_blups  <- spatsmodel1$coeff[geno_names]
blups_i     <- data.frame(LineID = Genotype,
                          BLUP   = as.numeric(geno_blups),
                          stringsAsFactors = FALSE)

#m = your spats model, LineID = genotype, geno_names = Sample names


########################################################################################
############################################################################################################################

