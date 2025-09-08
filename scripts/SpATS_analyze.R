# Install and load SpATS packages 

install.packages("SpATS")

library(SpATS)

# Load Data Run 1

run1_pheno <- read.csv("Data/Phenotype/Run1_Row_Column_AvgWeight.csv")
head(run1_pheno)

run1_map <- read.csv("Data/Chamber_Maps/Run1_GenoMap.csv")
head(run1_map)

# Checking number of columns and rows 

ncol(run1_map)
nrow(run1_map)

#switch A-H values to numerical values 

run1_pheno[run1_pheno=="A"] <- 8
run1_pheno[run1_pheno=="B"] <- 7
run1_pheno[run1_pheno=="C"] <- 6
run1_pheno[run1_pheno=="D"] <- 5
run1_pheno[run1_pheno=="E"] <- 4
run1_pheno[run1_pheno=="F"] <- 3
run1_pheno[run1_pheno=="G"] <- 2
run1_pheno[run1_pheno=="H"] <- 1

# changing Row to numerical values 

str(run1_pheno)
run1_pheno$Row <- as.numeric(as.character(run1_pheno$Row))

# SpATS model run

spatsmodel1 <- SpATS(response="Average_Weight",
                     spatial=~SAP(Row, Column, nseg = c(2,4), degree = 3, pord = 2),
                     genotype = "Genotype",
                     data = run1_pheno,
                     control = list(tolerance = 1e-3, maxit = 500),
                     genotype.as.random = TRUE
                     )
print(spatsmodel1)

# Make some plots 

plot(spatsmodel1)


plot(SpATS::variogram(spatsmodel1))

#Make corrected values dataframe 

run1_corrected_pheno <- as.data.frame(spatsmodel1$fitted)

# Add column called Run to run1_pheno and run1_corrected_pheno

run1_corrected_pheno$Run <- 1
run1_pheno$Run <- 1

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


#########################################################################################33
############################################################################################################################

# Load Data Run 2

run2_pheno <- read.csv("Data/Phenotype/Run2_Row_Column_AvgWeight.csv")
head(run2_pheno)

run2_map <- read.csv("Data/Chamber_Maps/Run2_GenoMap.csv")
head(run2_map)

# Checking number of columns and rows 

ncol(run2_map)
nrow(run2_map)

#switch A-H values to numerical values 

run2_pheno[run2_pheno=="A"] <- 8
run2_pheno[run2_pheno=="B"] <- 7
run2_pheno[run2_pheno=="C"] <- 6
run2_pheno[run2_pheno=="D"] <- 5
run2_pheno[run2_pheno=="E"] <- 4
run2_pheno[run2_pheno=="F"] <- 3
run2_pheno[run2_pheno=="G"] <- 2
run2_pheno[run2_pheno=="H"] <- 1

# changing Row to numerical values 

str(run2_pheno)
run2_pheno$Row <- as.numeric(as.character(run2_pheno$Row))

# SpATS model run

spatsmodel2 <- SpATS(response="Average_Weight",
                     spatial=~SAP(Row, Column, nseg = c(2,4), degree = 3, pord = 2),
                     genotype = "Genotype",
                     data = run2_pheno,
                     control = list(tolerance = 1e-3, maxit = 500),
                     genotype.as.random = TRUE
)
print(spatsmodel2)

# Make some plots 

plot(spatsmodel2)


plot(SpATS::variogram(spatsmodel2))

#Make corrected values dataframe 

run2_corrected_pheno <- as.data.frame(spatsmodel2$fitted)

# Add column called Run to run2_pheno and run1_corrected_pheno

run2_corrected_pheno$Run <- 2
run2_pheno$Run <- 2

#In run2_corrected_pheno add Genotype, Row, Column

run2_corrected_pheno[c("Genotype", "Row", "Column")] <- run2_pheno[c("Genotype", "Row", "Column")]

#save updated spatial run data

write.csv(run2_corrected_pheno, "C:/Users/zehta/Github/PTxB73-Cold-Experiment/Results/run2_corrected_pheno.csv", row.names = FALSE)


############################################################################################################################

# Load Data Run 3

run3_pheno <- read.csv("Data/Phenotype/Run3_Row_Column_AvgWeight.csv")
head(run3_pheno)

run3_map <- read.csv("Data/Chamber_Maps/Run3_GenoMap.csv")
head(run3_map)

# Checking number of columns and rows 

ncol(run3_map)
nrow(run3_map)

#switch A-H values to numerical values 

run3_pheno[run3_pheno=="A"] <- 8
run3_pheno[run3_pheno=="B"] <- 7
run3_pheno[run3_pheno=="C"] <- 6
run3_pheno[run3_pheno=="D"] <- 5
run3_pheno[run3_pheno=="E"] <- 4
run3_pheno[run3_pheno=="F"] <- 3
run3_pheno[run3_pheno=="G"] <- 2
run3_pheno[run3_pheno=="H"] <- 1

# changing Row to numerical values 

str(run3_pheno)
run3_pheno$Row <- as.numeric(as.character(run3_pheno$Row))

# SpATS model run

spatsmodel3 <- SpATS(response="Average_Weight",
                     spatial=~SAP(Row, Column, nseg = c(2,4), degree = 3, pord = 2),
                     genotype = "Genotype",
                     data = run3_pheno,
                     control = list(tolerance = 1e-3, maxit = 500),
                     genotype.as.random = TRUE
)
print(spatsmodel3)

# Make some plots 

plot(spatsmodel3)


plot(SpATS::variogram(spatsmodel3))

#Make corrected values dataframe 

run3_corrected_pheno <- as.data.frame(spatsmodel3$fitted)

# Add column called Run to run3_pheno and run3_corrected_pheno

run3_corrected_pheno$Run <- 3
run3_pheno$Run <- 3

#In run3_corrected_pheno add Genotype, Row, Column

run3_corrected_pheno[c("Genotype", "Row", "Column")] <- run3_pheno[c("Genotype", "Row", "Column")]

#save updated spatial run data

write.csv(run3_corrected_pheno, "C:/Users/zehta/Github/PTxB73-Cold-Experiment/Results/run3_corrected_pheno.csv", row.names = FALSE)


############################################################################################################################

# Load Data Run 6

run6_pheno <- read.csv("Data/Phenotype/Run6_Row_Column_AvgWeight.csv")
head(run6_pheno)

run6_map <- read.csv("Data/Chamber_Maps/Run6_GenoMap.csv")
head(run6_map)

# Checking number of columns and rows 

ncol(run6_map)
nrow(run6_map)

#switch A-H values to numerical values, there isa lower case g so added that in as well 

run6_pheno[run6_pheno=="A"] <- 8
run6_pheno[run6_pheno=="B"] <- 7
run6_pheno[run6_pheno=="C"] <- 6
run6_pheno[run6_pheno=="D"] <- 5
run6_pheno[run6_pheno=="E"] <- 4
run6_pheno[run6_pheno=="F"] <- 3
run6_pheno[run6_pheno=="G"] <- 2
run6_pheno[run6_pheno=="g"] <- 2
run6_pheno[run6_pheno=="H"] <- 1

# changing Row to numerical values 

str(run6_pheno)
run6_pheno$Row <- as.numeric(as.character(run6_pheno$Row))

# SpATS model run

spatsmodel6 <- SpATS(response="Average_Weight",
                     spatial=~SAP(Row, Column, nseg = c(2,4), degree = 3, pord = 2),
                     genotype = "Genotype",
                     data = run6_pheno,
                     control = list(tolerance = 1e-3, maxit = 500),
                     genotype.as.random = TRUE
)
print(spatsmodel6)

# Make some plots 

plot(spatsmodel6)


plot(SpATS::variogram(spatsmodel6))

#Make corrected values dataframe 

run6_corrected_pheno <- as.data.frame(spatsmodel6$fitted)

# Add column called Run to run6_pheno and run6_corrected_pheno

run6_corrected_pheno$Run <- 6
run6_pheno$Run <- 6

#In run6_corrected_pheno add Genotype, Row, Column

run6_corrected_pheno[c("Genotype", "Row", "Column")] <- run6_pheno[c("Genotype", "Row", "Column")]

#save updated spatial run data

write.csv(run6_corrected_pheno, "C:/Users/zehta/Github/PTxB73-Cold-Experiment/Results/run6_corrected_pheno.csv", row.names = FALSE)


############################################################################################################################

# Load Data Run 4

run4_pheno <- read.csv("Data/Phenotype/Run4_Row_Column_AvgWeight.csv")
head(run4_pheno)

run4_map <- read.csv("Data/Chamber_Maps/Run4_GenoMap.csv")
head(run4_map)

# Checking number of columns and rows 

ncol(run4_map)
nrow(run4_map)

#switch A-H values to numerical values, there isa lower case g so added that in as well 

run4_pheno[run4_pheno=="A"] <- 8
run4_pheno[run4_pheno=="B"] <- 7
run4_pheno[run4_pheno=="C"] <- 6
run4_pheno[run4_pheno=="D"] <- 5
run4_pheno[run4_pheno=="E"] <- 4
run4_pheno[run4_pheno=="F"] <- 3
run4_pheno[run4_pheno=="G"] <- 2
run4_pheno[run4_pheno=="g"] <- 2
run4_pheno[run4_pheno=="H"] <- 1

# changing Row to numerical values 

str(run4_pheno)
run4_pheno$Row <- as.numeric(as.character(run4_pheno$Row))

# SpATS model run

spatsmodel4 <- SpATS(response="Average_Weight",
                     spatial=~SAP(Row, Column, nseg = c(2,4), degree = 3, pord = 2),
                     genotype = "Genotype",
                     data = run4_pheno,
                     control = list(tolerance = 1e-3, maxit = 500),
                     genotype.as.random = TRUE
)
print(spatsmodel4)

# Make some plots 

plot(spatsmodel4)


plot(SpATS::variogram(spatsmodel4))

#Make corrected values dataframe 

run4_corrected_pheno <- as.data.frame(spatsmodel4$fitted)

# Add column called Run to run4_pheno and run4_corrected_pheno

run4_corrected_pheno$Run <- 4
run4_pheno$Run <- 4

#In run4_corrected_pheno add Genotype, Row, Column

run4_corrected_pheno[c("Genotype", "Row", "Column")] <- run4_pheno[c("Genotype", "Row", "Column")]

#save updated spatial run data

write.csv(run4_corrected_pheno, "C:/Users/zehta/Github/PTxB73-Cold-Experiment/Results/run4_corrected_pheno.csv", row.names = FALSE)



############################################################################################################################

# Load Data Run 5

run5_pheno <- read.csv("Data/Phenotype/Run5_Row_Column_AvgWeight.csv")
head(run5_pheno)

run5_map <- read.csv("Data/Chamber_Maps/Run5_GenoMap.csv")
head(run5_map)

# Checking number of columns and rows 

ncol(run5_map)
nrow(run5_map)

#switch A-H values to numerical values, there isa lower case g so added that in as well 

run5_pheno[run5_pheno=="A"] <- 8
run5_pheno[run5_pheno=="B"] <- 7
run5_pheno[run5_pheno=="C"] <- 6
run5_pheno[run5_pheno=="D"] <- 5
run5_pheno[run5_pheno=="E"] <- 4
run5_pheno[run5_pheno=="F"] <- 3
run5_pheno[run5_pheno=="G"] <- 2
run5_pheno[run5_pheno=="g"] <- 2
run5_pheno[run5_pheno=="H"] <- 1

# changing Row to numerical values 

str(run5_pheno)
run5_pheno$Row <- as.numeric(as.character(run5_pheno$Row))

# SpATS model run

spatsmodel5 <- SpATS(response="Average_Weight",
                     spatial=~SAP(Row, Column, nseg = c(2,4), degree = 3, pord = 2),
                     genotype = "Genotype",
                     data = run5_pheno,
                     control = list(tolerance = 1e-3, maxit = 500),
                     genotype.as.random = TRUE
)
print(spatsmodel5)

# Make some plots 

plot(spatsmodel5)


plot(SpATS::variogram(spatsmodel5))

#Make corrected values dataframe 

run5_corrected_pheno <- as.data.frame(spatsmodel5$fitted)

# Add column called Run to run5_pheno and run5_corrected_pheno

run5_corrected_pheno$Run <- 5
run5_pheno$Run <- 5

#In run5_corrected_pheno add Genotype, Row, Column

run5_corrected_pheno[c("Genotype", "Row", "Column")] <- run5_pheno[c("Genotype", "Row", "Column")]

#save updated spatial run data

write.csv(run5_corrected_pheno, "C:/Users/zehta/Github/PTxB73-Cold-Experiment/Results/run5_corrected_pheno.csv", row.names = FALSE)


?SpATs::SpAts()


