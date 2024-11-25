setwd("C:/Users/mosuw/Desktop/Caivil")
getwd()

#Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")
BiocManager::install("DESeq2")
BiocManager::install("data.table")

#Libraries
library(tximport)
library(GenomicFeatures)
library(readr)
library(DESeq2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)

#reading in count date
gct_data <- fread("GSE110137_counts.gct", skip = 2)

#setting expression matrix from count data
gene_names <- gct_data$NAME
expression_matrix <- as.matrix(gct_data[, -c(1, 2), with = FALSE])
rownames(expression_matrix) <- gene_names



library(GEOquery) # was not working look below
library(airway)

# Visit the Rtools for Windows page. Download the version of Rtools that matches your R version (since you are using R 4.4, download Rtools44).
system('g++ --version')
install.packages("BiocManager")
aBiocManager::install("GEOquery")

library(GEOquery)

#get metadata
gse <- getGEO(GEO = "GSE110137", GSEMatrix = T)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

#reading in count date
gct_data <- fread("GSE110137_counts.gct", skip = 2)

#find the positions of sample and MOA
data.frame(Position = 1:ncol(metadata), Column_Name = colnames(metadata))

MOA_samples <- metadata %>%
  select(1,12) %>%
  mutate(characteristics_ch1.2 = gsub("moa:", "", characteristics_ch1.2)) 


#filter out not available data sets
gct_data_f <- gct_data %>%
  filter(Description != "-") %>%
  filter(Description != "not available") %>%
  distinct()

#checking how many were filtered out
nrow(gct_data)
nrow(gct_data_f)

#check how many colomuns i have in each position
data.frame(Position = 1:ncol(gct_data_f), Column_Name = colnames(gct_data_f))


########################TPM normaisation#####################
# Load necessary library
library(dplyr)

# Extract gene lengths and sample data
g
ene_lengths <- gct_data_f[[1]] # First column contains gene lengths
counts_data <- gct_data_f[, 3:59] # Columns 3 to 59 contain the counts data

# Calculate RPK (Reads Per Kilobase)
rpk <- counts_data / (gene_lengths / 1000)

# Calculate per-sample scaling factor (sum of RPK values)
scaling_factors <- colSums(rpk)

# Calculate TPM
tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6

# Recombine TPM with original gene length and description
tpm_data <- cbind(gct_data_f[, 1:2], tpm)

# View resulting TPM data
head(tpm_data)

# Calculate the sum of each column (samples from column 3 to 59)
column_sums <- colSums(gct_data_f[, 3:59])

print(column_sums)

#remove NAME column
tpm_data <- tpm_data[, -1]


duplicates_sum <- metadata %>%
  dplyr::select(1, 8, 12) %>%
  mutate(characteristics_ch1.2 = gsub("moa:","",characteristics_ch1.2))

###########################renaming colunms of sample to MoA###########################
# Create a named vector mapping samples to MOA
samples <- c("PDD_P2_01", "PDD_P2_02", "PDD_P2_03", "PDD_P2_04", "PDD_P2_05", "PDD_P2_06", 
             "PDD_P2_07", "PDD_P2_09", "PDD_P2_10", "PDD_P2_11", "PDD_P2_12", "PDD_P2_14", 
             "PDD_P2_15", "PDD_P2_18", "PDD_P2_20", "PDD_P2_21", "PDD_P2_23", "PDD_P2_24", 
             "PDD_P2_26", "PDD_P2_27", "PDD_P2_28", "PDD_P2_29", "PDD_P2_30", "PDD_P2_31", 
             "PDD_P2_32", "PDD_P2_33", "PDD_P2_35", "PDD_P2_36", "PDD_P2_37", "PDD_P2_38", 
             "PDD_P2_40", "PDD_P2_41", "PDD_P2_44", "PDD_P2_46", "PDD_P2_47", "PDD_P2_49", 
             "PDD_P2_50", "PDD_P2_52", "PDD_P2_53", "PDD_P2_54", "PDD_P2_55", "PDD_P2_56", 
             "PDD_P2_57", "PDD_P2_58", "PDD_P2_59", "PDD_P2_61", "PDD_P2_62", "PDD_P2_63", 
             "PDD_P2_64", "PDD_P2_66", "PDD_P2_67", "PDD_P2_70", "PDD_P2_72", "PDD_P2_73", 
             "PDD_P2_75", "PDD_P2_76", "PDD_P2_78")

moa <- c("DMSO", "DMSO", "DNA damage", "cell wall synthesis inhibitor", "Protein synthesis", 
         "DNA replication inhibitor", "Protein synthesis", "Protein synthesis", "DNA replication inhibitor", 
         "cell wall synthesis inhibitor", "Membrane perturbation", "cell wall synthesis inhibitor", 
         "Membrane perturbation", "cell wall synthesis inhibitor / lipoprotein", "B-01(act)", 
         "DNA replication inhibitor", "cell wall synthesis inhibitor / lipoprotein", "DNA replication inhibitor", 
         "Folic Acid synthesis inhibitor", "DMSO", "DMSO", "DNA damage", "cell wall synthesis inhibitor", 
         "Protein synthesis", "DNA replication inhibitor", "Protein synthesis", "Protein synthesis", 
         "Protein synthesis", "DNA replication inhibitor", "cell wall synthesis inhibitor", 
         "Membrane perturbation", "Membrane perturbation", "Membrane perturbation", 
         "cell wall synthesis inhibitor / lipoprotein", "B-01(act)", "DNA replication inhibitor", 
         "cell wall synthesis inhibitor / lipoprotein", "DNA replication inhibitor", "Folic Acid synthesis inhibitor", 
         "DMSO", "DMSO", "DNA damage", "cell wall synthesis inhibitor", "Protein synthesis", 
         "DNA replication inhibitor", "Protein synthesis", "Protein synthesis", "Protein synthesis", 
         "DNA replication inhibitor", "cell wall synthesis inhibitor", "Membrane perturbation", 
         "cell wall synthesis inhibitor", "Membrane perturbation", "Membrane perturbation", 
         "cell wall synthesis inhibitor / lipoprotein", "B-01(act)", "DNA replication inhibitor")

# Create a named vector mapping samples to MOA
names(moa) <- samples

# Assuming 'tpm_data' is your data frame, and samples are in columns 3 to 59, we will rename those columns based on the MOA
colnames(tpm_data)[3:59] <- moa[match(colnames(tpm_data)[3:59], samples)]

# View the updated column names
head(colnames(tpm_data)[3:59])

#check how many colomuns i have in each position
data.frame(Position = 1:ncol(tpm_data), Column_Name = colnames(tpm_data))

#################################################save as CSV.file###################################################
#save tpm_data file to personal computer
write.csv(tpm_data, 
          "~/bioinfo_Masters_Project/Project/understanding MOA/ML/tpm_data.csv", 
          row.names = T)
#########################################import CSV.file after removing "NAME" column#######################
library(readr)
tpm_data <- read_csv("tpm_data.csv")

#check how many colomuns i have in each position
data.frame(Position = 1:ncol(tpm_data), Column_Name = colnames(tpm_data))
################################################rows to be columns and columns to rows ####################################################################
############################Random Forest########################
transposed_data_m <- t(tpm_data)
write.csv(transposed_data_m, 
          "~/bioinfo_Masters_Project/Project/understanding MOA/ML/transposed_data_m.csv", 
          row.names = T)

transposed_data <- read_csv("transposed_data.csv")


write.csv(T_gct_data, 
          "~/bioinfo_Masters_Project/Project/understanding MOA/ML/T_gct_data.csv", 
          row.names = T)
T_gct_data <- t(tpm_data)
