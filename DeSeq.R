#######################################Differential gene expression#################################

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

# finding individual drug treatments to make a dataset for each.
extraction <- metadata %>%
  select(8)
compounds <- data.frame(extraction$source_name_ch1)
compounds <- compounds %>%
  distinct() %>%
  filter(extraction.source_name_ch1 != "DMSO")

#make compound libraries metadata to know which compound correlates to which sample in the counts data
DMSO <- metadata %>%
  filter(source_name_ch1 == "DMSO")%>%
  select(1)

LolCDE <- metadata %>%
  filter(source_name_ch1 == "AZ-LolCDE") %>%
  select(1)

B01 <- metadata %>%
  filter(source_name_ch1 == "B-01") %>%
  select(1)

Ceftriaxone <- metadata %>%
  filter(source_name_ch1 == "Ceftriaxone") %>%
  select(1)

Chloramphenicol <- metadata %>%
  filter(source_name_ch1 == "Chloramphenicol") %>%
  select(1)

Ciprofloxacin <- metadata %>%
  filter(source_name_ch1 == "Ciprofloxacin") %>%
  select(1)

Clarithromycin <- metadata %>%
  filter(source_name_ch1 == "Clarithromycin") %>%
  select(1)

Colistin <- metadata %>%
  filter(source_name_ch1 == "Colistin") %>%
  select(1)

Doxycycline <- metadata %>%
  filter(source_name_ch1 == "Doxycycline") %>%
  select (1)

Globomycin <- metadata %>%
  filter(source_name_ch1 == "Globomycin") %>%
  select (1)
	
Levofloxacin <- metadata %>%
  filter(source_name_ch1 == "Levofloxacin") %>%
  select (1)

Mecillinam <- metadata %>%
  filter(source_name_ch1 == "Mecillinam") %>%
  select (1)

Meropenem <- metadata %>%
  filter(source_name_ch1 == "Meropenem") %>%
  select (1)

Nitrofurantoin <- metadata %>%
  filter(source_name_ch1 == "Nitrofurantoin") %>%
  select (1)

Nitroxolin <- metadata %>%
  filter(source_name_ch1 == "Nitroxolin") %>%
  select (1)

Norfloxacin <- metadata %>%
  filter(source_name_ch1 == "Norfloxacin") %>%
  select (1)

PolymyxineB <- metadata %>%
  filter(source_name_ch1 == "PolymyxineB") %>%
  select (1)

Trimethoprim <- metadata %>%
  filter(source_name_ch1 == "Trimethoprim") %>%
  select (1)

#data frame containing names and postions of columns 
data.frame(Position = 1:ncol(expression_matrix), Column_Name = colnames(expression_matrix))

# Separating matrix to each compounds: [, c(1,2,20,21,39,40)] is the controls
LolCDE_matrix <- expression_matrix[, c(1,2,20,21,39,40,52,14,33)]
B01_matrix <- expression_matrix[, c(1,2,20,21,39,40,53,15,34)]
Ceftriaxone_matrix <- expression_matrix[, c(1,2,20,21,39,40,48,10,29)]
Chloramphenicol_matrix  <- expression_matrix[, c(1,2,20,21,39,40,43,5,24)]
Ciprofloxacin_matrix <- expression_matrix[, c(1,2,20,21,39,40,54,16,35)]
Clarithromycin_matrix <- expression_matrix[, c(1,2,20,21,39,40,45,7,26)]
Colistin_matrix <- expression_matrix[, c(1,2,20,21,39,40,51,13,32)]
Doxycycline_matrix <- expression_matrix[, c(1,2,20,21,39,40,46,8,27)]
Globomycin_matrix <- expression_matrix[, c(1,2,20,21,39,40,55,17,36)]
Levofloxacin_matrix <- expression_matrix[, c(1,2,20,21,39,40,47,9,28)]
Mecillinam_matrix <- expression_matrix[, c(1,2,20,21,39,40,50,12,31)]
Meropenem_matrix <- expression_matrix[, c(1,2,20,21,39,40,42,4,23)]
Nitrofurantoin_matrix <- expression_matrix[, c(1,2,20,21,39,40,41,3,22)]
Nitroxolin_matrix <- expression_matrix[, c(1,2,20,21,39,40,44,6,25)]
Norfloxacin_matrix <- expression_matrix[, c(1,2,20,21,39,40,56,18,37)]
PolymyxineB_matrix <- expression_matrix[, c(1,2,20,21,39,40,49,11,30)]
Trimethoprim_matrix <- expression_matrix[, c(1,2,20,21,39,40,57,19,38)]


#save metadata file to personal computer
write.csv(metadata, "~/bioinfo_Masters_Project/Project/understanding MOA/ML", row.names = FALSE)

#differentiating treatment and control samples 
sample_info_LolCDE <- data.frame(
  row.names = colnames(LolCDE_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_B01 <- data.frame(
  row.names = colnames(B01_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Ceftriaxone <- data.frame(
  row.names = colnames(Ceftriaxone_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Chloramphenicol <- data.frame(
  row.names = colnames(Chloramphenicol_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Ciprofloxacin <- data.frame(
  row.names = colnames(Ciprofloxacin_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Clarithromycin <- data.frame(
  row.names = colnames(Clarithromycin_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Colistin<- data.frame(
  row.names = colnames(Colistin_matrix ),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Doxycycline <- data.frame(
  row.names = colnames(Doxycycline_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Globomycin <- data.frame(
  row.names = colnames(Globomycin_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Levofloxacin<- data.frame(
  row.names = colnames(Levofloxacin_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Mecillinam <- data.frame(
  row.names = colnames(Mecillinam_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Meropenem<- data.frame(
  row.names = colnames(Meropenem_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Nitrofurantoin <- data.frame(
  row.names = colnames(Nitrofurantoin_matrix ),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Nitroxolin<- data.frame(
  row.names = colnames(Nitroxolin_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Norfloxacin <- data.frame(
  row.names = colnames(Norfloxacin_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_PolymyxineB <- data.frame(
  row.names = colnames(PolymyxineB_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)
sample_info_Trimethoprim <- data.frame(
  row.names = colnames(Trimethoprim_matrix),
  condition = c("control", "control", "control","control", "control", "control","treated","treated","treated")
)



#converted to a round matrix for DESeq2)
LolCDE_matrix <- round(LolCDE_matrix)
B01_matrix <- round(B01_matrix)
Ceftriaxone_matrix <- round(Ceftriaxone_matrix)
Chloramphenicol_matrix <- round(Chloramphenicol_matrix)
Ciprofloxacin_matrix <- round(Ciprofloxacin_matrix)
Clarithromycin_matrix <- round(Clarithromycin_matrix)
Colistin_matrix <- round(Colistin_matrix)
Doxycycline_matrix <- round(Doxycycline_matrix)
Globomycin_matrix <- round(Globomycin_matrix)
Levofloxacin_matrix <- round(Levofloxacin_matrix)
Mecillinam_matrix <- round(Mecillinam_matrix)
Meropenem_matrix <- round(Meropenem_matrix)
Nitrofurantoin_matrix <- round(Nitrofurantoin_matrix )
Nitroxolin_matrix <- round(Nitroxolin_matrix)
Norfloxacin_matrix<- round(Norfloxacin_matrix)
PolymyxineB_matrix <- round(PolymyxineB_matrix)
Trimethoprim_matrix <- round(Trimethoprim_matrix)

# Verify the conversion
head(LolCDE_matrix)
head(sample_info_LolCDE)

# Create DESeqDataSet object
LolCDE_dds <- DESeqDataSetFromMatrix(countData = LolCDE_matrix,
                                  colData = sample_info_LolCDE,
                                  design = ~ condition)
LolCDE_dds<- DESeq(LolCDE_dds)
results_LolCDE <- results(LolCDE_dds)
###########################################################
B01_dds <- DESeqDataSetFromMatrix(countData = B01_matrix,
                                  colData = sample_info_B01,
                                  design = ~ condition)
B01_dds<- DESeq(B01_dds)
results_B01 <- results(B01_dds)
###########################################################
Ceftriaxone_dds <- DESeqDataSetFromMatrix(countData = Ceftriaxone_matrix,
                              colData = sample_info_Ceftriaxone,
                              design = ~ condition)
Ceftriaxone_dds<- DESeq(Ceftriaxone_dds)
results_Ceftriaxone <- results(Ceftriaxone_dds)
############################################################
Chloramphenicol_dds <- DESeqDataSetFromMatrix(countData = Chloramphenicol_matrix,
                                  colData = sample_info_Chloramphenicol,
                                  design = ~ condition)
Chloramphenicol_dds<- DESeq(Chloramphenicol_dds)
results_Chloramphenicol <- results(Chloramphenicol_dds)
############################################################
Ciprofloxacin_dds <- DESeqDataSetFromMatrix(countData = Ciprofloxacin_matrix,
                                  colData = sample_info_Ciprofloxacin,
                                  design = ~ condition)
Ciprofloxacin_dds<- DESeq(Ciprofloxacin_dds)
results_Ciprofloxacin <- results(Ciprofloxacin_dds)
###########################################################
Clarithromycin_dds <- DESeqDataSetFromMatrix(countData = Clarithromycin_matrix,
                                  colData = sample_info_Clarithromycin,
                                  design = ~ condition)
Clarithromycin_dds<- DESeq(Clarithromycin_dds)
results_Clarithromycin <- results(Clarithromycin_dds)
###############################################################
Colistin_dds <- DESeqDataSetFromMatrix(countData = Colistin_matrix,
                                  colData = sample_info_Colistin,
                                  design = ~ condition)
Colistin_dds<- DESeq(Colistin_dds)
results_Colistin <- results(Colistin_dds)
#################################################################
Doxycycline_dds <- DESeqDataSetFromMatrix(countData = Doxycycline_matrix,
                                  colData = sample_info_Doxycycline,
                                  design = ~ condition)
Doxycycline_dds<- DESeq(Doxycycline_dds)
results_Doxycycline <- results(Doxycycline_dds)
#################################################################
Globomycin_dds <- DESeqDataSetFromMatrix(countData = Globomycin_matrix,
                                  colData = sample_info_Globomycin,
                                  design = ~ condition)
Globomycin_dds<- DESeq(Globomycin_dds)
results_Globomycin <- results(Globomycin_dds)
#################################################################
Levofloxacin_dds <- DESeqDataSetFromMatrix(countData = Levofloxacin_matrix,
                                  colData = sample_info_Levofloxacin,
                                  design = ~ condition)
Levofloxacin_dds<- DESeq(Levofloxacin_dds)
results_Levofloxacin <- results(Levofloxacin_dds)
###################################################################
Mecillinam_dds <- DESeqDataSetFromMatrix(countData = Mecillinam_matrix,
                                  colData = sample_info_Mecillinam,
                                  design = ~ condition)
Mecillinam_dds<- DESeq(Mecillinam_dds)
results_Mecillinam <- results(Mecillinam_dds)
#################################################################
Meropenem_dds <- DESeqDataSetFromMatrix(countData = Meropenem_matrix,
                                  colData = sample_info_Meropenem,
                                  design = ~ condition)
Meropenem_dds<- DESeq(Meropenem_dds)
results_Meropenem <- results(Meropenem_dds)
##################################################################
Nitrofurantoin_dds <- DESeqDataSetFromMatrix(countData = Nitrofurantoin_matrix,
                                  colData = sample_info_Nitrofurantoin,
                                  design = ~ condition)
Nitrofurantoin_dds<- DESeq(Nitrofurantoin_dds)
results_Nitrofurantoin <- results(Nitrofurantoin_dds)
##################################################################
Norfloxacin_dds <- DESeqDataSetFromMatrix(countData = Norfloxacin_matrix,
                                  colData = sample_info_Norfloxacin,
                                  design = ~ condition)
Norfloxacin_dds<- DESeq(Norfloxacin_dds)
results_Norfloxacin <- results(Norfloxacin_dds)
#####################################################################
Nitroxolin_dds <- DESeqDataSetFromMatrix(countData = Nitroxolin_matrix,
                                  colData = sample_info_Nitroxolin,
                                  design = ~ condition)
Nitroxolin_dds<- DESeq(Nitroxolin_dds)
results_Nitroxolin <- results(Nitroxolin_dds)
####################################################################
PolymyxineB_dds <- DESeqDataSetFromMatrix(countData = PolymyxineB_matrix,
                                  colData = sample_info_PolymyxineB,
                                  design = ~ condition)
PolymyxineB_dds<- DESeq(PolymyxineB_dds)
results_PolymyxineB <- results(PolymyxineB_dds)
#######################################################################
Trimethoprim_dds <- DESeqDataSetFromMatrix(countData = Trimethoprim_matrix,
                                          colData = sample_info_Trimethoprim,
                                          design = ~ condition)
Trimethoprim_dds<- DESeq(Trimethoprim_dds)
results_Trimethoprim <- results(Trimethoprim_dds)

# This will display the adjusted p-values for each gene
#sig_res <- results[!is.na(results$padj) & results$padj < 0.05, ] (source code)
#head (sig_res)

LolCDE_sig_res <- results_LolCDE[!is.na(results_LolCDE$padj) & results_LolCDE$padj < 0.05, ]
B01_sig_res <- results_B01[!is.na(results_B01$padj) & results_B01$padj < 0.05, ]
Ceftriaxone_sig_res <- results_Ceftriaxone[!is.na(results_Ceftriaxone$padj) & results_Ceftriaxone$padj < 0.05, ]
Chloramphenicol_sig_res <- results_Chloramphenicol[!is.na(results_Chloramphenicol$padj) & results_Chloramphenicol$padj < 0.05, ]
Ciprofloxacin_sig_res <- results_Ciprofloxacin[!is.na(results_Ciprofloxacin$padj) & results_Ciprofloxacin$padj < 0.05, ]
Clarithromycin_sig_res <- results_Clarithromycin[!is.na(results_Clarithromycin$padj) & results_Clarithromycin$padj < 0.05, ]
Colistin_sig_res <- results_Colistin[!is.na(results_Colistin$padj) & results_Colistin$padj < 0.05, ]
Doxycycline_sig_res <- results_Doxycycline[!is.na(results_Doxycycline$padj) & results_Doxycycline$padj < 0.05, ]
Globomycin_sig_res <- results_Globomycin[!is.na(results_Globomycin$padj) & results_Globomycin$padj < 0.05, ]
Levofloxacin_sig_res <- results_Levofloxacin[!is.na(results_Levofloxacin$padj) & results_Levofloxacin$padj < 0.05, ]
Mecillinam_sig_res <- results_Mecillinam[!is.na(results_Mecillinam$padj) & results_Mecillinam$padj < 0.05, ]
Meropenem_sig_res <- results_Meropenem[!is.na(results_Meropenem$padj) & results_Meropenem$padj < 0.05, ]
Nitrofurantoin_sig_res <- results_Nitrofurantoin[!is.na(results_Nitrofurantoin$padj) & results_Nitrofurantoin$padj < 0.05, ]
Norfloxacin_sig_res <- results_Norfloxacin[!is.na(results_Norfloxacin$padj) & results_Norfloxacin$padj < 0.05, ]
Nitroxolin_sig_res <- results_Nitroxolin[!is.na(results_Nitroxolin$padj) & results_Nitroxolin$padj < 0.05, ]
PolymyxineB_sig_res <- results_PolymyxineB[!is.na(results_PolymyxineB$padj) & results_PolymyxineB$padj < 0.05, ]
Trimethoprim_sig_res <- results_Trimethoprim[!is.na(results_Trimethoprim$padj) & results_Trimethoprim$padj < 0.05, ]

# Convert results to a data frame
#results_df <- as.data.frame(results) (source code)

LolCDE_df <- as.data.frame(results_LolCDE)
B01_df <- as.data.frame(results_B01)
Ceftriaxone_df <- as.data.frame(results_Ceftriaxone)
Chloramphenicol_df <- as.data.frame(results_Chloramphenicol)
Ciprofloxacin_df <- as.data.frame(results_Ciprofloxacin)
Clarithromycin_df <- as.data.frame(results_Clarithromycin)
Colistin_df <- as.data.frame(results_Colistin)
Doxycycline_df <- as.data.frame(results_Doxycycline)
Globomycin_df <- as.data.frame(results_Globomycin)
Levofloxacin_df <- as.data.frame(results_Levofloxacin)
Mecillinam_df <- as.data.frame(results_Mecillinam)
Meropenem_df <- as.data.frame(results_Meropenem)
Nitrofurantoin_df <- as.data.frame(results_Nitrofurantoin)
Norfloxacin_df <- as.data.frame(results_Norfloxacin)
Nitroxolin_df <- as.data.frame(results_Nitroxolin)
PolymyxineB_df <- as.data.frame(results_PolymyxineB)
Trimethoprim_df <- as.data.frame(results_Trimethoprim)

#filter out " -" and "not available"
#results_df_filter <-  results_df%>%
  #filter(pvalue != "NA") (sources code)

LolCDE_df_filter <- LolCDE_df %>%
  filter(!is.na(pvalue))

B01_df_filter <- B01_df %>%
  filter(!is.na(pvalue))

Ceftriaxone_df_filter <- Ceftriaxone_df %>%
  filter(!is.na(pvalue))

Chloramphenicol_df_filter <- Chloramphenicol_df %>%
  filter(!is.na(pvalue))

Ciprofloxacin_df_filter <- Ciprofloxacin_df %>%
  filter(!is.na(pvalue))

Clarithromycin_df_filter <- Clarithromycin_df %>%
  filter(!is.na(pvalue))

Colistin_df_filter <- Colistin_df %>%
  filter(!is.na(pvalue))

Doxycycline_df_filter <- Doxycycline_df %>%
  filter(!is.na(pvalue))

Globomycin_df_filter <- Globomycin_df %>%
  filter(!is.na(pvalue))

Levofloxacin_df_filter <- Levofloxacin_df %>%
  filter(!is.na(pvalue))

Mecillinam_df_filter <- Mecillinam_df %>%
  filter(!is.na(pvalue))

Meropenem_df_filter <- Meropenem_df %>%
  filter(!is.na(pvalue))

Nitrofurantoin_df_filter <- Nitrofurantoin_df %>%
  filter(!is.na(pvalue))

Norfloxacin_df_filter <- Norfloxacin_df %>%
  filter(!is.na(pvalue))

Nitroxolin_df_filter <- Nitroxolin_df %>%
  filter(!is.na(pvalue))

PolymyxineB_df_filter <- PolymyxineB_df %>%
  filter(!is.na(pvalue))

Trimethoprim_df_filter <- Trimethoprim_df %>%
  filter(!is.na(pvalue))

#number of rows
#nrow(results_df)
#nrow(results_df_filter) (source code)

nrow(LolCDE_df)  # Unfiltered
nrow(LolCDE_df_filter)  # Filtered (example but not required)

###############################################################################################################################################################
###############################################WE ARE USING THIS OUTPUT#############################################################
##############################################################################################################################
# Add a column for significance (adjust as needed)
#results_df_filter$significant <- ifelse(results_df_filter$padj < 0.05 & abs(results_df_filter$log2FoldChange) > 1, "Significant", "Not Significant")

LolCDE_df$significant <- ifelse(LolCDE_df$padj < 0.05 & abs(LolCDE_df$log2FoldChange) > 1, "Significant", "Not Significant")

B01_df$significant <- ifelse(B01_df$padj < 0.05 & abs(B01_df$log2FoldChange) > 1, "Significant", "Not Significant")

Ceftriaxone_df$significant <- ifelse(Ceftriaxone_df$padj < 0.05 & abs(Ceftriaxone_df$log2FoldChange) > 1, "Significant", "Not Significant")

Chloramphenicol_df$significant <- ifelse(Chloramphenicol_df$padj < 0.05 & abs(Chloramphenicol_df$log2FoldChange) > 1, "Significant", "Not Significant")

Ciprofloxacin_df$significant <- ifelse(Ciprofloxacin_df$padj < 0.05 & abs(Ciprofloxacin_df$log2FoldChange) > 1, "Significant", "Not Significant")

Clarithromycin_df$significant <- ifelse(Clarithromycin_df$padj < 0.05 & abs(Clarithromycin_df$log2FoldChange) > 1, "Significant", "Not Significant")

Colistin_df$significant <- ifelse(Colistin_df$padj < 0.05 & abs(Colistin_df$log2FoldChange) > 1, "Significant", "Not Significant")

Doxycycline_df$significant <- ifelse(Doxycycline_df$padj < 0.05 & abs(Doxycycline_df$log2FoldChange) > 1, "Significant", "Not Significant")

Globomycin_df$significant <- ifelse(Globomycin_df$padj < 0.05 & abs(Globomycin_df$log2FoldChange) > 1, "Significant", "Not Significant")

Levofloxacin_df$significant <- ifelse(Levofloxacin_df$padj < 0.05 & abs(Levofloxacin_df$log2FoldChange) > 1, "Significant", "Not Significant")

Mecillinam_df$significant <- ifelse(Mecillinam_df$padj < 0.05 & abs(Mecillinam_df$log2FoldChange) > 1, "Significant", "Not Significant")

Meropenem_df$significant <- ifelse(Meropenem_df$padj < 0.05 & abs(Meropenem_df$log2FoldChange) > 1, "Significant", "Not Significant")

Nitrofurantoin_df$significant <- ifelse(Nitrofurantoin_df$padj < 0.05 & abs(Nitrofurantoin_df$log2FoldChange) > 1, "Significant", "Not Significant")

Norfloxacin_df$significant <- ifelse(Norfloxacin_df$padj < 0.05 & abs(Norfloxacin_df$log2FoldChange) > 1, "Significant", "Not Significant")

Nitroxolin_df$significant <- ifelse(Nitroxolin_df$padj < 0.05 & abs(Nitroxolin_df$log2FoldChange) > 1, "Significant", "Not Significant")

PolymyxineB_df$significant <- ifelse(PolymyxineB_df$padj < 0.05 & abs(PolymyxineB_df$log2FoldChange) > 1, "Significant", "Not Significant")

Trimethoprim_df$significant <- ifelse(Trimethoprim_df$padj < 0.05 & abs(Trimethoprim_df$log2FoldChange) > 1, "Significant", "Not Significant")


#filter only sig values (not yet tho)
#Sig_df <-  results_df_filter%>%
  #filter(significant == "Significant")

LolCDE_Sig_df <- LolCDE_df %>%
  filter(significant == "Significant")

B01_Sig_df <- B01_df %>%
  filter(significant == "Significant")

Ceftriaxone_Sig_df <- Ceftriaxone_df %>%
  filter(significant == "Significant")

Chloramphenicol_Sig_df <- Chloramphenicol_df %>%
  filter(significant == "Significant")

Ciprofloxacin_Sig_df <- Ciprofloxacin_df %>%
  filter(significant == "Significant")

Clarithromycin_Sig_df <- Clarithromycin_df %>%
  filter(significant == "Significant")

Colistin_Sig_df <- Colistin_df %>%
  filter(significant == "Significant")

Doxycycline_Sig_df <- Doxycycline_df %>%
  filter(significant == "Significant")

Globomycin_Sig_df <- Globomycin_df %>%
  filter(significant == "Significant")

Levofloxacin_Sig_df <- Levofloxacin_df %>%
  filter(significant == "Significant")

Mecillinam_Sig_df <- Mecillinam_df %>%
  filter(significant == "Significant")

Meropenem_Sig_df <- Meropenem_df %>%
  filter(significant == "Significant")

Nitrofurantoin_Sig_df <- Nitrofurantoin_df %>%
  filter(significant == "Significant")

Norfloxacin_Sig_df <- Norfloxacin_df %>%
  filter(significant == "Significant")

Nitroxolin_Sig_df <- Nitroxolin_df %>%
  filter(significant == "Significant")

PolymyxineB_Sig_df <- PolymyxineB_df %>%
  filter(significant == "Significant")

Trimethoprim_Sig_df <- Trimethoprim_df %>%
  filter(significant == "Significant")

#number of rows
#nrow(results_df_filter)
#nrow(Sig_df)
nrow(LolCDE_df)  # Number of rows in the LolCDE data frame
nrow(LolCDE_Sig_df)  # Number of rows in the significant LolCDE data frame

nrow(B01_df)  # Number of rows in the B01 data frame
nrow(B01_Sig_df)  # Number of rows in the significant B01 data frame

nrow(Ceftriaxone_df)  # Number of rows in the Ceftriaxone data frame
nrow(Ceftriaxone_Sig_df)  # Number of rows in the significant Ceftriaxone data frame

nrow(Chloramphenicol_df)  # Number of rows in the Chloramphenicol data frame
nrow(Chloramphenicol_Sig_df)  # Number of rows in the significant Chloramphenicol data frame

nrow(Ciprofloxacin_df)  # Number of rows in the Ciprofloxacin data frame
nrow(Ciprofloxacin_Sig_df)  # Number of rows in the significant Ciprofloxacin data frame

nrow(Clarithromycin_df)  # Number of rows in the Clarithromycin data frame
nrow(Clarithromycin_Sig_df)  # Number of rows in the significant Clarithromycin data frame

nrow(Colistin_df)  # Number of rows in the Colistin data frame
nrow(Colistin_Sig_df)  # Number of rows in the significant Colistin data frame

nrow(Doxycycline_df)  # Number of rows in the Doxycycline data frame
nrow(Doxycycline_Sig_df)  # Number of rows in the significant Doxycycline data frame

nrow(Globomycin_df)  # Number of rows in the Globomycin data frame
nrow(Globomycin_Sig_df)  # Number of rows in the significant Globomycin data frame

nrow(Levofloxacin_df)  # Number of rows in the Levofloxacin data frame
nrow(Levofloxacin_Sig_df)  # Number of rows in the significant Levofloxacin data frame

nrow(Mecillinam_df)  # Number of rows in the Mecillinam data frame
nrow(Mecillinam_Sig_df)  # Number of rows in the significant Mecillinam data frame

nrow(Meropenem_df)  # Number of rows in the Meropenem data frame
nrow(Meropenem_Sig_df)  # Number of rows in the significant Meropenem data frame

nrow(Nitrofurantoin_df)  # Number of rows in the Nitrofurantoin data frame
nrow(Nitrofurantoin_Sig_df)  # Number of rows in the significant Nitrofurantoin data frame

nrow(Norfloxacin_df)  # Number of rows in the Norfloxacin data frame
nrow(Norfloxacin_Sig_df)  # Number of rows in the significant Norfloxacin data frame

nrow(Nitroxolin_df)  # Number of rows in the Nitroxolin data frame
nrow(Nitroxolin_Sig_df)  # Number of rows in the significant Nitroxolin data frame

nrow(PolymyxineB_df)  # Number of rows in the PolymyxineB data frame
nrow(PolymyxineB_Sig_df)  # Number of rows in the significant PolymyxineB data frame

nrow(Trimethoprim_df)  # Number of rows in the Trimethoprim data frame
nrow(Trimethoprim_Sig_df)  # Number of rows in the significant Trimethoprim data frame

############################################################################################################################
#############################Making ML data matrix##################################################

# MOA == Compounds
data.frame(Position = 1:ncol(metadata), Column_Name = colnames(metadata))

MOA_compound <- metadata %>%
  select(8,12) %>%
  distinct() %>%
  filter(characteristics_ch1.2 != "moa:control")
