# 2024-06-16 Generation_Lab_proj, Yujie Chen
## DNA methylation data (run by EPICv1/850k chip) from public datasets and internal dataset
## Goal:  Ensure consistency across all data

#Step 1: Load necessary libraries
library(minfi) 
library(sva)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggplot2)
library(dplyr)

# save sessionInfo
writeLines(capture.output(sessionInfo()), "Generation_lab_proj_sessionInfo.txt")


#Step 2: Load Data
# Load IDAT files 
public_baseDir <- "~/Downloads/public_datasets_with_metadata"
internal_baseDir <- "~/Downloads/internal_dataset_with_metadata"

#"GSM4444086_201496850037_R03C01_Grn.idat" is removed, no paired Red.idat file
public_rgSet <- read.metharray.exp(base = public_baseDir, force=TRUE) 
internal_rgSet <- read.metharray.exp(base = internal_baseDir)

# Load metadata, reformat metadata, consider each internal user data and each publication as a batch
public_metadata <- read.csv("~/Downloads/public_dataset_metadata.csv") %>% 
  select(idat_match, sex, Age, GSE_ID) %>%
  rename(sample_id = idat_match, sex = sex, age = Age, batch = GSE_ID)

internal_metadata <- read.csv("~/Downloads/internal_dataset_metadata.csv") %>%
  rename(sample_id = Sample_ID, sex = Sex, age = Age, batch = Batch) %>%
  mutate(sex = ifelse(sex == "M", "male", "female"))

# Combine metadata
combined_metadata <- rbind(public_metadata, internal_metadata)


#Step 3: Quatlity control
#Combine rgSets
combined_rgSet <- combineArrays(public_rgSet, internal_rgSet, outType = "IlluminaHumanMethylationEPIC")

#Step 4: Quality Control
# Perform quality control
qcReport(combined_rgSet, pdf = "./QC_Report.pdf")
qcReport(public_rgSet, pdf = "./QC_Report_public.pdf")
qcReport(internal_rgSet, pdf = "./QC_Report_internal.pdf")

# Detect poor quality samples
options(matrixStats.useNames.NA = "deprecated")
detP <- detectionP(combined_rgSet,type="m+u")

failedSamples <- colMeans(detP > 0.01) > 0.05
combined_rgSet <- combined_rgSet[, !failedSamples]

# Identify and remove poor quality probes
failedProbes <- rowMeans(detP > 0.01) > 0.01
combined_rgSet <- combined_rgSet[!failedProbes, ]


# Step 5: Normalization
# Use Noob for Correction (Tested SWAN, Noob works better than SWAN based on normalized beta density distribution)
noobSet <- preprocessNoob(combined_rgSet) 

# Check density distribution across samples after normalization
phenoData <- pData(noobSet)
densityPlot(noobSet)
ggsave("Beta_Value_density_plot_after_NoobNorm.pdf")


# Step 6: Batch Effect Correction
# Extract beta values
#beta_values <- getBeta(swanSet)
beta_values <- getBeta(noobSet)

# Prepare batch and covariate information
combined_metadata.match <- combined_metadata[combined_metadata$sample_id %in% colnames(beta_values), ]
batch <- combined_metadata.match$batch
covariates <- model.matrix(~ sex + age, data = combined_metadata.match)

# Perform batch correction
beta_corrected <- ComBat(dat = beta_values, batch = batch, mod = covariates)


#Step 7: Further Quality Checks
# Visualize data before and after batch correction using PCA
plotPCA <- function(beta_values, group_var, title, file_name) {
  pca <- prcomp(t(beta_values))
  pca_data <- data.frame(SampleID = rownames(pca$x), 
                         PC1 = pca$x[,1], 
                         PC2 = pca$x[,2], 
                         Group = group_var)
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) + 
    geom_point() + 
    ggtitle(title)
  ggsave(file_name, plot = p)
}

# Extract metadata for matched samples
combined_metadata.match <- combined_metadata[combined_metadata$sample_id %in% colnames(beta_values), ]
batch <- combined_metadata.match$batch
sex <- combined_metadata.match$sex
age <- combined_metadata.match$age
sample <- combined_metadata.match$sample_id

# PCA plots before and after batch correction
plotPCA(beta_values, batch, "PCA Before Batch Correction (Noob) - Batch", "PCA_Before_Batch_Correction_Noob_Batch.pdf")
plotPCA(beta_corrected, batch, "PCA After Batch Correction (Noob) - Batch", "PCA_After_Batch_Correction_Noob_Batch.pdf")

plotPCA(beta_values, sex, "PCA Before Batch Correction (Noob) - Sex", "PCA_Before_Batch_Correction_Noob_Sex.pdf")
plotPCA(beta_corrected, sex, "PCA After Batch Correction (Noob) - Sex", "PCA_After_Batch_Correction_Noob_Sex.pdf")

plotPCA(beta_values, age, "PCA Before Batch Correction (Noob) - Age", "PCA_Before_Batch_Correction_Noob_Age.pdf")
plotPCA(beta_corrected, age, "PCA After Batch Correction (Noob) - Age", "PCA_After_Batch_Correction_Noob_Age.pdf")

plotPCA(beta_values, sample, "PCA Before Batch Correction (Noob) - Sample", "PCA_Before_Batch_Correction_Noob_Sample.pdf")
plotPCA(beta_corrected, sample, "PCA After Batch Correction (Noob) - Sample", "PCA_After_Batch_Correction_Noob_Sample.pdf")


# Step 8: Save Corrected Data
## Split the corrected data back into public and internal datasets
public_corrected <- beta_corrected[, 1:ncol(public_rgSet)]
internal_corrected <- beta_corrected[, (ncol(public_rgSet) + 1):ncol(beta_corrected)]

# Save to CSV
write.csv(public_corrected, "public_beta_values_corrected_Noob.csv")
write.csv(internal_corrected, "internal_beta_values_corrected_Noob.csv")

# Save to CSV
write.csv(public_corrected, "public_beta_values_corrected_Noob.csv")
write.csv(internal_corrected, "internal_beta_values_corrected_Noob.csv")

# Save each sample's beta values into separate CSV files
saveIndividualCSVs <- function(beta_matrix, output_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  for (sample_id in colnames(beta_matrix)) {
    sample_beta_values <- beta_matrix[, sample_id, drop = FALSE]
    write.csv(sample_beta_values, file.path(output_dir, paste0(sample_id, "_beta_values.csv")), row.names = TRUE)
  }
}

# Save individual CSVs for public and internal corrected data
saveIndividualCSVs(public_corrected, "public_corrected_samples")
saveIndividualCSVs(internal_corrected, "internal_corrected_samples")
