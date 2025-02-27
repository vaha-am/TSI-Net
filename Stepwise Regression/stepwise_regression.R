library(tidyverse)
library(MASS)
library(dplyr)
library(here)
library(stringr)
library(fastDummies)
library(hexbin)
library(ggplot2)
library(reshape2)

# We will input everything here. This includes metabolite/lipid data, covariates, and phenotype
# Currently, covariates are name hypothetically: {varA}, {varB}, {varC}, ..., {varN}
# In the TSI-Net paper, the covariates used for adjustment with stepwise regression include Age, Age^(2), sex, clinical field centers, medication data, smoking status, and metabolom/lipidome principle components
# In this code, we assume we have an input for 1. metabolome peak areas, 2. principle components, 3. medication, 4. age/sex/field centers, and 5. smoking
metabolite_peak_df <- read.csv("/path/to/raw/metabolomic/data.csv")
principle_components_df <- read.csv("/path/to/metabolomic/pc.csv")
medication_data_df <- read.csv("/path/to/medication/csv")
age_sex_field_centers_df <- read.csv("/path/to/age_sex_fc.csv")
smoking_df <- read.csv("/path/to/smoking.csv")

# Selection of relevant columns

# All metabolomic data
metabolite_peak_df <- metabolite_peak_df %>% dplyr::select(everything())

# Log transformation for all columns of the metabolomic data, except metadata columns (e.g., subject IDs)
metadata_columns <- c("subjectID", "{meta_data2}", "...", "{meta_dataN}")
cols_to_transform <- setdiff(colnames(metabolite_peak_df), metadata_columns)
metabolite_peak_df[cols_to_transform] <- lapply(metabolite_peak_df[cols_to_transform], log)

# Covariate selection
principle_components_df <- principle_components_df %>% dplyr::select({meta_data1}, {varA}:{varB})
medication_data_df <- medication_data_df %>% dplyr::select({meta_data1}, {varC}, {varD}, {varE}, {varF}, {varG})
age_sex_field_centers_df <- age_sex_field_centers_df %>% dplyr::select({meta_data1}, {varH}, {varI}, {varJ}, {varK}, {varL})
smoking_df <- smoking_df %>% dplyr::select({meta_data1}, {varM})

# Combining dataframes by subjects, removing missing values
dfs_list <- list(metabolite_peak_df, principle_components_df, age_sex_field_centers_df, smoking_df, medication_data_df)
combined_df <- dfs_list %>% purrr::reduce(full_join, by = "subjectID")
combined_df[combined_df == ""] <- NA
combined_df <- combined_df %>% drop_na()

# Creating a column for age^(2)
combined_df$age_2 <- combined_df${varK} ^ 2

rownames(combined_df) <- combined_df${meta_data1}
subject_ids <- rownames(combined_df)

# Output the combined_df as a CSV file for inspection
write.csv(combined_df, "/output/path/combined_df.csv", row.names = FALSE)

# Selecting metabolite columns
metabolite_colnames <- colnames(metabolite_peak_df)[-1] 

# Determining the positions of metabolite columns in the combined dataframe
metabolite_positions <- which(colnames(combined_df) %in% metabolite_colnames)

# Printing out the positions of metabolite columns in the combined df
print(metabolite_positions)

# Printing the number of metabolite columns
print(paste("Number of columns: ", ncol(combined_df)))

# Defining base covariates for stepwise regression
base_covariates_dfs <- list(age_sex_field_centers_df)
base_covariates <- purrr::map(base_covariates_dfs, function(df) setdiff(colnames(df), "subjectID")) %>% unlist()
# Adding "age_2" to the base covariates
base_covariates <- c(base_covariates, "age_2")
base_covariates_positions <- which(colnames(combined_df) %in% base_covariates)

# Number of base covariates
print(paste("Number of base covariates: ", length(base_covariates_positions)))
print(base_covariates)

# Defining additional covariates for stepwise regression
additional_covariates_dfs <- list(principle_components_df, smoking_df, medication_data_df)
additional_covariates <- purrr::map(additional_covariates_dfs, function(df) setdiff(colnames(df), "subjectID")) %>% unlist()
additional_covariates_positions <- which(colnames(combined_df) %in% additional_covariates)

# Number of additional covariates
print(paste("Number of additional covariates: ", length(additional_covariates_positions)))
print(additional_covariates)

# Total number of covariates (base + additional)
print(paste("Total number of covariates: ", length(base_covariates_positions) + length(additional_covariates_positions)))

# Initializing residuals storage vectors based on metabolite peaks
residuals_normality_shapiro_p_values <- vector(mode = 'list', length = length(metabolite_positions))
residuals_normality_shapiro_w_values <- vector(mode = 'list', length = length(metabolite_positions))

# Initializing residuals_df for storing outputs
residuals_df = data.frame(matrix(ncol = length(metabolite_positions), nrow = nrow(combined_df)))
colnames(residuals_df) <- colnames(combined_df)[metabolite_positions]

# Defining "subject" column in residuals_df
residuals_df${meta_data1} <- combined_df${meta_data1}

# Setting up residual adjustment using stepwise regression

# summary directory will store summary statistics of the stepwise regression for each metabolite as a txt file - 1 txt file per metabolite 
summary_directory <- "/path/to/summary/files"
dir.create(summary_directory, recursive = TRUE, showWarnings = FALSE)

# Stepwise regression loop
for (metabolite_position in metabolite_positions) {
  # Creating a copy of combined_df
  temp_df <- combined_df
  temp_df$metabo <- as.numeric(temp_df[, metabolite_position])  # Changed from 'gene' to 'metabo'
  
  # Performing stepwise regression
  fixed_effects <- as.formula(paste("metabo", paste(c(base_covariates, additional_covariates), collapse = "+"), sep = "~"))
  model <- lm(fixed_effects, data = temp_df)
  model_step <- step(model, scope = list(upper = as.formula(paste("~", paste(c(base_covariates, additional_covariates), collapse = "+"))), 
                                         lower = as.formula(paste("~", paste(base_covariates, collapse = "+")))), trace = FALSE)
  temp_residuals <- resid(model_step)
  
  # Saving the summary of the model to a text file named after the metabolite
  metabolite_name <- colnames(combined_df)[metabolite_position]
  summary_file_path <- file.path(summary_directory, paste0(metabolite_name, ".txt"))
  sink(summary_file_path)
  print(summary(model_step))
  sink()  # Stop redirecting output to file
  
  # Assigning residuals to the corresponding metabolite column in residuals_df
  residuals_df[, metabolite_name] <- temp_residuals
}

# Reordering columns so that subject is the first column
residuals_df <- dplyr::select(residuals_df, "subjectID", dplyr::everything())

# Residuals output
write.csv(residuals_df, "/path/to/processed/metabolomic_file.csv", row.names = FALSE)
