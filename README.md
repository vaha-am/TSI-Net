# TSI-Net
Multi-modal network inference framework to construct Transcriptome-Small Molecule Interaction Network (TSI-Net)

This repository contains code implementations for the paper "Construction of Multi-Modal Transcriptome-Small Molecule Interaction Networks from High-Throughput Measurements to Study Human Complex Traits" (doi: https://doi.org/10.1101/2025.01.22.634403)

# Short Description
The project is divided into three parts: 1. A protocol for processing liquid chromatography/mass spectrometry (LC/MS) data, 2. A semi-supervised method for constructing multi-modal networks, and 3. An analysis of the multi-modal TSI-Nets with respect to human complex traits. Code implementations for parts 1 and 2 are available in this repository.

# LC/MS data processing protocol
The LC/MS data processing protocol consists of three steps: 1. Batch effect correction using a random forest model, 2. principal component analysis (PCA), and 3. covariate adjustment with a stepwise regression model

## Batch effect correction
Batch effect correction was implemented from https://colab.research.google.com/drive/1CLB6WNPN8JJAwezuKQpkbNztR4o6qqFK?usp=sharing based on the paper from Stancliffe et al[1].
To correct for batch effects, a random forest-based approach was applied to the raw peak areas. After log2 transformation of peak areas, batch number and QC sample run order were used to model each small molecule’s (SM) intensity deviation in each QC sample from the mean intensity of that molecule across all QC samples, with one model per small molecule. The model was trained on QC samples, and applied on research samples to predict the deviation of SM peak areas in research samples. The predicted deviation was subtracted from the peak area for each SM, yielding a semi-processed data

## PCA
After batch effect correction, principal components (PCs) were calculated based on the log₂-transformed, semi-processed LC/MS data to obtain the top 10 PCs. These PCs are used in the covariate adjustment step. The PCA code implementation can be found in the [PCA directory](./PCA/).

## Covariate adjustment
We use a stepwise regression model to adjust the input matrix (in this case, the semi-processed log-transformed LC/MS data) for a predefined set of variables. The stepwise regression model consists of:
  1. Base variables that are always included in the model, regardless of their relationship with the outcome variable.
  2. Additional variables that are included only if they meet a significance threshold for the outcome variable.

A separate stepwise regression model is applied to each feature in the input matrix (i.e., each semi-processed LC/MS measurement).
The same methodology is also applied to the transformed RNA-seq data. For more details, please refer to the paper.
The stepwise regression code implementation can be found in the
