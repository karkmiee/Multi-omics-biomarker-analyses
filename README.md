# Multi-omics Approaches to Uncover Liquid-Based Cancer-Predicting Biomarkers in Lynch Syndrome

## Background  
Lynch Syndrome is a genetic cancer-predisposing syndrome caused by pathogenic mutations in DNA mismatch repair (path_MMR) genes. Due to the elevated cancer risk, novel screening methods are crucial to enhance cancer risk stratification. This project demonstrates how multi-omics integration can pinpoint cancer-predicting biomarkers in Lynch Syndrome. Specifically, we investigated blood-based circulating microRNAs and metabolites that predict cancer occurrence within a 5.8-year prospective surveillance period.

## Study methods and analytical workflow  

The study cohort included 116 Lynch Syndrome carriers who were healthy at the time of sampling, 17 of whom developed cancer during surveillance. The analyses included the following:

1. **Principal Coordinate Analysis (PCoA)**:  
To explore omics-level differences between pathogenic mismatch repair gene variant carriers  and healthy and future cancer groups. 

2. **Partial Least Squares (PLS) regression and sparse PLS (sPLS)**:  
To uncover associations between circulating microRNAs and metabolites.

3. **Weighted Correlation Network Analysis (WGCNA)**:  
To identify omics-level co-expression modules and their associations with cancer incidence or path_MMR variant.

4. **Lasso Cox Regression**:  
To identify cancer-predicting biomarkers and internally validate the predictive model using 5-fold cross-validation.




## Repository Contents  

This repository contains R scripts used in the analyses described above. Below is a brief description of each script:

- **`cmiR_filtering&normalization.R`**: Scripts for microRNA data cleaning, normalization, and handling batch effect.
- **`Boxcox_transformation.R`**: Scripts for boxcox-transformation of metabolomics data.
- **`PcoA_cmiR.R`**: Scripts for performing Principal Coordinate Analysis.  
- **`PcoA_metabolomics.R`**: Scripts for performing Principal Coordinate Analysis.  
- **`WGCNA_cmiR.R`**: Scripts for Weighted Correlation Network Analysis (microRNAs).  
- **`WGCNA_metabolomics.R`**: Scripts for Weighted Correlation Network Analysis (metabolomics).
- **`MixOmics_miR_metab.Rmd`**:   Scripts for Partial Least Squares (PLS) regression and sparse PLS (sPLS).
- **`multiomics_split 1.py`**: Scripts for randomized datasplit.
- **`Lasso_Cox_regression.R`**: Scripts for Lasso Cox regression and model validation.  


For questions or collaboration inquiries, please contact Minta Kärkkäinen (minta.e.m.karkkainen@jyu.fi) or corresponding author, Research director, PhD. Tiina Jokela (tiina.a.jokela@jyu.fi). Note: Due to ethical standards, not all phenotypic data required for analyses can be publicly shared. 
