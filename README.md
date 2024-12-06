# Genetic Covariance and Variance-Covariance Matrix Estimation for Heterogeneous Stock Rats

This repository provides resources and scripts to estimate genetic covariance and variance-covariance matrices required for the **Genomic SEM R package**. The estimation is performed using **[MPH (Multi-component Penalized Heritability)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11093526/)**, customized specifically for research involving **Heterogeneous Stock (HS) rats**.

### Scope of the Code
This repository is designed for cases where **traits are measured in different cohorts**, rather than in a single cohort with multiple traits. Specifically:
- **Traits A, B, and C are measured in cohorts 1, 2, and 3, respectively.**  
- The goal is to estimate the genetic covariance matrix between the three traits, accounting for the cohort structure and ensuring accurate variance-covariance modeling.
