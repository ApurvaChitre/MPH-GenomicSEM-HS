# Genetic Covariance and Variance-Covariance Matrix Estimation for Heterogeneous Stock Rats

This repository provides resources and scripts to estimate genetic covariance and variance-covariance matrices required for the **Genomic SEM R package**. The estimation is performed using **[MPH (MINQUE for Partitioning Heritability)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11093526/)**, customized specifically for research involving **Heterogeneous Stock (HS) rats**.

### Scope of the Code
This repository is designed for cases where **traits are measured in different cohorts**, rather than in a single cohort with multiple traits. In this context, cohort does not refer to a phenotyping batch. Specifically:
- **Cohorts represent distinct Heterogeneous Stock rats (HS) projects arising from collaborations with the Palmer lab, with each project focusing on a unique set of traits.**
- For example, **Traits A, B, and C are measured in cohorts 1, 2, and 3, respectively.**

The goal is to estimate the genetic covariance matrix between the three traits, accounting for the cohort structure and ensuring accurate variance-covariance modeling.

### Goal of the Analysis

The objective is to compute two key matrices for use in the **Genomic SEM R package**:

1. **Matrix S (Genetic Covariance Matrix):**
   - A symmetric $n \times n$ matrix for $n$ traits.
   - **Diagonal**: Genetic variance for each trait.
   - **Off-diagonal**: Genetic covariances between traits.

2. **Matrix V (Variance-Covariance Matrix of Parameter Estimates):**
   - A symmetric $\frac{n(n+1)}{2} \times \frac{n(n+1)}{2}$ matrix.
   - **Diagonal**: Squared standard errors of genetic variances.
   - **Off-diagonal**: Covariances of standard errors for genetic variances and covariances.

This method overcomes limitations of **LDSC (Linkage Disequilibrium Score Regression)**, which provides out-of-bound estimates for HS rats.

### Correspondence Highlights with the author of MPH, Jicai Jiang

The function for the problem has not been hard-coded in MPH, but MPH does have a flexible (though requiring more data manipulation) approach to doing so. You can stack the phenotypes across cohorts as a single trait. By customizing genomic relationship matrices, the estimation problem for multiple cohorts is reduced to a single-variate multi-component model.


Suppose **G** is a genomic relationship matrix across three cohorts (1-3):

$$
\mathbf{G} =
\begin{pmatrix}
G_{11} & G_{12} & G_{13} \\
G_{21} & G_{22} & G_{23} \\
G_{31} & G_{32} & G_{33}
\end{pmatrix}
$$

Suppose $\mathbf{y}_1, \mathbf{y}_2, \mathbf{y}_3$ are phenotypes for cohorts 1-3. 

$$
\begin{pmatrix} 
\mathbf{y}_1 \\ 
\mathbf{y}_2 \\ 
\mathbf{y}_3 
\end{pmatrix} 
$$

The **variance-covariance matrix** for phenotypes measured across three cohorts is as follows:

![Variance-Covariance Matrix](formula.jpg)

- **Matrix on the Left**: Represents genetic variance ($\sigma^2_g$) and genetic covariance between traits across cohorts.
  - **Diagonal elements**: Genetic variance for each cohort.
  - **Off-diagonal elements**: Genetic covariance between cohorts.

- **Matrix on the Right**: Represents residual variance ($\sigma^2_e$).
  - **Diagonal elements**: Residual variance for each cohort.
  - **Off-diagonal elements**: Zero, indicating no residual covariance between cohorts.

This can be expanded as follows:
![Expanded Variance-Covariance Matrix](expanded_formula.jpg)


### Overview of steps for the analysis (3 Cohorts)

The whole analysis (for 3 cohorts) includes a few major steps:

1. **Compute G Across Cohorts**  
   Compute the genomic relationship matrix (G) across all cohorts.

2. **Create Custom Matrices**  
   - Generate the six genomic matrices and the 3  residual matrices (shown above).  
   - Only 3 residual matrices are required because MPH automatically includes a diagonal matrix as the sixth one.  
   - Use the provided R scripts for manipulating GRMs [GRM Scripts](https://jiang18.github.io/mph/scripts/#grm-inputoutput)
   - This step should be quick if the total number of rats across cohorts is around 10k.

3. **List Custom Relationship Matrices**  
   - Include all 9 custom relationship matrices in the GRM list file for MPH.

4. **Run MPH**  
   - Execute MPH, which will output the matrices **S** and **V** in the `.mq.vc.csv` file.

### Detailed Steps and Commands

Detailed instructions to process and prepare files for MPH.


#### 1. Prepare Phenotype and Covariate Files

This step combines phenotype data from multiple files and prepares cohort-specific covariates.

**Requirements:**  
- **Phenotype Files Directory:** Directory containing individual phenotype files.  
- **Cohort Mapping File:** Maps cohort numbers to phenotype files.  
  - 📄 **Available here:** `/tscc/projects/ps-palmer/apurva/locomotor/mph/individual_projects_seq/pheno_cohort_project_dict.csv`  

**Modify File Paths:**  
Update these paths in the script:

```r
pheno_dir <- "/path/to/phenotype/files"
pheno_cohort_file <- "/path/to/cohort_mapping.csv"
output_pheno_file <- "/path/to/output/combined_phenotype.csv"
output_cov_file <- "/path/to/output/indicator_covariates.csv"
```

**Covariate File Structure:**  
The covariate file acts as an **incidence matrix**, where each cohort is represented as a column, and IDs indicate whether an individual belongs to a specific cohort.

Example structure:

```csv
ID,cohort1,2,3,4,5,6
A,1,0,0,0,0,0
B,1,0,0,0,0,0
C,0,1,0,0,0,0
```

- **ID:** Unique identifier for each sample.  
- **cohortX:** Binary indicator (1 = belongs to cohort, 0 = does not belong).  

**Run the Script:**  
Save the script as `stack_pheno_create_cov.R` and execute it:

```bash
Rscript stack_pheno_create_cov.R
```

**Output Files:**  
- `combined_phenotype.csv`: Combined phenotype data.  
- `indicator_covariates.csv`: Cohort indicator covariates.


#### 2. Make GRM (Genomic Relationship Matrix)

**📝 Note:** Ensure that the **RFIDs/Samples** overlap between the **genotype** and **phenotype** data. MPH will error out even if the genotype file contains additional samples not present in the phenotype data, even if there is overlap. This approach is intended for autosomes only. 

### **Step 1: Extract IDs from the Covariate File**
```bash
cut -d',' -f1 /tscc/projects/ps-palmer/apurva/locomotor/mph/individual_projects_no_dups/indicator_covariates.csv | tail -n +2 > /tscc/projects/ps-palmer/apurva/locomotor/mph/genotypes/subset_ids.txt
```

### **Step 2: Subset Genotype Files**
```bash
plink --bfile /tscc/projects/ps-palmer/apurva/locomotor/individual_projects_no_dups/genotypes/genotypes \
      --keep-fam /tscc/projects/ps-palmer/apurva/locomotor/mph/genotypes/subset_ids.txt \
      --chr 1-20 \
      --make-bed \
      --out /tscc/projects/ps-palmer/apurva/locomotor/mph/genotypes/subset_genotypes
```


Generate the GRM using the following commands:

```bash
perl -e 'print "SNP\n"; while(<>){@c=split /\s+/; print "$c[1]\n"}' < /tscc/projects/ps-palmer/apurva/locomotor/mph/genotypes/subset_genotypes.bim > snp_info.csv
```

Run the mph command to create the GRM:
**Note:** It is best to submit this as a **batch job** on the cluster since it is computationally intensive. See the example batch script: `grm_job.sh`.
```bash
mph --make_grm --binary_genotype output_gen_hard_calls_subset_mph.qc --snp_info snp_info.csv --num_threads 14 --out gg
```

#### 3. Generate All GRMs
To generate all GRMs, use the following commands:
```bash
mkdir grm
mkdir reml

Rscript make_cohort_grms.R
```

#### 4. Run REML Analysis

Perform the REML (Restricted Maximum Likelihood) analysis using mph:
```bash
mph --reml --grm_list gg.list.txt --phenotype combined_phenotype.csv --trait trait_value --covariate_file indicator_covariates.csv --covariate_names all --output reml/all
```

### Output Directory

The results of the analysis are in the `output/` directory, which contains the following key files:

- **`all.mq.*`**: The primary output files from the analysis, including the estimates and variances.
- **`cohort_VCs.xlsx`**: A spreadsheet highlighting the estimates and (co)variances of the estimates of interest.

For reproducing the results, refer to the accompanying `README` file and the `make_cohort_grms.R` script.

### File Formats

For details about the required file formats, please refer to:

- The **[MPH Manual](https://jiang18.github.io/mph/options/)**, which provides comprehensive documentation on all input and output file structures.
- The Palmer Lab TSCC directory at `/tscc/projects/ps-palmer/s3/data/genomic_sem`, which contains example files and datasets for genomic SEM analyses.

These resources will provide all necessary guidance on preparing and using the required files for analysis.


### Acknowledgments

Special thanks to **[Jicai Jiang](https://cals.ncsu.edu/animal-science/people/jicai-jiang/)** for his invaluable guidance on using the **MPH framework** and providing the code constructing custom matrices.

For more details, visit the **[MPH manual](https://jiang18.github.io/mph/)**.









