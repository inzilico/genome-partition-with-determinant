# Genome Partitioning Using Determinant of LD Matrix 

## Description

The repo contains the scripts applied in the [article](DOI: 10.1093/nargab/lqaf024) 

## Use cases

### Genome-wide association studies (GWAS)

Following multimarker approaches in GWAS, we proposed the workflow: 

1. Get an LD-matrix with $r^2$ values. 

2. Transform LD-matrix into block-diagonal one. 

3. Restore haplotypes within each block

4. Test each block and identify those that are significantly associated with the phenotype.

### Other potential usages

* tagSNP selection
* dimensionality reduction for machine learning approaches 


