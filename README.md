# Genome Partitioning Using Determinant of LD Matrix 

## Description

The repo contains the scripts applied in the article (DOI: 10.1093/nargab/lqaf024) 

## Use cases

### Genome-wide association studies (GWAS) based on testing of haplotypes

Following multimarker approaches in GWAS, we proposed the following workflow. Assume the original files are given as Plink 1.9 binary format. 

1. Get an LD-matrix with $r^2$ values. 

```bash
plink --bfile path/to/prefix --r2 square --out path/to/prefix 

```

2. Transform LD-matrix into block-diagonal one. 

3. Restore haplotypes within each block

4. Test each block and identify those that are significantly associated with the phenotype.

### Other potential usages

* tagSNP selection
* dimensionality reduction for machine learning approaches 


