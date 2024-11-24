## GWAS 

### Univariate GWAS
For univariate GWAS we employed [plink2](https://www.cog-genomics.org/plink/2.0/).

#### Binary trait 
```
plink2 \
  --threads 6 \
  --out out/gwas_sumstats \
  --keep indivs2keep.txt \
  --pfile imputed/pgen.files/genotype \
  --pheno phenotypes.binary.tsv \
  --chr 1-22,X \
  --covar covariates.tsv \
  --covar-variance-standardize \
  --glm firth-fallback hide-covar omit-ref no-x-sex
```
#### Continuous trait
```plink2 \
  --threads 6 \
  --out out/gwas_sumstats \
  --keep indivs2keep.txt \
  --pfile imputed/pgen.files/genotype \
  --pheno phenotypes.quantitative.tsv \
  --pheno-name cbcl_externalizing,cbcl_internalizing,liability_Xception_without_cbcl,liability_Xception_with_CBCL,liability_Xception_without_cbcl_v2,liability_Xception_with_CBCL_v2,XGB_without_cbcl,XGB_with_cbcl,XGB_without_cbcl_v2,XGB_with_cbcl_v2 \
  --chr 1-22,X \
  --covar covariates.tsv \
  --covar-variance-standardize \
  --glm hide-covar omit-ref no-x-sex
```

### Multivariate GWAS
For multivariate GWAS we employed [mvgwas-nf](https://github.com/dgarrimar/mvgwas-nf).
```
nextflow run mvgwas.nf --l 1000 --geno all.chr.vcf.gz --pheno phenotypes.tsv --cov covariates.tsv --out mvgwas.tsv -resume -with-singularity -with-trace -bg -with-mpi
```

For all other code requests please email mark@gersteinlab.org
