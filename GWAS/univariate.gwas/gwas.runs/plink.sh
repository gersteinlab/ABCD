# ******
# GWAS |
# ******

# univariate GWAS for binary trait
plink2 \
  --threads 6 \
  --out out/gwas_sumstats \
  --keep case.ctrl.gwas/adhd/binary.quantitative.gwas/gwas.runs/1191/indivs2keep.txt \
  --pfile imputed/pgen.files/genotype \
  --pheno case.ctrl.gwas/adhd/binary.quantitative.gwas/gwas.runs/1191/phenotypes.binary.tsv \
  --chr 1-22,X \
  --covar covariates.allCovs.tsv \
  --covar-variance-standardize \
  --glm firth-fallback hide-covar omit-ref no-x-sex

# univariate GWAS for quantitative (continuous) traits
plink2 \
  --threads 6 \
  --out out/gwas_sumstats \
  --keep case.ctrl.gwas/adhd/binary.quantitative.gwas/gwas.runs/1191/indivs2keep.txt \
  --pfile imputed/pgen.files/genotype \
  --pheno case.ctrl.gwas/adhd/binary.quantitative.gwas/gwas.runs/1191/phenotypes.quantitative.tsv \
  --pheno-name cbcl_externalizing,cbcl_internalizing,liability_Xception_without_cbcl,liability_Xception_with_CBCL,liability_Xception_without_cbcl_v2,liability_Xception_with_CBCL_v2,XGB_without_cbcl,XGB_with_cbcl,XGB_without_cbcl_v2,XGB_with_cbcl_v2 \
  --chr 1-22,X \
  --covar covariates.allCovs.tsv \
  --covar-variance-standardize \
  --glm hide-covar omit-ref no-x-sex

