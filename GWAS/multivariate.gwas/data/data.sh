# 1. prepare phenotype and covariate matrix for the 14 heart rate static features
Rscript get.pheno.covs.matrices.R -f formatted.matrices/v1.tsv -o ./v1

# 2. prepare phenotype and covariate matrix for the 258 wearable-derived static features
Rscript get.pheno.covs.matrices.R -f formatted.matrices/v1.v2.tsv -o ./v1.v2

# 3. compute top 5 Principal Components based on the 258 wearable-derived static features
# these 5 PCs constitute the phenotype matrix
# Note that for this GWAS we used the covariate file generated in point 2.
Rscript v1.v2.get.5PCs.R

# 4. cluster the 258 static features into 7 clusters
Rscript v1.v2.get.7clusters.R

# 5. prepare the phenotype matrix for each of the seven clusters
# Note that for this GWAS we used the covariate file generated in point 2. 
Rscript v1.v2.get.clusters.phenotypes.R

