# 1. obtain clinical covariates
Rscript clinical.covariates.R

# 2. integrate clinical covariates with other types of covariates
# (demographic, genotyping batch, sex, etc)
Rscript all.covariates.R
