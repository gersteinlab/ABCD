

#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(manta)
library(car)

#********
# BEGIN *
#********

# 1. input folder:
inFolder="data4gwas/data/"


#---------------
# v1.v2.wrCovs |
#---------------

pheno <- read.delim(paste0(inFolder, "v1.v2/phenotypes.tsv"), h=T, stringsAsFactors = F)
covs <- read.delim(paste0(inFolder, "v1.v2/covariates.tsv"), h=T, stringsAsFactors = F)

pca.pheno <- prcomp(pheno[, 2:ncol(pheno)], scale = TRUE)

pheno.PCs <- as.data.frame(pca.pheno$x[, 1:5])
colnames(pheno.PCs) <- paste0("pheno_", colnames(pheno.PCs))
pheno.PCs$ID <- covs$ID
pheno.PCs <- pheno.PCs[, c(ncol(pheno.PCs), 1:(ncol(pheno.PCs)-1))]
stopifnot(identical(pheno.PCs$ID, covs$ID))

write.table(pheno.PCs, "v1.v2.5PCs/PC.phenotypes.tsv",
            row.names = F, col.names = T, sep="\t", quote = F) 

