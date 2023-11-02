#!/usr/bin/env Rscript

#************
# LIBRARIES *
#************

library(data.table)



#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  
  make_option( c ( "-f", "--features" ), default=NULL,
               help = "The features matrix used to construct the phenotypes."),
  
  make_option( c( "-o", "--output_folder" ), default = NULL,
               help = "Output folder where to save phenotype.tsv & covariates.tsv [default = %default]." )
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nPrepares phenotypes and covariates matrices for gwas analysis."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#********
# BEGIN *
#********

# opt <- list()
# opt$features <- "/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/bb926/manta/data/JL.features.matrices/formatted.matrices/v1.tsv"


# 1. read df w/ wearables features matrix by JL 
m <- read.delim(opt$features, stringsAsFactors = F)

# 2. read metadata info
covs <- read.delim("/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/bb926/metadata/metadata.tsv", stringsAsFactors = F)

# 3. common indivs.
sel.indvs <- intersect(m$src_subject_id, covs$src_subject_id)
m <- m[m$src_subject_id %in% sel.indvs, ]
covs <- covs[covs$src_subject_id %in% sel.indvs, ]


# 4. reorder rownames and make sure the order is the same in m and covs
m$src_subject_id <- factor(m$src_subject_id, levels = sel.indvs)
m <- m[order(m$src_subject_id), ]
covs$src_subject_id <- factor(covs$src_subject_id, levels = sel.indvs)
covs <- covs[order(covs$src_subject_id), ]
stopifnot(identical(m$src_subject_id, covs$src_subject_id))


# 5. add col ID (RUID + src_subject_id) to features matrix
m$ID <- paste(covs$RUID, covs$src_subject_id, sep="_") 
m$src_subject_id <- NULL
m <- m[, c(ncol(m), 1:(ncol(m)-1))]


# 6. add col ID + reduce number of covariates
covs$ID <- paste(covs$RUID, covs$src_subject_id, sep="_") 
covs$RUID <- NULL
covs$src_subject_id <- NULL
sel.covs <- c("acs_raked_propensity_score", 
              "sex", 
              "interview_age",
              "BATCH",
              "PC1",
              "PC2",
              "PC3",
              "PC4",
              "PC5",
              "cbcl_externalizing", 
              "cbcl_internalizing", 
              "dsm_externalizing", 
              "dsm_internalizing", 
              "famhis_bip", 
              "famhis_schiz", 
              "famhis_antisocial", 
              "famhis_nerves",
              "famhis_treatment", 
              "famhis_hospital", 
              "famhis_suicide", 
              "parent_div_cat", 
              "parent_grade", 
              "family_income", 
              "adopted")

covs <- covs[, colnames(covs) %in% c("ID", sel.covs)]              
covs <- covs[, c(25, 1:24)]              
              
 

# 10. save phenotypes and covariates matrices
write.table(m, file = paste0(opt$output_folder, "/phenotypes.tsv"),
            row.names = F, col.names = T, sep="\t", quote = F)

write.table(covs, file = paste0(opt$output_folder, "/covariates.tsv"),
            row.names = F, col.names = T, sep="\t", quote = F)

              
