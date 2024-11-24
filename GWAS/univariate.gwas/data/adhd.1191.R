#************
# LIBRARIES *
#************

library(data.table)


#********
# BEGIN *
#********

outfolder="ABCD/bb926/data4gwas/data/adhd.binary.quantitative/1191/"


#----------------------------
# read all metadata sources |
#----------------------------


# 1. metadata1 (ancestry file)
# covariates: sex, (ethnicity), propensity score
# we're not using this ethnicity, but just in case we keep this info in the metadata
# we, instead, use the first 5 PCs computed on the genotype data (see metadata5)

metadata1 <- fread("ABCD/bb926/metadata/acspsw03.txt", 
                   data.table = F, header = T, stringsAsFactors = F)
metadata1 <- metadata1[2:nrow(metadata1), ]
m1.cols <- c("src_subject_id", 
             "sex", 
             "race_ethnicity",
             "genetic_af_african", 
             "genetic_af_european", 
             "genetic_af_east_asian", 
             "genetic_af_american",
             "acs_raked_propensity_score")
metadata1 <- metadata1[metadata1$eventname == "baseline_year_1_arm_1", m1.cols]


# 2. metadata2 
# covariates: interview_age, sex, cohorts groups
metadata2 <- read.csv("ABCD/metadata/ABCD_metadata_03152022.csv", h=T)
m2.cols <- c("src_subject_id", 
             "interview_age",
             "sex",
             "adhd_group",
             "anxiety_disorder",
             "bipolar_group",
             "ocd_group",
             "severe_disorder",
             "panic_group",
             "pstd_group",
             "depression_group",
             "eating_group",
             "sleep_group",
             "nonclinical_controls",
             "nonsevere_controls",
             "internalizing_comparison_controls",
             "no_externalizing")
metadata2 <- metadata2[, m2.cols]


# 3. metadata3
# covariates: externalizing / internalizing gradients, family history
metadata3 <- fread("ABCD/bb926/metadata/wr.folder/clinical.covariates.tsv", data.table = F, header = T, stringsAsFactors = F)
m3.cols <- c("src_subject_id",
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
metadata3 <- metadata3[, m3.cols]


# 4. metadata 4 (sequencing info)
# covariates: batch, axiom_plate
# + RUID info
metadata4 <- fread("ABCD/ABCD_genotype/ABCD_release3.0_.batch_info.txt", data.table = F, header = T, stringsAsFactors = F)
colnames(metadata4)[2] <- "src_subject_id"


# 5. metadata 5 (ancestry PCs, computed on genotype data)
# this metadata only includes indivs. that passed the QC on the genotype data
# (8,816 indivs. )
metadata5 <- read.delim("ABCD/bb926/genotype/abcd.eigenvec", h=T, sep="\t")
metadata5$X.FID <- NULL
colnames(metadata5)[1] <- "src_subject_id"


# 6. merge all metadata files
metadata <- merge(metadata2, metadata3, by = "src_subject_id")
metadata <- merge(metadata, metadata1, by = "src_subject_id")
stopifnot(identical(metadata$sex.x, metadata$sex.y))
metadata$sex.y <- NULL
colnames(metadata)[which(colnames(metadata) == "sex.x")] <- "sex"
# this indv (NDAR_INVP2G0PXCM) is missing from metadata1
metadata <- merge(metadata, metadata4, by = "src_subject_id")
metadata <- merge(metadata, metadata5, by = "src_subject_id")

# compute NAs per covariate
# round(sort(apply(metadata, 2, function(x){sum(is.na(x))/length(x)}), decreasing = T), 2)
# round(sort(apply(metadata, 1, function(x){sum(is.na(x))/length(x)}), decreasing = T), 2)


#----------------------------
# intersect with JL dataset |
#----------------------------

# 7. read in JL dataset

## 7.1. xception
jl.subset1 <- read.csv("ABCD/liability/Xception_liability_w.wo.CBCL_04062023.csv", h=T)
jl.subset1$X <- NULL

## 7.2. xgboost
jl.subset2 <- read.csv("ABCD/liability/XGBoost_liability_w.wo.CBCL_04272023.csv", h=T)
jl.subset2$X <- NULL

## 7.3. merge the two dataframes above
jl.subset <- merge(jl.subset1, jl.subset2, by = "subject")


# 8. list of individuals w/ wearables and genetic info
# n = 1,397
a <- intersect(metadata$src_subject_id, jl.subset$subject)
test <- intersect(metadata5$src_subject_id, jl.subset$subject) # sanity check against indivs. that pass QC for genotyped snps


# 9. get the final metadata file
# this is after removing individuals that have NAs in some covariates
metadata.final <- metadata[complete.cases(metadata), ] # inverting these two lines doesn't change the result
metadata.final <- metadata.final[, c(39, 1:38, 40:46)]


# 10. list of individuals w/ wearables, genetics and complete covariate info
# n = 1,191
b <- intersect(metadata.final$src_subject_id, jl.subset$subject)


# 11. checking the numbers of cases and controls
# for adhd, there are 137 cases and 1,042 controls
# for any psych disorder, there are 12 individuals that are both cases and controls
# for now we keep all of them
metadata.jl <- metadata.final[metadata.final$src_subject_id %in% b, ]
table(metadata.jl[, c("adhd_group", "nonclinical_controls")])

metadata.jl$case <- apply(metadata.jl[, 5:14], 1, sum)
metadata.jl$case <- ifelse(metadata.jl$case >= 1, 1, 0)
metadata.jl$ctrl <- metadata.jl$nonclinical_controls
table(metadata.jl[, c("case", "ctrl")])



#----------------------------------
# prepare final dataframe with    |
# response (binary/quantitative)  |
# and covariates                  |
#----------------------------------


metadata.jl$ID <- paste(metadata.jl$RUID, metadata.jl$src_subject_id, sep = "_")
metadata.jl$RUID <- metadata.jl$ID 

all.covs <- c("acs_raked_propensity_score", 
              "sex", 
              "interview_age",
              "BATCH",
              "PC1",
              "PC2",
              "PC3",
              "PC4",
              "PC5",
              # "cbcl_externalizing", 
              # "cbcl_internalizing", 
              # "dsm_externalizing", 
              # "dsm_internalizing", 
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

elasticnet.covs <- c("interview_age",
                     "sex",
                     "famhis_bip",
                     "famhis_schiz",
                     "famhis_nerves",
                     "famhis_treatment",
                     "famhis_hospital",
                     "parent_div_cat",
                     "parent_grade",
                     "family_income",
                     "BATCH", "PC2", "PC4")


PGC.covs <- c("interview_age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5")



# create covariates dataframe

## all covariates
covariates.allCovs.df <- metadata.jl[, colnames(metadata.jl) %in% c("RUID", all.covs)]
colnames(covariates.allCovs.df)[1] <- "#IID"
# write.table(covariates.allCovs.df, file=paste0(outfolder, "/covariates.allCovs.tsv"), 
#             row.names = F, col.names = T, sep = "\t", quote = F)

## elasticnet covariates
covariates.elasticNet.df <- metadata.jl[, colnames(metadata.jl) %in% c("RUID", elasticnet.covs)]
colnames(covariates.elasticNet.df)[1] <- "#IID"
# write.table(covariates.elasticNet.df, file=paste0(outfolder, "/covariates.elasticNet.tsv"), 
#             row.names = F, col.names = T, sep = "\t", quote = F)

## PGC covariates
covariates.PGC.df <- metadata.jl[, colnames(metadata.jl) %in% c("RUID", PGC.covs)]
colnames(covariates.PGC.df)[1] <- "#IID"
# write.table(covariates.PGC.df, file=paste0(outfolder, "/covariates.PGC.tsv"), 
#             row.names = F, col.names = T, sep = "\t", quote = F)




# create response dataframe
# this df contains both the binary and quantitative response
response.df <- merge(metadata.jl[, c("RUID", "src_subject_id", 
                                     "adhd_group", "cbcl_externalizing", "cbcl_internalizing")], 
                     jl.subset, by.x = "src_subject_id", by.y = "subject")
response.df$src_subject_id <- NULL
colnames(response.df)[1] <- "#IID"

# convert cases-ctrls to 1/2
response.df$adhd_group <- ifelse(response.df$adhd_group == 0, 1, 2)
table(response.df$adhd_group)



# generate new liability scores setting cases to 1 (JL suggestion)

## xception
response.df$liability_Xception_without_cbcl_v2 <- ifelse(response.df$adhd_group == 2, 1, response.df$liability_Xception_without_cbcl)
response.df$liability_Xception_with_CBCL_v2 <- ifelse(response.df$adhd_group == 2, 1, response.df$liability_Xception_with_CBCL)

## xgboost
response.df$XGB_without_cbcl_v2 <- ifelse(response.df$adhd_group == 2, 1, response.df$XGB_without_cbcl)
response.df$XGB_with_cbcl_v2 <- ifelse(response.df$adhd_group == 2, 1, response.df$XGB_with_cbcl)


# save file
# write.table(response.df, file=paste0(outfolder, "/phenotypes.tsv"), row.names = F, col.names = T, sep = "\t", quote = F)

