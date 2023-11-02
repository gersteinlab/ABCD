
#************
# LIBRARIES *
#************

library(data.table)


#********
# BEGIN *
#********


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


# 2. metadata2 (JL metadata)
# covariates: interview_age, sex, cohorts groups
# need to make sure sex in metadata1 & 2 coincide
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


# 3. metadata3 (clinical covariates by WR)
# covariates: externalizing / internalizing gradients, family history
metadata3 <- fread("ABCD/bb926/metadata/wr.folder/clinical.covariates.tsv",
                   data.table = F, header = T, stringsAsFactors = F)
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
metadata4 <- fread("ABCD/ABCD_genotype/ABCD_release3.0_.batch_info.txt",
                   data.table = F, header = T, stringsAsFactors = F)
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
round(sort(apply(metadata, 2, function(x){sum(is.na(x))/length(x)}), decreasing = T), 2)
round(sort(apply(metadata, 1, function(x){sum(is.na(x))/length(x)}), decreasing = T), 2)

metadata <- metadata[complete.cases(metadata), ] # inverting these two lines doesn't change the result
metadata <- metadata[, c(39, 1:38, 40:46)]



# 7. save output
write.table(metadata, file = "ABCD/bb926/metadata/metadata.tsv",
            row.names = F, col.names = T, sep = "\t", quote = F)
