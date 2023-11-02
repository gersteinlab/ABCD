# 1. read batch info metadata 
batch.info <- read.delim("/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/ABCD_genotype/ABCD_release3.0_.batch_info.txt", h=T, sep="\t")


# 2. read covariates metadata
df <- read.delim("/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/bb926/metadata/acspsw03.txt", h=T, sep="\t")

## 2.1. removing second header from metadata
df <- df[2:nrow(df), ]
## 2.2. keep info only for individuals for which we have genotype info
df <- df[df$src_subject_id %in% batch.info$abcd.id_redcap, ]
## 2.3. we will keep only rows corresponding for baseline_year
df <- df[df$eventname == "baseline_year_1_arm_1", ]


# 3. read results of gender check performed by plink
gender <- read.table("/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/bb926/genotype/QC/plink.sexcheck", header=T,as.is=T)


# 4. compare gender detected by plink with gender present in metadata
# 35 indivs. are missing
gender <- merge(gender, df[, c("src_subject_id", "sex")], by.x = "IID", by.y = "src_subject_id")

# 5. all indivs. with available metadata 
# have metadata gender matching with genotyped gender
table(gender[, c("SNPSEX", "sex")])