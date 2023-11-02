#************
# LIBRARIES *
#************

library(data.table)
library(tidyverse)


#********
# BEGIN *
#********

# 1. set wd
setwd("ABCD/bb926/metadata/wr.folder")

# 2. read tables
df_fitbit <- fread("abcd_fbdpas01.txt", stringsAsFactors = F, data.table = F)
df_fitbit <- df_fitbit[-1, ]

df_ksad <- fread("abcd_ksad01.txt", stringsAsFactors = F, data.table = F) ## read in ksad datafile 
df_ksad <- df_ksad[-1, ]

df_famhis <- fread("fhxp201.txt", stringsAsFactors = F, data.table = F) ##read in family history datafile
df_famhis <- df_famhis[-1, ]

df_cbcl <- fread("abcd_cbcls01.txt", stringsAsFactors = F, data.table = F) ##read in the CBCL dataset
df_cbcl <- df_cbcl[-1, ]

df_demo <- fread("pdem02.txt", stringsAsFactors = F, data.table = F) ###demographic
df_demo <- df_demo[-1, ]

df_sst <- fread("abcd_sst02.txt", stringsAsFactors = F, data.table = F) ### stop signal task performance
df_sst <- df_sst[-1, ]

df_nihtb <- fread("abcd_tbss01.txt", stringsAsFactors = F, data.table = F) ### cognitive performance
df_nihtb <- df_nihtb[-1, ]

df_demo2 <- fread("acspsw03.txt", stringsAsFactors = F, data.table = F) ### other demographic questionnaire
df_demo2 <- df_demo2[-1, ]

# 3. filter data to only include observations from year 2 follow_up_arm_1
time_filter <- function(x){
  filtered <- filter(x, eventname == '2_year_follow_up_y_arm_1')
  return(filtered)
}

df_fitbit <- time_filter(df_fitbit)
df_ksad <- time_filter(df_ksad)
df_cbcl <- time_filter(df_cbcl)
df_sst <- time_filter(df_sst)
df_nihtb <- time_filter(df_nihtb)

# 4. pull family history from baseline eval
df_famhis <- filter(df_famhis, visit == 'baseline_year_1_arm_1')
df_demo2 <- filter(df_demo2, eventname == '1_year_follow_up_y_arm_1')

df_list <- list(df_ksad, df_cbcl, df_demo, df_sst, df_nihtb, df_famhis, df_demo2)
df_union <- reduce(df_list,full_join, by = 'src_subject_id') 

# 5. create vector of unique subject identifiers and filter ksad 
fb_unique <- unique(df_fitbit$src_subject_id)
df_union <- filter(df_union, src_subject_id %in% fb_unique)

internalizing_cluster <- c('ksads_1_843_p', 
                           'ksads_1_844_p',
                           'ksads_1_840_p', 
                           'ksads_1_841_p',
                           'ksads_1_846_p', 
                           'ksads_5_857_p', 
                           'ksads_5_906_p', 
                           'ksads_5_857_p',
                           'ksads_6_908_p',
                           'ksads_6_859_p', 
                           'ksads_7_861_p',
                           'ksads_7_909_p',
                           'ksads_8_863_p',
                           'ksads_8_911_p',
                           'ksads_9_868_p',
                           'ksads_10_913_p', 
                           'ksads_10_869_p',
                           'ksads_25_915_p',
                           'ksads_25_865_p'
)

externalizing_cluster <- c('ksads_14_856_p',
                           'ksads_14_855_p',
                           'ksads_14_853_p',
                           'ksads_15_901_p',
                           'ksads_16_897_p',
                           'ksads_16_898_p'
)

# 6. create diagnostic clusters for classification testing. All have corresponding genetic correlates per 
group_define <- function(cluster){
  group <- apply(df_union[cluster], 1, function(x) any(x %in% 1)) %>% as.numeric
  return(group)
}

# 7. behavioral phenotypes from the diagnostic perspective 
df_union$internalizing_dx <- group_define(internalizing_cluster)
df_union$externalizing_dx <- group_define(externalizing_cluster)

# 8. extracting important variables and dropping the rest
df_union <- select(df_union, 
                   src_subject_id = src_subject_id,
                   ###STARTING DIAGNOSTIC INFO SECTION
                   cbcl_externalizing = cbcl_scr_syn_external_t,##continuous measure externalizing symptoms
                   cbcl_internalizing = cbcl_scr_syn_internal_t, ##continuous measure internalizing symptoms
                   dsm_externalizing = externalizing_dx, ##binary indicator for externalizing diagnosis
                   dsm_internalizing = internalizing_dx, ##binary indicator for internalizing diagnosis
                   ###STARTING COGNITIVE FUNCTIONING VARIABLES
                   stopsigreactiontime = tfmri_sst_all_beh_total_meanrt, ##inhibitory reaction time
                   stopsignalgort = tfmri_sst_all_beh_crgo_mrt, ##activation reaction time
                   standevgort = tfmri_sst_all_beh_crgo_stdrt, ## standard deviation of reaction time
                   picvocab = nihtbx_picvocab_agecorrected, ##picture vocabulary scaled score
                   flanker = nihtbx_flanker_agecorrected, ##flanker task scaled score
                   processspeed = nihtbx_pattern_agecorrected, ##processing speed scaled score
                   picmemory = nihtbx_picture_agecorrected, ##picture sequence memory test scaled score
                   readingscore = nihtbx_reading_agecorrected, ##oral reading recongition test scaled score
                   ##START DEMOGRAPHIC/PARENTAL COVARIATES
                   famhis_bip = fam_history_7_yes_no, ##any reported family history of mania
                   famhis_schiz = fam_history_8_yes_no, ##any reported family history of schizophrenia
                   famhis_antisocial = fam_history_9_yes_no, ##any reported family history of antisocial behavior
                   famhis_nerves = fam_history_10_yes_no, ##any family history of anxiety disorder or "nervous experiences"
                   famhis_treatment = fam_history_11_yes_no, ##any family history of psych treatment
                   famhis_hospital = fam_history_12_yes_no, ##any history of psych hospitalization in family
                   famhis_suicide = fam_history_13_yes_no, ##any history of suicide attempt in family
                   parent_div_cat = demo_prnt_marital_v2, ### marital status of parents ##CATEGORICAL##
                   parent_grade = demo_prnt_ed_v2, ### highest grade completed respondent parent
                   family_income = demo_comb_income_v2, ##combined income for full family
                   race_eth_cat = race_ethnicity, ##categorical race/eth ###CATEGORICAL###
                   ####ADD ADOPTED INDICATOR VARIABLE
                   adopted = demo_prim
                   
)

# 9. convert character columns to integer / numeric columns
df_union[, 2:3] <- apply(df_union[, 2:3], 2, function(x){x <- as.integer(x)})
df_union[, 6:13] <- apply(df_union[, 6:13], 2, function(x){x <- as.numeric(x)})
df_union[, 14:25] <- apply(df_union[, 14:25], 2, function(x){x <- as.integer(x)})

for ( i in 1:ncol(df_union) ) {print(class(df_union[, i]))}

df_union_safe <- df_union

# 10. CONFIRM MISSINGNESS REMOVED
# replacing missing value indicators with missing 
df_union <- mutate(df_union, across(where(is.integer), ~recode(., 
                                                               `777` = NA_integer_, 
                                                               `999` = NA_integer_)
)
)


df_union <- mutate(df_union, across(starts_with('famhis_'), ~recode(., `7` = NA_integer_)))                  

# 11. recode adopted as indicator
df_union$adopted <- ifelse(df_union$adopted == 3, 1, 0)


# 12. save output
write.table(df_union, file = "clinical.covariates.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)                  

