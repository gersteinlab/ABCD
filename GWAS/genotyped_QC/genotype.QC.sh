#~~~~~~~~~~~~~~~~~~~~~~~
# This QC was performed following the tutorial
# at: https://github.com/MareesAT/GWA_tutorial/blob/master/1_QC_GWAS.zip
# (see publication: https://pubmed.ncbi.nlm.nih.gov/29484742/)
#~~~~~~~~~~~~~~~~~~~~~~~

############################
### Step 1 ###

bfile='ABCD/ABCD_genotype/genotype_QCed/ABCD_release_3.0_QCed'

#*****
# 1.1. Investigate missingness per individual and per SNP and make histograms.
plink --bfile $bfile --missing
#*****

#*****
# 1.2. Delete SNPs and individuals with high levels of missingness, explanation of this and all following steps can be found in box 1 and table 1 of the article mentioned in the comments of this script.
# The following two QC commands will not remove any SNPs or individuals. However, it is good practice to start the QC with these non-stringent thresholds.  

# 1.2.1. Delete SNPs with missingness >0.2.
plink --bfile $bfile --geno 0.2 --make-bed --out res_2

# 1.2.2. Delete individuals with missingness >0.2.
plink --bfile res_2 --mind 0.2 --make-bed --out res_3

# 1.2.3. Delete SNPs with missingness >0.02.
plink --bfile res_3 --geno 0.02 --make-bed --out res_4

# 1.2.4. Delete individuals with missingness >0.02.
plink --bfile res_4 --mind 0.02 --make-bed --out res_5
#*****


############################
### Step 2 ###
# Check for sex discrepancy.
# Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. 
# This F value is based on the X chromosome inbreeding (homozygosity) estimate.
# Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.
# NOTE: phenotypes are missing in our data, so all indivs. will be flagged as "PROBLEM"

#*****
# 2.1. infer sex
plink --bfile res_5 --check-sex 
#*****

#*****
# 2.2. Check consistency between gender declared in the metadata and 
# gender detected by PLINK (35/10660 indivs. don't have available metadata)
Rscript check.sex.R
#*****

#*****
# 2.3. Because all indivs. w/ available metadata 
# show consistent gender, we will impute gender for the 35 indivs.
# w/ missing metadata info
plink --bfile res_5 --impute-sex --make-bed --out res_6
#*****



############################
### Step 3 ###
# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).

#*****
# 3.1. Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' res_6.bim > snp_1_22.txt
plink --bfile res_6 --extract snp_1_22.txt --make-bed --out res_7
#*****

#*****
# 3.2. Remove SNPs with a low MAF frequency.
# 427704 SNPs are left
plink --bfile res_7 --maf 0.01 --make-bed --out res_8
#*****



############################
### Step 4 ###

#*****
# 4.1. Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
# Check the distribution of HWE p-values of all SNPs.
plink --bfile res_8 --hardy
#*****

#*****
# 4.2. By default the --hwe option in plink only filters for controls.
# Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data.
# The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE. 
# This second HWE step only focusses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed
# NOTE: phenotypes are not specified in our data,
# so plink cannot distinguish between case and controls
# we apply a more relaxed threshold to all snps
plink --bfile res_8 --hwe 1e-10 --hwe-all --make-bed --out res_9
#*****



############################
### Step 5 ###

# Generate a plot of the distribution of the heterozygosity rate of your subjects.
# And remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

#*****
# 5.1. Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly) correlated SNPs, we prune the SNPs using the command --indep-pairwise
# The parameters 50 5 0.2 stand respectively for: 
# the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.
plink --bfile res_9 --indep-pairwise 50 5 0.2 --out indepSNP
#*****

#*****
# 5.2. check heterozyg. rate for every indiv. based on pruned snps 
plink --bfile res_9 --extract indepSNP.prune.in --het --out R_check
#*****

#*****
# 5.3. Plot the heterozygosity rate distribution
# and generate a list of individuals who deviate more than 3 standard deviations from the mean het rate
# 193 indivs. are filtered out
Rscript check.hetRate.R
#*****

#*****
# 5.4. Output of the command above: fail-het-qc.txt .
# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print $1, $2}'> het_fail_ind.txt
#*****

#*****
# 5.5. Remove heterozygosity rate outliers.
plink --bfile res_9 --remove het_fail_ind.txt --make-bed --out res_10
#*****



############################
### Step 6 ###

# It is essential to check datasets you analyse for cryptic relatedness
# We are going to exclude all individuals above the pihat threshold of 0.2

#******
# 6.1. Check for relationships between individuals with a pihat > 0.2
# This step takes ~10'
# The file 'pihat_min0.2.genome' contains individuals with pihat > 0.2
plink --bfile res_10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
#******

#******
# 6.2. Double-check that the groups of related indivs. are consistent with the info reported in the metadata

# 6.2.1. groups of related indivs. according to pihat
awk '{print $2"\t"$4}' pihat_min0.2.genome | tail -n+2 | awk 'BEGIN{FS=OFS="\t"}{if (!($2 in a)){print $0; a[$2]="checked"; a[$1]="checked"}}' | ~/bin/make.column.list.py > pihat.0.2.related.indivs.tsv
awk '{print $1";"$2}' pihat.0.2.related.indivs.tsv | awk '{split($0, a, ";"); asort(a); for(i=1;i<length(a);i++) printf("%s;",a[i]); printf(a[length(a)]"\n")}' > tmp; mv tmp pihat.0.2.related.indivs.tsv

# 6.2.2. groups of related indivs. according to metadata
sed 's/\"//g' ../../ancestry/acspsw03.txt | awk 'BEGIN{FS=OFS="\t"}$9=="baseline_year_1_arm_1"{print $5, $30, $31, $32, $17}' | awk 'BEGIN{FS=OFS="\t"}{if (!($1 in b)){b[$1]="checked"; for (i=2;i<=NF;i++){if ($i!=""){if (!($1 in a)){a[$1]=$i} else{a[$1]=a[$1]";"$i}; if (!($i in b)){b[$i]="checked"}}}}}END{for (k in a){print k, a[k]}}' > metadata.related.indivs.tsv
awk '{print $1";"$2}' metadata.related.indivs.tsv | awk '{split($0, a, ";"); asort(a); for(i=1;i<length(a);i++) printf("%s;",a[i]); printf(a[length(a)]"\n")}' > tmp; mv tmp metadata.related.indivs.tsv 

# NOTE: 1,506 groups of related individs consistent between pihat calculation and metadata info
# 88 groups of related indivs. found by pihat that are inconsistent with metadata info
# multiple reasons: there are groups annotated in the metadata as a family but then they don't have a pihat, and the computed pihat is not > 0.2
# other reasons: annotation of families sometimes doesn't incorporate in the same nucleus couples of cousins
#******

#******
# 6.3. For each group of 'related' individuals with a pihat > 0.2, we recommend to remove the individual(s) with the lowest call rate. 

# 6.3.1. re-compute missingness per individual
plink --bfile res_10 --missing

# 6.3.2. There are 1,594 groups of related individuals; 
# there's a problematic indivd. (XXX anonymized for GitHub) that is related to two different people
# which are unrelated to each other (cousin of a cousin?) - w/o this indiv. we have 1,592 groups of related indivds.
wc -l pihat.0.2.related.indivs.tsv
grep -v "XXX" pihat.0.2.related.indivs.tsv | wc -l

# 6.3.3. These 1,594 (1,592) groups comprise 3,245 (3,242) indivs. 
cat pihat.0.2.related.indivs.tsv | while read group; do echo $group | awk '{n=split($1, a, ";"); for (i=1;i<=n;i++){print a[i]}}'; done | sort -u | wc -l
grep -v "XXX" pihat.0.2.related.indivs.tsv | while read group; do echo $group | awk '{n=split($1, a, ";"); for (i=1;i<=n;i++){print a[i]}}'; done | wc -l

# 6.3.4. Obtain a list of 1,651 (1,650) individuals to be removed
cat pihat.0.2.related.indivs.tsv | while read group; do grep -Ff <(echo $group | awk '{n=split($1, a, ";"); for (i=1;i<=n;i++){print a[i]}}') plink.imiss | awk '{print $2"\t"$6}' | sort -k2,2g | tail -n+2; done | sort -u | cut -f1 > 0.2_low_call_rate_pihat.txt
grep -v "XXX" pihat.0.2.related.indivs.tsv | while read group; do grep -Ff <(echo $group | awk '{n=split($1, a, ";"); for (i=1;i<=n;i++){print a[i]}}') plink.imiss | awk '{print $2"\t"$6}' | sort -k2,2g | tail -n+2; done | cut -f1 > test 
#******

#******
# 6.4. Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
grep -Ff 0.2_low_call_rate_pihat.txt res_10.fam | awk '{print $1"\t"$2}' > tmp; mv tmp 0.2_low_call_rate_pihat.txt
plink --bfile res_10 --remove 0.2_low_call_rate_pihat.txt --make-bed --out res_11
#******


##########################
### Step 7 #

#*******
# pruning was already performed as part of the QC (see command 5.1)
plink2 --bfile QC/res_11 --extract QC/indepSNP.prune.in --out abcd.pruned --make-bed --threads 1
#*******

#*******
# PCA based on genotypes
# "approx" can be a good idea when you have >5000 samples, and is almost required once you have >50000
# we have 8,816
plink2 --bfile abcd.pruned --pca 5 "approx" --out ./abcd --threads 1
#*******

