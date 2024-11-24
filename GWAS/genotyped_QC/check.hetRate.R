#********
# BEGIN *
#********

# 1. set working directory
setwd("/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/bb926/genotype/QC")

# 2. read heterozygosity rate
het <- read.table("R_check.het", head=TRUE)

# 3. compute distribution
pdf("../plots/genotype.QC.protocol/heterozygosity.pdf")
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
dev.off()

# 4. filter out snps that deviate more than 3 sd from the mean het rate
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)

