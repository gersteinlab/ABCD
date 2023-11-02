
#************
# LIBRARIES *
#************

library(readr)


#********
# BEGIN *
#********

# 1. read phenotypes.tsv file 
m <- read.delim("v1.v2/phenotypes.tsv", h=T, stringsAsFactors = F)


# 2. read featurs.tsv 
f <- read.delim("v1.v2.7clusters/features.clusters.tsv",  h=T, stringsAsFactors = F)


# 3. get the first cluster 
cluster1 <- subset(f, cluster=="1")
cluster1_pheno <- m %>% select(one_of(dput(as.character(cluster1$feature))))
cluster1_pheno$ID <- m$ID
cluster1_pheno <- cluster1_pheno[,c(ncol(cluster1_pheno),1:(ncol(cluster1_pheno)-1))]
write_tsv(cluster1_pheno, path = "v1.v2.7clusters/v1.v2.cluster1/cluster1.phenotypes.tsv")

# 4. get the second cluster 
cluster2 <- subset(f, cluster=="2")
cluster2_pheno <- m %>% select(one_of(dput(as.character(cluster2$feature))))
cluster2_pheno$ID <- m$ID
cluster2_pheno <- cluster2_pheno[,c(ncol(cluster2_pheno),1:(ncol(cluster2_pheno)-1))]
write_tsv(cluster2_pheno, path = "v1.v2.7clusters/v1.v2.cluster2/cluster2.phenotypes.tsv")


# 5. get the third cluster 
cluster3 <- subset(f, cluster=="3")
cluster3_pheno <- m %>% select(one_of(dput(as.character(cluster3$feature))))
cluster3_pheno$ID <- m$ID
cluster3_pheno <- cluster3_pheno[,c(ncol(cluster3_pheno),1:(ncol(cluster3_pheno)-1))]
write_tsv(cluster3_pheno, path = "v1.v2.7clusters/v1.v2.cluster3/cluster3.phenotypes.tsv")


# 6. get the forth cluster 
cluster4 <- subset(f, cluster=="4")
cluster4_pheno <- m %>% select(one_of(dput(as.character(cluster4$feature))))
cluster4_pheno$ID <- m$ID
cluster4_pheno <- cluster4_pheno[,c(ncol(cluster4_pheno),1:(ncol(cluster4_pheno)-1))]
write_tsv(cluster4_pheno, path = "v1.v2.7clusters/v1.v2.cluster4/cluster4.phenotypes.tsv")


# 7. get the fifth cluster 
cluster5 <- subset(f, cluster=="5")
cluster5_pheno <- m %>% select(one_of(dput(as.character(cluster5$feature))))
cluster5_pheno$ID <- m$ID
cluster5_pheno <- cluster5_pheno[,c(ncol(cluster5_pheno),1:(ncol(cluster5_pheno)-1))]
write_tsv(cluster5_pheno, path = "v1.v2.7clusters/v1.v2.cluster5/cluster5.phenotypes.tsv")


# 8. get the sixth cluster 
cluster6 <- subset(f, cluster=="6")
cluster6_pheno <- m %>% select(one_of(dput(as.character(cluster6$feature))))
cluster6_pheno$ID <- m$ID
cluster6_pheno <- cluster5_pheno[,c(ncol(cluster6_pheno),1:(ncol(cluster6_pheno)-1))]
write_tsv(cluster6_pheno, path = "v1.v2.7clusters/v1.v2.cluster6/cluster6.phenotypes.tsv")


# 9. get the 7th cluster 
cluster7 <- subset(f, cluster=="7")
cluster7_pheno <- m %>% select(one_of(dput(as.character(cluster7$feature))))
cluster7_pheno$ID <- m$ID
cluster7_pheno <- cluster7_pheno[,c(ncol(cluster7_pheno),1:(ncol(cluster7_pheno)-1))]
write_tsv(cluster7_pheno, path = "v1.v2.7clusters/v1.v2.cluster7/cluster7.phenotypes.tsv")

