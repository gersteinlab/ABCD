

#************
# LIBRARIES *
#************

library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

#********
# BEGIN *
#********

# 1.compute correlation between phenotype features 
#we want to define clusters of correlated features   
m <- read.delim("v1.v2/phenotypes.tsv")
m <- m[, 2:259]
x <- as.data.frame(cor(m))


# 2.perform hierarchical clustering on the correlation matrix 
dist <- dist(x, diag=TRUE)
hc <- hclust(dist, method = "complete")


# 3. elbow plot on kmeans clusters
set.seed(123)
res <- c()
for (k in 1:20){
  res <- rbind(res, cbind(k, J = kmeans(x, centers = k, nstart = 10)$tot.withinss))
}
ggplot(as.data.frame(res), aes(k, -J)) +
  geom_path() +
  geom_point() +
  scale_x_continuous(breaks = 1:20) +
  theme_bw()


# 4. extract 7 clusters
annotation_clusters <- data.frame(cluster = cutree(hc, k = 7))
rownames(annotation_clusters) <- rownames(x)
annotation_clusters$cluster <- as.factor(annotation_clusters$cluster)

df.clusters <- data.frame(stringsAsFactors = F)
cl <- c(2, 7, 5, 6, 3, 4, 1)
start = 1
for (k in cl) {

  l = nrow(annotation_clusters[annotation_clusters$cluster == k, , drop=F ])
  end = start + l - 1
  df.clusters <- rbind(df.clusters,
                       data.frame(feature = colnames(x)[hc$order[start:end]],
                                  cluster = rep(k, l)))
  start = end + 1
}


# 5. save features per cluster matrix
write.table(df.clusters, file = "v1.v2.7clusters/features.clusters.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)


# 6. heatmap of correlation of phenotype features across clusters
pheatmap(x, show_rownames = F, show_colnames = F,
         clustering_method = "complete",
         clustering_distance_rows = dist,
         clustering_distance_cols = dist,
         annotation_row = annotation_clusters,
         annotation_col = annotation_clusters,
         annotation_colors = list(cluster = c("1" = "#F8766D",
                                              "2" = "#C49A00",
                                              "3" = "#53B400",
                                              "4" = "#00C094",
                                              "5" = "#00B6EB",
                                              "6"  = "#A58AFF",
                                              "7" =  "#FB61D7")),
         main = "Pearson r between features (v1+v2)")

