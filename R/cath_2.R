library(gplots)

all_ids <- unique(subset_cath$stable_id)
cormat <- matrix(, nrow=length(all_ids), ncol=length(all_ids))
rownames(cormat) <- all_ids
colnames(cormat) <- all_ids

for (i in 1:nrow(subset_cors)) {
  cormat[subset_cors[i, "id_1"], subset_cors[i, "id_2"]] <- 1-subset_cors[i, "cor"]
  cormat[subset_cors[i, "id_2"], subset_cors[i, "id_1"]] <- 1-subset_cors[i, "cor"]  
}
diag(cormat) <- 0
sum(is.na(cormat))
cormat[is.na(cormat)] <- 0

cluster <- hclust(as.dist(cormat), method = "ward.D2")

annotated <- subset(results, name_1006 == "detection of chemical stimulus involved in sensory perception of smell" |
                      name_1006 == "sensory perception of taste")$ensembl_peptide_id
colours <- rep("black", length(all_ids))
colours[all_ids %in% annotated] <- "red"

annotated_taste <- subset(results, name_1006 == "detection of chemical stimulus involved in sensory perception of taste")$ensembl_peptide_id

for (i in 1:nrow(subset_cors)) {
  cormat[subset_cors[i, "id_1"], subset_cors[i, "id_2"]] <- subset_cors[i, "cor"]
  cormat[subset_cors[i, "id_2"], subset_cors[i, "id_1"]] <- subset_cors[i, "cor"]  
}
diag(cormat) <- 1
heatmap.2(cormat,
          # symm=T,
          trace="none",
          col=redgreen(75),
          # ColSideColors=colours,
          # RowSideColors=colours,
          Rowv = as.dendrogram(cluster),
          Colv = as.dendrogram(cluster))

other_ids <- all_ids[!(all_ids %in% annotated)]
print(paste(other_ids, collapse=" "))

cluster_colours <- rep("black", length(all_ids))
cluster_ids <- cutree(cluster, k=4)[all_ids]
names(cluster_ids) == all_ids
cluster_colours[cluster_ids == 2 | cluster_ids == 1] <- "blue"
heatmap.2(cormat,
          symm=T,
          trace="none",
          ColSideColors=cluster_colours,
          RowSideColors=colours,
          Rowv = as.dendrogram(cluster),
          Colv = as.dendrogram(cluster))

print(paste(all_ids[cluster_ids == 2 | cluster_ids == 1], collapse=" "))

setdiff(all_ids[cluster_ids == 2 | cluster_ids == 1], all_ids[all_ids %in% annotated])
"ENSP00000314992" %in% all_ids
