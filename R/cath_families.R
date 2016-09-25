library(biomaRt)

cath_id   N
133 1.20.1070.10 539
106  1.10.510.10 213
675  3.40.50.300 190
347   2.60.40.10 170
257  2.130.10.10 164
724   3.80.10.10 133
222   1.25.10.10 127
50   1.10.238.10 117
154 1.20.1250.20 115
559  3.30.70.330 112
667 3.40.50.1820  92
225   1.25.40.10  89
663  3.40.50.150  68
230   1.25.40.20  66
296   2.30.30.40  64
521   3.30.40.10  62
360  2.60.40.150  60
680  3.40.50.720  58
397  3.10.100.10  56
754  3.90.190.10  56
44   1.10.20.10 53
116 1.10.630.10 53
242  1.50.40.10 51
39  1.10.150.50 44
302  2.30.42.10 44
112 1.10.565.10 39
110 1.10.555.10 38
204 1.20.58.390 38
574 3.30.710.10 37
794 4.10.280.10 37

subset_cath <- subset(cath_all, cath_id == "1.20.1070.10")
subset_cors <- subset(cath.cors, cath_id == "1.20.1070.10")
nrow(subset_cors)

family_ids <- unique(c(subset_cors$id_1, subset_cors$id_2))
length(family_ids)
subset(cath_all, stable_id %in% family_ids & omega > 1)

pdf("cors_hist.pdf")
hist(-log10(as.numeric(subset_cors$pval)), breaks=100)
hist(as.numeric(subset_cors$cor), breaks=100)
dev.off()

head(subset_cors)
subset_cors <- ddply(subset_cors, c("cath_id", "id_1", "id_2", "cor", "pval"), function(df) {
  data.frame(dataset_1=as.character(id2dataset[[ df$id_1[1] ]]), dataset_2=as.character(id2dataset[[ df$id_2[1] ]]), stringsAsFactors=F)
})


cath.cors.diff <- ddply(subset_cors, c("id_1", "id_2"), function(r) {
  if (strsplit(r$dataset_1[1], "_", fixed=T)[[1]][1] != strsplit(r$dataset_2[1], "_", fixed=T)[[1]][1]) {
    return(data.frame())
  } else {
    return(r)
  }
})

cath.cors.same <- ddply(subset_cors, c("id_1", "id_2"), function(r) {
  if (strsplit(r$dataset_1[1], "_", fixed=T)[[1]][1] == strsplit(r$dataset_2[1], "_", fixed=T)[[1]][1]) {
    return(data.frame())
  } else {
    return(r)
  }
})

hist(cath.cors.same$pval, breaks=100)
hist(cath.cors.diff$pval, breaks=100)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_peptide_id", "go_id", "name_1006", "definition_1006", "go_linkage_type", "namespace_1003"),
                 filters = list(ensembl_peptide_id=family_ids),
                 mart = mart)
results <- subset(results, namespace_1003 == "molecular_function")

evidence_types <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP")
sum(results$go_linkage_type %in% evidence_types)

results_subset <- subset(results, go_linkage_type %in% evidence_types)

id_combs <- data.frame(t(combn(family_ids, 2, simplify=T)))
colnames(id_combs) <- c("id_1", "id_2")

common_go <- function(r) {
  length(intersect(subset(results_subset, ensembl_peptide_id==r$id_1[1])$name_1006,
            subset(results_subset, ensembl_peptide_id==r$id_2[1])$name_1006)) / 
    length(union(subset(results_subset, ensembl_peptide_id==r$id_1[1])$name_1006,
                     subset(results_subset, ensembl_peptide_id==r$id_2[1])$name_1006))
}
go_dists <- ddply(subset_cors, c("cath_id", "id_1", "id_2", "cor", "pval", "dataset_1", "dataset_2"), common_go)
plot(go_dists$V1, -log(go_dists$pval))
plot(go_dists$V1, go_dists$cor)

# For Cytoscape
subset_cors$p.adjusted <- p.adjust(subset_cors$pval, method="BH")
subset_cors_sign <- subset(subset_cors, p.adjusted < 1e-05)

cors_sif <- data.frame(P1=subset_cors_sign$id_1, int="cor", P2=subset_cors_sign$id_2)
write.table(cors_sif, file="data/cath/cytoscape/1.20.1070.10.sif", sep="\t", quote=F, row.names=F)

cors_ann <- data.frame(int=paste(subset_cors$id_1, "(cor)", subset_cors$id_2), pval=subset_cors$pval)
write.table(cors_ann, file="data/cath/cytoscape/1.20.1070.10_ann.tab", sep="\t", quote=F, row.names=F)

# Number of occurrences of each GO term
go_numbers <- ddply(results, "name_1006", function(df) {
  data.frame(n_genes=length(unique(df$ensembl_peptide_id)))
})
go_numbers <- go_numbers[order(go_numbers$n_genes, decreasing=T), ]

all_go

unique(results$go_linkage_type)
