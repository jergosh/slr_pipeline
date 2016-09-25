library(biomaRt)
library(org.Hs.eg.db)
library(topGO)
library(GO.db)

all_genes_names <- unique(slr_all$stable_id)
all_genes <- as.integer(!(all_genes_names %in% pdb_master_globular$stable_id))
names(all_genes) <- all_genes_names
all_genes <- all_genes[names(all_genes) %in% all_tm]

ensMart <- useMart("ensembl", "hsapiens_gene_ensembl")
ens_map <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), filters =
                      "ensembl_peptide_id", values = names(all_genes), mart = ensMart)

nrow(ens_map) == length(all_genes)
names(all_genes) <- ens_map$ensembl_gene_id

GOdata <- new("topGOdata", ontology = "CC", allGenes = all_genes,
              geneSel = function(p) p < 1e-2, description = "Test", annot = annFUN.org,
              mapping="org.Hs.eg.db", ID="Ensembl")

resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

GenTable(GOdata, classicFisher = resultFisher, topNodes = 100)

