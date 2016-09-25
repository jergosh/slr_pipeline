require(ReactomePA)
library(biomaRt)
library(GOstats)

mart <- useMart("ensembl", "hsapiens_gene_ensembl")

genes <- getBM(c("ensembl_gene_id", "entrezgene"), filters="ensembl_peptide_id", values=interesting_genes, mart)$entrezgene
universe_genes <- getBM(c("ensembl_gene_id", "entrezgene"), filters="ensembl_peptide_id", values=unique(slr_all$stable_id), mart)$entrezgene
x <- enrichPathway(gene = genes[!is.na(genes)], pvalueCutoff = thr, readable = T, universe=universe_genes)
summary(x)

top_genes_ann <- getBM(c("ensembl_gene_id", "ensembl_peptide_id", "description"),
                       filters="ensembl_peptide_id", values=top_genes, mart)

rownames(top_genes_ann) <- top_genes_ann$ensembl_gene_id


GOparams <- new("GOHyperGParams", geneIds=genes, 
                universeGeneIds=universe_genes,
                annotation="org.Hs.eg.db", ontology="MF", pvalueCutoff=0.05,
                conditional=TRUE, testDirection="over" )

GO <- hyperGTest(GOparams)

write.table(file=stdout(), cbind(summary(GO)$GOMFID, summary(GO)$Pvalue), quote=F, sep="\t", row.names=F)

## Use REVIGO to summarise
