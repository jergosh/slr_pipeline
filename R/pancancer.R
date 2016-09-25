library(biomaRt)
library(dplyr)

ensembl.db <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
                      host="aug2014.archive.ensembl.org")
grep("hgnc", listAttributes(ensembl.db)[,1], value=T)

pancan <- read.table("~/Downloads/PanCan.maf", sep="\t", stringsAsFactors=F, header=T)

id.mapping <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol",
                     values=unique(pancan$gene), mart=ensembl.db)
id.mapping$gene <- id.mapping$hgnc_symbol
sum(unique(pancan$gene) %in% id.mapping$hgnc_symbol)
length(unique(pancan$gene))

pancan.ids <- inner_join(pancan, id.mapping, type="inner", copy=TRUE)
pancan$stable_id <- rep(NA, nrow(pancan))

pancan$stable_id[grep("ENSG0", pancan$gene)] <- grep("ENSG0", pancan$gene, value=T)

# Need to obtain gene start position and strand from Biomart

# If we get the sequence from Biomart we can verify if the sequences match
# Having the NT position witih the gene, we can work out the AA position
# Additionally, having the sequence, I can work out what the codons are and verify
# if the amino-acids agree as well