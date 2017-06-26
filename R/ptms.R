library(ReactomePA)
library(Biostrings)
library(plyr)
library(biomaRt)

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
                      host="feb2014.archive.ensembl.org")

PTM_table <- read.table("data/PTMfunc/PTM_lists/psites_spectral_counts.tab", header=T, fill=T, stringsAsFactors=F)
PTM_table$PTMpos = PTM_table$PTMpos + 1 # Coordinate offset
PTM_human <- subset(PTM_table, Species == "HS")
nrow(PTM_human)

head(subset(PTM_human, is.na(ExperimentCounts)))
sum(is.na(PTM_human$ExperimentCounts))
# Filter out the rows with missing peptides -- for now
PTM_human <- subset(PTM_human, !is.na(PTM_human$ExperimentCounts))

w = 5 # Window size
common_ids <- intersect(slr_all$stable_id, PTM_human$ProteinID)
seqs <- getBM(attributes=c("coding", "ensembl_peptide_id"), filters="ensembl_peptide_id", values=unique(common_ids), mart=mart)

aa_seqs <- list()
for (stable_id in common_ids) {
  ptms <- subset(PTM_human, ProteinID == stable_id)
  if (! nrow(ptms)) {
    next
  }
  if (!stable_id %in% seqs$ensembl_peptide_id) {
    next
  }
  aa_seqs[[stable_id]] <- as.character(translate(DNAString(subset(seqs, ensembl_peptide_id == stable_id)[1, "coding"])))
}


# Filter out the proteins which are not present in the SLR dataset
PTM_human_filtered <- subset(PTM_human, ProteinID %in% common_ids)

rownames(seqs) <- seqs$ensembl_peptide_id

# seqs_matching <- seqs[PTM_human_filtered$ProteinID, ]

ens_context <- apply(PTM_human_filtered, 1, function(x) { pos = as.numeric(x["PTMpos"]); substr(aa_seqs[[x["ProteinID"]]], pos-w, pos+w) })
PTM_human_matching <- subset(PTM_human_filtered, ens_context == Peptide)

# ptm_nbhood <- ddply(PTM_human_matching, c("ProteinID", "PTMpos"),
#                     function(x) { slr_all[slr_all$stable_id == x$ProteinID, ][(x$PTMpos-w):(x$PTMpos+w), ] }) 

slr_split <- split(slr_all, slr_all$stable_id)
ptm_nbhood <- ddply(PTM_human_matching, c("ProteinID", "PTMpos"),
                    function(x) { slr_split[[x$ProteinID]][(x$PTMpos-w):(x$PTMpos+w), ] }) 

# ptm_nbhood_doughnut <- ddply(PTM_human_matching, c("ProteinID", "PTMpos"),
#                     function(x) { slr_split[[x$ProteinID]][c((x$PTMpos-w):(x$PTMpos-1), (x$PTMpos+1):(x$PTMpos+w)), ] }) 

plot(density(slr_all$Omega, bw=0.005))
# lines(density(ptm_nbhood_doughnut$Omega), col="red")
lines(density(ptm_nbhood$Omega, bw=0.005), col="green")
dim(ptm_nbhood)
hist(slr_all$Omega, freq=F, breaks=200)
hist(ptm_nbhood$Omega, add=T, freq=F, breaks=200, col="green")

m <- matrix(c(sum(ptm_nbhood$Adj.Pval < thr & ptm_nbhood$Omega > 1), sum(pos_sel$Adj.Pval < thr & pos_sel$Omega > 1), nrow(ptm_nbhood), nrow(slr_all)), nrow=2)

fisher.test(m, alternative="greater")

# Enrichment

genes <- getBM(c("ensembl_gene_id", "entrezgene"), filters="ensembl_peptide_id",
               values=unique(subset(ptm_nbhood, Adj.Pval < thr & Omega > 1)$stable_id), mart)$entrezgene
ptm_universe <- getBM(c("ensembl_gene_id", "entrezgene"), filters="ensembl_peptide_id",
                      values=unique(PTM_human_matching$ProteinID), mart)$entrezgene

x <- enrichPathway(gene = genes, pvalueCutoff = 0.05, readable = T, universe=ptm_universe)
summary(x)
str(x@geneInCategory)
