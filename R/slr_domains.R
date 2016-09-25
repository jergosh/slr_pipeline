# TODO could calculate lengths by taking max(human_idx) 
slr_sign = subset(slr_all, Omega > 1 & Adj.Pval < thr)
slr_ranges = GRanges(seqnames=slr_sign$stable_id, ranges=IRanges(start=slr_sign$human_idx, end=slr_sign$human_idx))

# range_test = GRanges(seqnames="ENSP00000360158", ranges=IRanges(start=1, end=1200))
# countOverlaps(slr_ranges, range_test)

# Pfam domain stuff
feature_data <- read.table("~/Documents/projects/slr_pipeline/data/protein_features.tab", header=F, fill=T, sep="\t")
colnames(feature_data) <- c("stable_id", "start", "end", "type", "intepro_id", "domain_id", "domain_desc")

domain_overlaps <- function(feature_data, dtype) {
  domain_data <- subset(feature_data, type == dtype)
  domain_ranges <- with(domain_data, GRanges(seqnames=stable_id, ranges=IRanges(start=start, end=end)))
  rdomain_ranges <- reduce(domain_ranges)
  
  intersect(as.character(seqnames(slr_ranges)), as.character(seqnames(domain_ranges)))
  
  coverlaps <- countOverlaps(rdomain_ranges, slr_ranges)
  # overlaps <- findOverlaps(rdomain_ranges, slr_ranges)
  coverlaps
}

test_overlaps <- function(overlap, dtype) {
  domain_data <- subset(feature_data, type == dtype)
  print(head(domain_data))
  domain_ranges <- with(domain_data, GRanges(seqnames=stable_id, ranges=IRanges(start=start, end=end)))
  rdomain_ranges <- reduce(domain_ranges)
  
  N <- nrow(slr_all) # white + black
  n <- sum(width(rdomain_ranges)) # total white
  k <- sum(width(slr_ranges)) # total drawn
  x <- sum(overlap) # white drawn
  
  m <- matrix(c(x, k-x, n-x, N-(k+n)+x), nrow=2)
  print(m)
  # white drawn, total white, total black, total drawn
  # phyper(x, n, N-n, k, lower.tail=F)
  fisher.test(m, alternative="greater")
}

overlaps = list()
overlap_tests = list()
for (domain_type in unique(feature_data$type)) {
  overlaps[[domain_type]] <- domain_overlaps(feature_data, domain_type)
  overlap_tests[[domain_type]] <- test_overlaps(overlaps[[domain_type]], domain_type)
}

domain_data <- subset(feature_data, type == "pfam" & domain_id != "")
# Get domain breakdowns
# Possibilities include looking at the most common domain types
# Possibly excluding some promiscuous ones
# How do we decide which domains have more potential for adaptive evolution?
domain_counts <- daply(domain_data, "domain_id", nrow)
domain_counts <- domain_counts[!is.na(domain_counts)]
domain_ids <- names(domain_counts[domain_counts > 10])

domain_ranges <- with(domain_data, GRanges(seqnames=stable_id, ranges=IRanges(start=start, end=end)))
rdomain_ranges <- reduce(domain_ranges)
slr_domain_hits <- findOverlaps(slr_ranges, domain_ranges)
# slr_domain_ranges <- ranges(slr_domain_hits, slr_ranges, domain_ranges)

domain_type_tests <- list()

for (dom_id in domain_ids) {
  print(dom_id)
  domain_type_data <- subset(domain_data, domain_id == dom_id)
  domain_type_ranges <- with(domain_type_data, GRanges(seqnames=stable_id, ranges=IRanges(start=start, end=end)))
  rdomain_type_ranges <- reduce(domain_type_ranges)
  
  # intersect(as.character(seqnames(slr_ranges)), as.character(seqnames(domain_ranges)))
  
  coverlaps <- countOverlaps(rdomain_type_ranges, slr_ranges)
  print(sum(coverlaps))
  N <- sum(width(rdomain_ranges)) # white + black
  n <- sum(width(rdomain_type_ranges)) # total white
  k <- length(slr_domain_hits) # total drawn
  x <- sum(coverlaps) # white drawn
  
  m <- matrix(c(x, k-x, n-x, N-(k+n)+x), nrow=2)
  domain_type_tests[[dom_type]] <- fisher.test(m, alternative="greater")
  print(domain_type_tests[[dom_type]])
}

domain_type_pvals <- sapply(domain_type_tests, function(x) {x$p})
hist(domain_type_pvals, breaks=50)
domain_type_pvals_adj <- p.adjust(domain_type_pvals, "BH")
enriched_domain_types <- names(which(domain_type_pvals_adj < 0.05))
# Get the InterPro IDs for BiomaRt
enriched_domain_ids <- unique(subset(domain_data, domain_type %in% enriched_domain_types)$domain_id)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
enriched_domain_ann <- getBM(attributes=c("interpro", "interpro_description"), filters="interpro",
                             values=enriched_domain_ids, mart=mart)
# TODO: Add # of sites, number of proteins
enriched_domain_ann$Pval <- domain_type_pvals_adj[enriched_domain_types]
# enriched_domain_ann$NSites <- 


# Collagen looks suspicious
collagen_domains <- subset(domain_data, domain_type=="Collagen")
collagen_stable_ids <- unique(collagen_domains$stable_id)
collagen_pos <- subset(pos_sel, stable_id %in% collagen_stable_ids & Adj.Pval < thr)
collagen_pos_count <- daply(collagen_pos, "stable_id", nrow)
sort(collagen_pos_count)

dim(collagen_domains)
unique(collagen_domains$stable_id)


## Fisher's exact for individual genes
## NOTE: Is this really meaningful?
## Presumably it's better to combine pvalues using Fisher's method?
## But how do we do it only for 
n_pos <- sum(pos_sel$Adj.Pval < thr)
N <- nrow(slr_all)
genewise_fisher = list()
gene_lengths <- daply(slr_all, "stable_id", nrow)
for (gene in unique(slr_all$stable_id)) {
  if (is.na(pos_per_gene[gene])) {
    pos_per_gene[gene] <- 0
  }
  m <- matrix(c(pos_per_gene[gene], n_pos, gene_lengths[gene], N), nrow=2)
  genewise_fisher[[gene]] <- fisher.test(m, alternative="greater")
}

genewise_fisher_padjust <- p.adjust(sapply(genewise_fisher, function(x) { x$p}), "BH")
names(genewise_fisher_padjust) <- names(genewise_fisher)
sum(genewise_fisher_padjust < 0.05)
interesting_genes <- names(which(genewise_fisher_padjust < 0.05))

for (gene in interesting_genes) {
  print(pos_per_gene[gene])
}

hist(sapply(interesting_genes, function(x) { pos_per_gene[x] }), breaks=100)



extract_sequences <- function(df) {
  start <- df$start[1]
  end <- df$end[1]
  
  subset(slr_all, stable_id == df$stable_id[1] & human_idx >= start & human_idx <= end)
}

domain_omega <- ddply(subset(domain_data, domain_type == "VHS"),
                      .variables=c("stable_id", "start", "end"), extract_sequences)

ddply(domain_omega, .variables=c("stable_id"), function(df) { max(df$human_idx) - min(df$human_idx) })

ddply(domain_omega, .variables=c("stable_id"), function(df) { mean(df$Omega) })

ddply(domain_omega, .variables=c("stable_id"), function(df) { plot(df$human_idx, df$Omega, type="l") })

# Now we need an alignment
outdir <- "data/domains"

write.domain.df <- function(domain_df) {
  outfile <- paste(c(outdir, "/", as.character(domain_df$domain_id[1]), ".tab"), collapse="")
  print(outfile)
  domain_ann <- join(domain_df, slr_all, by="stable_id", match="all")
  domain_ann <- unique(domain_ann[, c("stable_id", "dataset", "domain_id", "start", "end")])
  write.table(domain_ann, file=outfile, quote=F, col.names=F, sep="\t", row.names=F)
}

ddply(domain_data, .variables=c("domain_id"), write.domain.df)

##
## Run bin/parse_domains.py
##

domain_cors <- read.table("data/domain_alns/PF06480.tab", header=T)
g1 <- subset(domain_cors, stable_id == "ENSP00000268704")$Omega
g2 <- subset(domain_cors, stable_id == "ENSP00000269143")$Omega

plot(g1, g2)
cor.test(g1, g2, use="complete.obs", method="spearman")

plot(subset(domain_cors, stable_id == "ENSP00000268704")$Omega, type="l")
lines(subset(domain_cors, stable_id == "ENSP00000269143")$Omega, col="red")
