bs_sites <- read.table("data/bs_sites.tab", header=F, sep="\t")
head(bs_sites)
colnames(bs_sites) <- c("dataset", "sample", "lnL", "pos", "res", "C1", "C2", "C3", "C4", "cat")

ddply(bs_sites, "dataset", function(df) {
  data.frame(N=length(unique(df$sample)))
})

slr_sites <- subset(slr_all, dataset == "7631_1")[, c("Site", "human_idx")]
rownames(slr_sites) <- slr_sites$Site
bs_subset <- subset(bs_sites, dataset == "7631_1" & pos %in% slr_sites$Site)
bs_subset$idx <- slr_sites[as.character(bs_subset$pos), "human_idx"]
bs_subset <- bs_subset[order(bs_subset$pos), ]
