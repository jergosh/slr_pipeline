library(plyr)

setwd("~/Documents/projects/slr_pipeline/")

slr_div<- read.table("data/slr_div.tab", sep="\t", header=T)
head(slr_div)
slr_div_1 <- subset(slr_div, grepl("_1$", dataset))
slr_div_2 <- subset(slr_div, grepl("_2$", dataset))

setdiff(unique(slr_div_1$stable_id), unique(slr_div_2$stable_id))
setdiff(unique(slr_div_2$stable_id), unique(slr_div_1$stable_id))

common_ids <- intersect(unique(slr_div_1$stable_id), unique(slr_div_2$stable_id))

omegas_1 <- split(slr_div_1, slr_div_1$stable_id)
omegas_2 <- split(slr_div_2, slr_div_2$stable_id)

sub.cors <- adply(common_ids, 1, function(id) {
  test <- cor.test(omegas_1[[id]]$Omega, subset(omegas_2[[id]], stable_id == id)$Omega, method="spearman")
  data.frame(p=test$p.value, cor=test$estimate)
})

pdf("subsample_cors.pdf")
hist(-log10(sub.cors$p), breaks=100)
hist(sub.cors$cor, breaks=100)
dev.off()
# Paste _1 and _2 after 

