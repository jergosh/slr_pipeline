library(plyr)
results_dir <- "results/2015-05-30"

head(pdb_master_globular)

pdb_master_pos <- subset(pdb_master_globular, omega > 1.0)

nbhood_indices <- rep(FALSE, nrow(pdb_master_globular))

ddply(pdb_master_pos, c("stable_id", "ens_pos"), function(df) {
  nbhood_indices[pdb_master_globular$stable_id == df$stable_id[1] & (
                   pdb_master_globular$ens_pos == (df$ens_pos-1) |
                   pdb_master_globular$ens_pos == (df$ens_pos+1))] <<- TRUE
})
sum(nbhood_indices)

pdb_master_globular$nbhood <- nbhood_indices
pdb_master_pos_nbhood <- subset(pdb_master_globular, nbhood_indices)

mean(subset(pdb_master_globular, buried == "buried")$omega)
mean(subset(pdb_master_globular, buried == "exposed")$omega)

mean(subset(pdb_master_pos_nbhood, buried == "buried")$omega)
mean(subset(pdb_master_pos_nbhood, buried == "exposed")$omega)

ggplot(pdb_master_globular, aes(x=omega, colour=nbhood, linetype=buried)) +
  geom_density()

ggplot(pdb_master_globular, aes(x=rsa, colour=nbhood)) +
  geom_density()

# Need to remember to exclude the sites which we used to define the neighbourhoods
# in the first places

# TODO make a plot of the distribution of omega depending on buriedness
# Once again it's possible that the 

##
## Structural neighbourhood
##
write.table(subset(pdb_master_globular, omega > 1), file="data/pdb_master_table_high_omega.tab",
            sep="\t", quote=F)
cols_nbhood <- c("stable_id", "ens_pos", "uniprot_id", "uniprot_pos", "pdb_id", "pdb_chain", "pdb_pos", "secondary", "rsa", "omega")

high_omega <- read.table("data/pdb_master_table_high_omega.tab", stringsAsFactors=F,
                         sep="\t", header=FALSE)
colnames(high_omega) <- cols_nbhood

high_omega_nbhood_all <- data.frame()
for (thr in c(4.5, 6.5, 8.5)) {
  high_omega_nbhood <- read.table(paste0("data/pdb_master_high_omega_nbhood_", thr, ".tab"),
                                  stringsAsFactors=F, sep="\t", header=FALSE)
  colnames(high_omega_nbhood) <- c("stable_id", "pdb_id", "pdb_chain", "pdb_pos")
  high_omega_merged <- join(high_omega_nbhood, pdb_master_globular, by=c("stable_id", "pdb_id", "pdb_chain", "pdb_pos"))
  mean(subset(high_omega_merged, complete.cases(high_omega_merged))$omega)
  
  high_omega_merged$threshold <- thr
  high_omega_nbhood_all <- rbind(high_omega_nbhood_all, high_omega_merged)
}

high_omega_nbhood_all$thr <- factor(high_omega_nbhood_all$thr, levels=c("4.5", "6.5", "8.5"))
ggplot(subset(high_omega_nbhood_all, omega < 4), aes(x=omega, colour=factor(thr))) + 
  stat_ecdf()
ggsave(paste(results_dir, "nbhood_ecdf.pdf", sep="/"), height=7, width=9)

ggplot(subset(high_omega_nbhood_all, omega < 4), aes(y=omega, x=thr)) + 
  geom_boxplot()
ggsave(paste(results_dir, "nbhood_boxplot.pdf", sep="/"), height=7, width=9)

ggplot(subset(high_omega_nbhood_all, omega < 4), aes(x=thr, y=omega)) + 
  geom_violin()
ggsave(paste(results_dir, "nbhood_violin.pdf", sep="/"), height=7, width=9)


  

##
# Neighbourhood of catalytic sites?
# TODO Dump the file with catalytic sites for finding the neighbourhood
pdb_master_csa_nbhood <- subset(pdb_master_csa, csa == "Catalytic site")[, cols_nbhood]
write.table(pdb_master_csa_nbhood, file="data/pdb_master_csa.tab", sep="\t", row.names=F, quote=F)

pdb_master_csa_nbhood <- read.table("data/pdb_master_csa_nbhood.tab", stringsAsFactors=F, sep="\t", header=F)
colnames(pdb_master_csa_nbhood) <- c("stable_id", "pdb_id", "pdb_chain", "pdb_pos")
nrow(pdb_master_csa_nbhood)
pdb_master_csa_nbhood_merged <- join(pdb_master_csa_nbhood, pdb_master_globular, by=c("stable_id", "pdb_id", "pdb_chain", "pdb_pos"), type="inner")
plot(density(pdb_master_csa_nbhood_merged$omega))
hist(pdb_master_csa_nbhood_merged$omega, breaks=1000)
mean(pdb_master_csa_nbhood_merged$omega)
mean(pdb_master_csa$omega)




high_omega_nbhood_all <- data.frame()
for (thr in c(4.5, 6.5, 8.5)) {
  high_omega_nbhood <- read.table(paste0("data/pdb_master_csa_nbhood_", thr, ".tab"),
                                  stringsAsFactors=F, sep="\t", header=FALSE)
  colnames(high_omega_nbhood) <- c("stable_id", "pdb_id", "pdb_chain", "pdb_pos")
  high_omega_merged <- join(high_omega_nbhood, pdb_master_globular, by=c("stable_id", "pdb_id", "pdb_chain", "pdb_pos"))
  mean(subset(high_omega_merged, complete.cases(high_omega_merged))$omega)
  
  high_omega_merged$threshold <- thr
  high_omega_nbhood_all <- rbind(high_omega_nbhood_all, high_omega_merged)
}

high_omega_nbhood_all$thr <- factor(high_omega_nbhood_all$thr, levels=c("4.5", "6.5", "8.5"))
ggplot(subset(high_omega_nbhood_all, omega < 4), aes(x=omega, colour=factor(thr))) + 
  stat_ecdf()
ggsave(paste(results_dir, "csa_nbhood_ecdf.pdf", sep="/"), height=7, width=9)

ggplot(subset(high_omega_nbhood_all, omega < 4), aes(y=omega, x=thr)) + 
  geom_boxplot()
ggsave(paste(results_dir, "csa_nbhood_boxplot.pdf", sep="/"), height=7, width=9)

ggplot(subset(high_omega_nbhood_all, omega < 4), aes(x=thr, y=omega)) + 
  geom_violin()
ggsave(paste(results_dir, "csa_nbhood_violin.pdf", sep="/"), height=7, width=9)

##
## Global test
## 
## Can we detect enrichment of 
##
nbhood_45 <- subset(high_omega_nbhood_all, thr == "4.5")
nbhood_hyperg <- matrix(c(sum(nbhood_45$omega > 1, na.rm=T),
         sum(pdb_master_globular$omega > 1, na.rm=T)-sum(nbhood_45$omega > 1, na.rm=T),
         sum(nbhood_45$omega <= 1, na.rm=T),
         sum(pdb_master_globular$omega <= 1, na.rm=T)-sum(nbhood_45$omega <= 1, na.rm=T)
  ), nrow=2)
  
fisher.test(nbhood_hyperg)
