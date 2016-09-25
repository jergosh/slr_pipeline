residue_depths <- read.table("data/residue_depths.tab", sep='\t', header=F)
colnames(residue_depths) <- c("stable_id", "pdb_id", "pdb_chain", "pdb_pos", "rd1", "rd2")
head(residue_depths)

pdb_master_rd <- join(pdb_master, residue_depths, type="right")
nrow(pdb_master_rd)
hist(pdb_master_rd$rd1, breaks=100, xlim=c(0, 30))

pdb_master_low_rd <- subset(pdb_master_rd, rd1 < 4)
nrow(pdb_master_low_rd)
sum(pdb_master_low_rd$omega > 1)/nrow(pdb_master_low_rd)

pdb_master_mid_rd <- subset(pdb_master_rd, rd1 > 4 & rd1 < 7)
nrow(pdb_master_mid_rd)
sum(pdb_master_mid_rd$omega > 1)/nrow(pdb_master_mid_rd)

pdb_master_high_rd <- subset(pdb_master_rd, rd1 > 5)
nrow(pdb_master_high_rd)
sum(pdb_master_high_rd$omega > 1)/nrow(pdb_master_high_rd)

buried_sites <- subset(pdb_master_high_rd, omega > 1)
buried_sites$human_idx <- buried_sites$ens_pos
buried_sites_full <- join(buried_sites, slr_all)
buried_sites_full[, c("stable_id", "ens_pos", "pdb_id", "pdb_chain", "pdb_pos", "lower", "Adj.Pval")]
# print(paste(unique(pdb_master_high_rd$stable_id), collapse=" "))

plot(pdb_master_rd$rd2, pdb_master_rd$rsa)
