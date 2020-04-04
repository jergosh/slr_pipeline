library(biomaRt)
library(plyr)
library(ggplot2)

# uni <- useMart("unimart", "uniprot", verbose=T)
# ec_numbers <- getBM(attributes = c("accession", "ec_number"), filters =
                               "accession", values = unique(pdb_master$uniprot_id), mart = uni)

# Catalytic Site Atlas
csa <- read.table("~/Documents/projects/slr_pipeline/data/CSA_2_0_121113.txt", sep=",", stringsAsFactors=FALSE, header=TRUE)
nrow(csa)
head(csa)
csa <- subset(csa, PDB.ID %in% unique(pdb_master$pdb_id))
length(unique(subset(csa, PDB.ID %in% unique(pdb_master$pdb_id))$PDB.ID))

pdb_master_csa <- subset(slr_structure_all, pdb_id %in% unique(csa$PDB.ID))
nrow(pdb_master_csa)

csa_indices <- rep(FALSE, nrow(pdb_master_csa))

ddply(csa, c("PDB.ID", "CHAIN.ID", "RESIDUE.NUMBER"), function(df) {
  csa_indices[pdb_master_csa$pdb_id == df$PDB.ID[1] &
                pdb_master_csa$pdb_pos == df$RESIDUE.NUMBER[1] &
                pdb_master_csa$pdb_chain == df$CHAIN.ID[1]] <<- TRUE
})
sum(csa_indices)

pdb_master_csa$csa <- factor("Other", levels=c("Catalytic site", "Other"))
pdb_master_csa$csa[csa_indices] <- "Catalytic site"

fractions.csa <- make.fractions(pdb_master_csa, c("buried", "csa"))
ggplot(fractions.csa, aes(x=buried, y=Fraction, fill=csa)) +
  scale_fill_manual(values=c("lightblue", "blue", "darkblue",
                             "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
                             "brown1", "brown3", "brown4")) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Positive selection and structure") +
  theme_bw(base_size=18)
ggsave(paste(results, "catalytic_sites.pdf", sep="/"), height=8, width=11)


csa_pos <- subset(pdb_master_csa, omega > 1 & csa == TRUE)
csa_pos$human_idx <- csa_pos$ens_pos

csa_pos_merged <- join(csa_pos, slr_all, by=c("stable_id", "human_idx"))
subset(csa_pos_merged, Adj.Pval < thr)

# Density plots of omega ----
ggplot(pdb_master_csa, aes(x=omega, linetype=buried, colour=sec_simple)) +
  geom_density()

# Cropped
ggplot(pdb_master_csa, aes(x=omega, linetype=buried, colour=sec_simple)) +
  coord_cartesian(xlim=c(0, 5), ylim=c(0, 10)) +
  geom_density()

# We really have a single site from the glutathione transferase.


# Breakdown plot ----
plot.csa.fractions <- function(fractions, main) {
  ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=is_domain)) +
  #   scale_fill_manual(values=c("lightblue", "blue", "darkblue",
  #                              "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
  #                              "brown1", "brown3", "brown4")) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle(main) +
  theme_bw(base_size=18)
}

"Positive selection in/out of domains in enzymes"
fractions.csa <- make.fractions(pdb_master_csa, c("sec_simple", "is_domain"))
plot.csa.fractions(fractions.csa, "Positive selection and domains (enzymes)")
fractions.csa.buried <- make.fractions(subset(pdb_master_csa, buried=="buried"), c("sec_simple", "is_domain"))
plot.csa.fractions(fractions.csa.buried, "Positive selection and domains (enzymes)")
fractions.csa.exposed <- make.fractions(subset(pdb_master_csa, buried=="exposed"), c("sec_simple", "is_domain"))
plot.csa.fractions(fractions.csa.exposed, "Positive selection and domains (enzymes)")

