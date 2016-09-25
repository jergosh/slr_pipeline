library(ggplot2)
library(GenomicRanges)
library(plyr)

# The overall SLR results


# Preparing the data frames with structural mapping results -------------------

head(slr_all)
setwd("~/Documents/projects/slr_pipeline")
# results <- "results/2015-04-08"
results <- "results/2016-09"

# Barplots of fraction of sites under positive selection
# This should be integrated into the main file with structural analysis
pdb_master <- read.table("data/pdb_master_table_DSSP_200916.tab", header=F,
                         #colClasses=c("character", "character", "numeric", "factor", "numeric", "numeric"),
                         na.strings=c("None", "NA"), stringsAsFactors=F)


colnames(pdb_master) <- c("stable_id", "ens_pos", "uniprot_id", "uniprot_pos", "pdb_id", "pdb_chain", "pdb_pos", "secondary", "rsa", "omega")
pdb_master$ens_pos <- pdb_master$ens_pos + 1


nrow(pdb_master)
sum(complete.cases(pdb_master))
# pdb_master <- subset(pdb_master, complete.cases(pdb_master))
head(pdb_master)

mapped_lengths <- ddply(pdb_master, c("stable_id", "pdb_id", "pdb_chain"), function(df) {
  nrow(df)
})

hist(mapped_lengths$V1, breaks=200)
hist(subset(mapped_lengths, V1 < 50)$V1)
sum(mapped_lengths$V1 >= 50)

pdb_master$buried <- factor("buried", levels=c("buried", "exposed"))
pdb_master$buried[pdb_master$rsa > 0.25] <- "exposed"
pdb_master$positive <- pdb_master$omega > 1
head(pdb_master)



# This presumably needs to be re-downloaded using the Perl script I wrote
write.table(unique(pdb_master$stable_id), file="data/all_ids.tab", quote=F, row.names=F, col.names=F)

feature_data <- read.table("~/Documents/projects/slr_pipeline/data/protein_features.tab",
                           header=F, fill=T, sep="\t", stringsAsFactors=F, quote="")
colnames(feature_data) <- c("stable_id", "start", "end", "type", "intepro_id", "domain_id", "domain_desc")
head(feature_data)
tmhmm_features <- subset(feature_data, type == "tmhmm")
head(tmhmm_features)
all_tm <- unique(tmhmm_features$stable_id)
pdb_master$tm <- pdb_master$stable_id %in% all_tm
summary(pdb_master$tm)
length(unique(subset(pdb_master, tm == TRUE)$stable_id))
tmhmm_features <- subset(tmhmm_features, stable_id %in% pdb_master$stable_id)

pdb_master$tm_sites <- FALSE

for (i in 1:nrow(tmhmm_features)) {
  stable_id <- tmhmm_features[i, "stable_id"]
  start <- tmhmm_features[i, "start"]
  end <- tmhmm_features[i, "end"]

  pdb_master[pdb_master$stable_id == stable_id &
               pdb_master$ens_pos >= start &
               pdb_master$ens_pos <= end, "tm_sites"] <- TRUE
}
sum(pdb_master$tm_sites)

## TODO look at the distribution of lengths of the TM regions
tm_lengths <- tmhmm_features$end-tmhmm_features$start

# Transmembrane proteins are those that have at least one TM region annotated
pdb_master$tm_full <- "Globular"
pdb_master$tm_full[pdb_master$tm == TRUE & pdb_master$tm_sites == TRUE] <- "Transmembrane (sites)"
pdb_master$tm_full[pdb_master$tm == TRUE & pdb_master$tm_sites == FALSE] <- "Transmembrane (protein)"
pdb_master$tm_full <- factor(pdb_master$tm_full, levels=c("Transmembrane (sites)", "Transmembrane (protein)", "Globular"))

# How many sites of each kind
summary(pdb_master$tm_full)


# Secondary structure
pdb_master$sec_simple <- FALSE
# pdb_master$sec_simple[pdb_master$secondary %in% c("G", "H", "I")] <- "Helix"
# pdb_master$sec_simple[pdb_master$secondary %in% c("B", "E")] <- "Beta sheet"
# pdb_master$sec_simple[pdb_master$secondary %in% c("S", "T", "-")] <- "Loop"

pdb_master$sec_simple[pdb_master$secondary %in% c("H")] <- "Helix"
pdb_master$sec_simple[pdb_master$secondary %in% c("E")] <- "Beta sheet"
pdb_master$sec_simple[pdb_master$secondary %in% c("B", "G", "I", "-")] <- "Coil"
pdb_master$sec_simple[pdb_master$secondary %in% c("S", "T")] <- "Turn"

sum(pdb_master$sec_simple == FALSE)
pdb_master$sec_simple <- factor(pdb_master$sec_simple, levels=c("Helix", "Beta sheet", "Coil", "Turn"))
pdb_master$structure <- factor(paste(pdb_master$tm_full, '-', pdb_master$sec_simple),
                                  levels=paste(rep(levels(pdb_master$tm_full), 4), '-', rep(levels(pdb_master$sec_simple), each=3)))

# Enzyme (EC) classification
write.table(unique(pdb_master$uniprot_id), file="~/Downloads/for_uniprot.tab", quote=F, row.names=F)

ec_numbers <- read.table("data/ec_uniprot.tab", sep="\t", stringsAsFactors=F, header=T)
ec_numbers$uniprot_id <- ec_numbers$Entry
ec_numbers <- subset(ec_numbers, EC.number != "")

pdb_master$is_enzyme <- factor("non-enzyme", levels=c("enzyme", "non-enzyme"))
pdb_master$is_enzyme[pdb_master$uniprot_id %in% ec_numbers$uniprot_id] <- "enzyme"
summary(pdb_master$is_enzyme)

pdb_master_globular <- subset(pdb_master, tm_full == "Globular")
# Mapping between the 'full' results and the structural coverage
pdb_master_globular$human_idx <- pdb_master_globular$ens_pos

slr_structure <- join(pdb_master_globular, slr_raw, by=c("stable_id", "human_idx"), type="inner")
slr_structure_raw <- slr_structure
nrow(subset(slr_structure, omega != Omega)[, c("omega", "Omega")])

# FIXME This seems wrong! Have I f u?
# slr_structure_pos <- slr_structure
# slr_structure_pos$Pval[slr_structure_pos$omega < 1] <- 1 - slr_structure_pos$Pval[slr_structure_pos$omega < 1]
# slr_structure_pos$Adj.Pval[slr_structure_pos$omega < 1] <- 1

# Join the PDB-mapped sites to the SLR results
# This is the strict one
pdb_master$human_idx <- pdb_master$ens_pos
slr_structure_all <- join(pdb_master, slr_raw, by=c("stable_id", "human_idx"), type="inner")
slr_structure_raw <- slr_structure_all
nrow(subset(slr_structure_all, omega != Omega)[, c("omega", "Omega")])

# This needs to be up to date
tree_stats <- read.csv("tree_stats_new.csv", header=T, stringsAsFactors=F)
tree_stats$ortholog_frac <- tree_stats$X..of.Paralogs / tree_stats$X..of.Leaves
datasets <- unlist(lapply(strsplit(tree_stats$Directory, "/", fixed=T), function(x) { x[length(x)] }))
tree_stats$dataset <- unlist(lapply(strsplit(datasets, ".", fixed=T), function(x) { x[1] }))

length(unique(slr_structure_raw$dataset))
sum(unique(slr_structure_raw$dataset) %in% unique(tree_stats$dataset)) # Are 

slr_structure_both <- slr_structure_raw
slr_structure_both$Adj.Pval <- p.adjust(slr_structure_raw$Pval, method="BH")
head(slr_structure_both)
slr_structure_both_exposed <- subset(slr_structure_both, rsa < 0.2)
nrow(slr_structure_both_exposed)
sum(slr_structure_both_exposed$Adj.Pval < 0.05)/nrow(slr_structure_both_exposed)

tree_stats_filtered <- subset(tree_stats, X..of.Human.seqs == 1 & ortholog_frac <= 0.1)

# This is the strict version
str_subset_long <- subset(mapped_lengths, V1 >= 100)
slr_structure_filtered <- join(slr_structure_raw, str_subset_long, type="inner")
slr_structure_all_strict <- slr_structure_filtered
# slr_structure_all_strict <- subset(slr_structure_all_strict, dataset %in% tree_stats_filtered$dataset)
length(unique(slr_structure_all_strict$stable_id))
 
# slr_structure_all_strict$Pval[slr_structure_all_strict$omega < 1.0] <- runif(sum(slr_structure_all_strict$omega < 1.0))

# Convert from a two-sided to one-sided test
slr_structure_all_strict$Pval[slr_structure_all_strict$omega < 1.0] <- 1 - slr_structure_all_strict$Pval[slr_structure_all_strict$omega < 1.0]/2
slr_structure_all_strict$Pval[slr_structure_all_strict$omega > 1.0] <- slr_structure_all_strict$Pval[slr_structure_all_strict$omega > 1.0]/2

hist(slr_structure_all_strict$Pval, breaks=100)

slr_structure_all_strict$Adj.Pval <- p.adjust(slr_structure_all_strict$Pval, "BH")
sum(slr_structure_all_strict$Adj.Pval < 0.2 & slr_structure_all_strict$omega > 1)

hist(slr_structure_all_strict$Pval, breaks=100)
write.table(subset(slr_structure_all_strict, ), file="../cluster/data/slr_structure_strict.tab", row.names=F, quote=F, sep="\t")

# This is the lax version
slr_structure_all <- subset(slr_structure_all, dataset %in% tree_stats_filtered$dataset)
length(unique(slr_structure_all$stable_id))

slr_structure_all$Adj.Pval <- p.adjust(slr_structure_all$Pval, "BH")
slr_structure_all$Pval[slr_structure_all$omega < 1] <- 1
slr_structure_all$Adj.Pval[slr_structure_all$omega < 1] <- 1
sum(slr_structure_all$Adj.Pval < 0.2 & slr_structure_all$omega > 1)
sum(slr_structure_all$Adj.Pval < 0.05 & slr_structure_all$omega > 1)
hist(slr_structure$Pval, breaks=100)
write.table(slr_structure_all, file="../cluster/data/slr_structure.tab", row.names=F, quote=F, sep="\t")

