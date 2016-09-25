library(ggplot2)
library(GenomicRanges)
library(plyr)

head(slr_all)
setwd("~/Documents/projects/slr_pipeline")

slr_subs <- subset(slr_all, stable_id == "ENSP00000259951")
nrow(slr_subs)

color_file <- function(slr_data) {
  
  colors <- data.frame(r=rep(0.5, nrow(slr_data)), g=0.5, b=0.5)
  #  c(0.5, 0.5, 0.5)
  colors[slr_data$lower > 1.0, ] <- c(1.0, 0.0, 0.0)
  colors[slr_data$upper < 1.0, ] <- c(0.0, 0.0, 1.0)

  colors
}

slr_mhc <- subset(slr_all, stable_id == "ENSP00000379873")
head(slr_mhc)
sum(slr_mhc$lower > 1)
mean(slr_mhc$Omega)
nrow(slr_mhc)

slr_mhc$human_idx <- 1:nrow(slr_mhc)
slr_mhc$xmin <- slr_mhc$human_idx - 0.5
slr_mhc$xmax <- slr_mhc$human_idx + 0.5
slr_mhc$ymin <- slr_mhc$lower
slr_mhc$ymax <- slr_mhc$upper

slr_mhc$colour <- "grey"
slr_mhc$colour[slr_mhc$lower > 1.0] <- "red"
slr_mhc$colour[slr_mhc$upper < 1.0] <- "blue"
slr_mhc$colour <- factor(slr_mhc$colour)

pdf("hla_talk.pdf", height=5, width=11)
ggplot(slr_mhc, aes(x=human_idx, y=Omega, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=colour)) +
  # scale_colour_identity(labels=levels(slr_mhc$colour)) +
  scale_fill_manual(values=c(grey="grey", red="red", blue="blue")) +
  geom_hline(yintercept=1.0) +
  coord_cartesian(xlim=c(0.5, nrow(slr_mhc)+0.5), ylim=c(0, 3)) +
  geom_rect() +
  geom_point(shape=15, size=1) +
  theme_bw(base_size=18) +
  theme(legend.position="none") +
  ylab(expression(omega)) +
  xlab("Position")
dev.off()

pdf("hla_talk_short.pdf", height=5, width=11)
ggplot(slr_mhc, aes(x=human_idx, y=Omega, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=colour)) +
  # scale_colour_identity(labels=levels(slr_mhc$colour)) +
  scale_fill_manual(values=c(grey="grey", red="red", blue="blue")) +
  geom_hline(yintercept=1.0) +
  coord_cartesian(xlim=c(0.5, nrow(slr_mhc)+0.5), ylim=c(0, 3)) +
  geom_rect() +
  geom_point(shape=15, size=1) +
  theme_bw(base_size=18) +
  theme(legend.position="none") +
  ylab(expression(omega)) +
  xlab("Position")
dev.off()


write.table(color_file(slr_mhc), file="mhc_color_file.tab", sep="\t", quote=FALSE, col.names=FALSE)

# Barplots of fraction of sites under positive selection
# This should be integrated into the main file with structural analysis
pdb_master <- subset(pdb_master, complete.cases(pdb_master))
head(pdb_master)
pdb_master$buried <- TRUE
pdb_master$buried[pdb_master$rsa > 0.25] <- FALSE
pdb_master$positive <- pdb_master$omega > 1

# Solvent accessibility
frac.buried <- sum(pdb_master$buried & pdb_master$omega > 1.0) / sum(pdb_master$buried)
frac.exposed <- sum((!pdb_master$buried) & pdb_master$omega > 1.0) / sum(!pdb_master$buried)

rsa.df <- data.frame(RSA=factor(c("buried", "exposed"), levels=c("buried", "exposed")),
                     Fraction=c(frac.buried, frac.exposed))

ggplot(rsa.df, aes(x=RSA, y=Fraction, fill=RSA)) +
  scale_fill_manual(values=c("indianred", "darkred")) +
  geom_bar(stat="identity") +
  ggtitle("Positive selection and solvent accessibility") +
  theme_bw(base_size=18) +
  theme(legend.position="none")

# Secondary structure
pdb_master$sec_simple <- FALSE
pdb_master$sec_simple[pdb_master$secondary %in% c("G", "H", "I")] <- "H"
pdb_master$sec_simple[pdb_master$secondary %in% c("B", "E")] <- "B"
pdb_master$sec_simple[pdb_master$secondary %in% c("S", "T", "-")] <- "L"
sum(pdb_master$sec_simple == FALSE)

frac.helix <- sum((pdb_master$sec_simple == "H") & pdb_master$omega > 1.0) / sum(pdb_master$sec_simple == "H")
frac.sheet <- sum((pdb_master$sec_simple == "B") & pdb_master$omega > 1.0) / sum(pdb_master$sec_simple == "B")
frac.loop <- sum((pdb_master$sec_simple == "L") & pdb_master$omega > 1.0) / sum(pdb_master$sec_simple == "L")

secondary.df <- data.frame(secondary=factor(c("Helix", "Beta sheet", "Other"),
                             levels=c("Helix", "Beta sheet", "Other")),
                     Fraction=c(frac.helix, frac.sheet, frac.loop))

pdf("secondary_barplot.pdf")
ggplot(secondary.df, aes(x=secondary, y=Fraction, fill=secondary)) +
  scale_fill_manual(values=c("lightblue", "blue", "darkblue")) +
  geom_bar(stat="identity") +
  ggtitle("Positive selection and secondary structure") +
  theme_bw(base_size=18) +
  theme(legend.position="none")
dev.off()

# Disorder predictions?

# SLR and domains
feature_data_slr <- subset(feature_data, stable_id %in% unique(slr_all$stable_id))
feature_data_slr <- subset(feature_data_slr, type == "pfam")

pfam.ranges.slr <- with(feature_data_slr,
                        GRanges(seqnames=stable_id, ranges=IRanges(start=start, end=end)))


slr_all <- ddply(slr_all, "stable_id", function(df) { df$human_idx <- 1:nrow(df); df })
slr.all.ranges <- with(slr_all, GRanges(seqnames=stable_id, ranges=IRanges(start=human_idx, end=human_idx)))
slr.all.ranges <- reduce(slr.all.ranges)
sum(width(slr.all.ranges))

# TODO Should figure out the overlaps between where
# SLR results are available and domains
pfam.ranges.reduced <- reduce(pfam.ranges.slr)
## sum(width(pfam.ranges.reduced))
## domain.slr.overlap <- intersect(slr.all.ranges, pfam.ranges.slr)
## sum(width(domain.slr.overlap))

thr = 0.05
slr_sign = subset(slr_all, lower > 1)

slr.ranges <- GRanges(seqnames=slr_sign$stable_id, ranges=IRanges(start=slr_sign$human_idx, end=slr_sign$human_idx))
               
slr.domain.overlap <- intersect(slr.ranges, pfam.ranges.reduced)
domain.overlap.observed <- sum(width(slr.domain.overlap)) / sum(width(slr.ranges))
domain.overlap.expected <- sum(width(pfam.ranges.reduced)) / nrow(slr_all)

# domain.overlap.observed <- sum(width(slr.domain.overlap)) / nrow(slr_all)
# domain.overlap.expected <- sum(width(pfam.ranges.reduced)) / nrow(slr_all) * sum(width(slr.ranges)) / nrow(slr_all)

slr.domain.df <- data.frame(Overlap=factor(c("Expected", "Observed"),
                             levels=c("Expected", "Observed")),
                     Fraction=c(domain.overlap.expected, domain.overlap.observed))

pdf("slr_domain_barplot.pdf")
ggplot(slr.domain.df, aes(x=Overlap, y=Fraction, fill=Overlap)) +
  scale_fill_manual(values=c("skyblue2", "skyblue4")) +
  geom_bar(stat="identity") +
  ggtitle("Positive selection and domain annotations") +
  theme_bw(base_size=18) +
  theme(legend.position="none")
dev.off()

snp.domain.overlap <- intersect(snp.ranges, pfam.ranges.reduced)
snp.domain.overlap.observed <- sum(width(snp.domain.overlap)) / sum(width(snp.ranges))
# snp.domain.overlap.expected <- sum(width(pfam.ranges.reduced)) / sum(seqlengths(pfam.ranges))
snp.domain.overlap.expected

snp.domain.df <- data.frame(Overlap=factor(c("Expected", "Observed"),
                             levels=c("Expected", "Observed")),
                     Fraction=c(snp.domain.overlap.expected, snp.domain.overlap.observed))

pdf("snp_domain_barplot.pdf")
ggplot(snp.domain.df, aes(x=Overlap, y=Fraction, fill=Overlap)) +
  scale_fill_manual(values=c("skyblue1", "skyblue3")) +
  geom_bar(stat="identity") +
  ggtitle("psSNPs and domain annotations") +
  theme_bw(base_size=18) +
  theme(legend.position="none")
dev.off()

# PTMs
# psites, acek, ubi
ptm_data <- read.table("data/PTMfunc/PTM_lists/acek_spectral_counts.tab",
                       header=TRUE, sep="\t")


head(ptm_data)
ptm_data_subset <- subset(ptm_data, ProteinID %in% unique(slr_all$stable_id))
nrow(ptm_data_subset)

slr_all_seqlengths <- ddply(slr_all, "stable_id", function(df) max(df$human_idx))
slr_seqlengths_map <- slr_all_seqlengths$V1
names(slr_seqlengths_map) <- slr_all_seqlengths$stable_id
head(slr_seqlengths_map)
  
ptm_data_subset$seqlength <- slr_seqlengths_map[as.character(ptm_data_subset$ProteinID)]
ptm_data_subset <- subset(ptm_data_subset, complete.cases(ptm_data_subset))
head(ptm_data_subset)

ptm.ranges <- with(ptm_data_subset, GRanges(seqnames=ProteinID,
                                            ranges=IRanges(start=pmax(1, PTMpos-7),
                                              end=pmin(ptm_data_subset$seqlength, PTMpos+7)),
                                            seqlengths=slr_seqlengths_map))
head(ptm.ranges)

ptm.ranges.reduced <- reduce(ptm.ranges)
ptm.slr.overlap <- intersect(slr.ranges, ptm.ranges.reduced)
ptm.overlap.observed <- sum(width(ptm.slr.overlap)) / sum(width(slr.ranges))
ptm.overlap.expected <- sum(width(ptm.ranges.reduced)) / nrow(slr_all)

ptm.slr.df <- data.frame(Overlap=factor(c("Expected", "Observed"),
                             levels=c("Expected", "Observed")),
                     Fraction=c(ptm.overlap.expected, ptm.overlap.observed))


pdf("slr_acek_barplot.pdf", height=10, width=7)
ggplot(ptm.slr.df, aes(x=Overlap, y=Fraction, fill=Overlap)) +
  scale_fill_manual(values=c("skyblue2", "skyblue4")) +
  geom_bar(stat="identity") +
  ggtitle("Positive selection and acetylation sites") +
  theme_bw(base_size=18) +
  theme(legend.position="none")
dev.off()


# Disorder and SLR
slr.ids <- getBM(attributes = c("peptide", "ensembl_transcript_id", "ensembl_peptide_id", 
                                     "transcript_biotype"), filters="ensembl_peptide_id",
                 values=unique(slr_all$stable_id), mart=ensembl.db)

write.table(snp.ann.unique, file="data/snp_ann.txt", quote=F, sep="\t")
