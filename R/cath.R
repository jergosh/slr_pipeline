library(plyr)
library(seqinr)
library(phylosim)

setwd("~/Documents/projects/slr_pipeline/")

is.letter <- function(str) {
  is.upper_RE <- "[A-Za-z]"
  grepl(pattern = is.upper_RE, x = str)
}

is.upper <- function(str) {
  is.upper_RE <- "[A-Z]"
  grepl(pattern = is.upper_RE, x = str)
}

is.lower <- function(str) {
  is.upper_RE <- "[a-z]"
  grepl(pattern = is.upper_RE, x = str)
}

cath_all <- read.table("data/cath_new.tab", sep="\t", stringsAsFactors=F, header=T)
cath_raw <- cath_all
cath_map <- read.table("data/cath/cath_map.tab", sep="\t", stringsAsFactors=F, header=F)
colnames(cath_map) <- c("stable_id", "coords", "cath_id", "pdb")
head(cath_map)
sort(table(cath_map$cath_id))
sum(table(cath_map$cath_id) > 10)
length(unique(cath_map$cath_id))

id2dataset = dlply(cath_all, .variables="stable_id", function(df) { df$dataset[1] })

sum(!is.na(aln$Omega)) / nrow(aln)
unique(aln$stable_id)

# Also do the matching on the basis of shared alignments?
# for (stable_id in unique(slr  )

# FIXME This is actually filtering on the fraction of length of the gene, not the domain length/available data
# Could redo this to contain only domain regions
filter.cath <- function(df) {
  # if (sum(!is.na(df$omega)) / nrow(df) > 0.8) {
  if (sum(!is.na(df$omega)) > 50) {
    return(df)
  } else {
    return(data.frame())
  }
}

cath_filtered <- ddply(cath_all, .variables=c("stable_id"), filter.cath)
nrow(cath_filtered)
summary(cath_filtered$omega)
# write.table(cath_filtered, file="data/cath/cath_filtered.tab", sep="\t", quote=F)

## Get summary stats for domain families
## to get a lay of the land
family_sizes <- ddply(cath_filtered, c("cath_id"), function(df) {
  data.frame(N=length(unique(df$stable_id)))
})
hist(family_sizes$N, breaks=100)
subset(family_sizes, N > 100)

domain.cors <- function(df) {
  print(paste(df$cath_id[1], nrow(df)))
  if (length(unique(df$stable_id)) < 2) {
    return(data.frame())
  }
  cath_id <- df$cath_id[1]
  res <- data.frame(cath_id=character(), id_1=character(), id_2=character(), corr=numeric(), pval=numeric(), stringsAsFactors=F)
  for (comb in combn(as.character(unique(df$stable_id)), 2, simplify=F)) {
    if (comb[1] == comb[2]) {
      print("IDs identical!")
    }
    omega_1 <- subset(df, stable_id == comb[1])$omega
    omega_2 <- subset(df, stable_id == comb[2])$omega
    
    # Make sure there is actually stuff to correlate
    if (sum((!is.na(omega_1)) & (!is.na(omega_2))) < 10) {
      res[nrow(res)+1, ] <- c(cath_id, comb[1], comb[2], NA, NA)
    } else {
      corr <- cor.test(omega_1, omega_2, use="complete.obs", method="spearman")
      res[nrow(res)+1, ] <- c(cath_id, comb[1], comb[2], corr$estimate, corr$p.value)
    }
  }
  
  res
}

# cath.cors.spearman <- ddply(cath_filtered, .variables="cath_id", .fun=domain.cors)
cath.cors.spearman <- read.table("cath_cors_py_2.tab", sep="\t", header=T, stringsAsFactors=F)
head(cath.cors.spearman)
str(cath.cors.spearman)
hist(cath.cors.spearman$pval, breaks=100)
hist(-log10(cath.cors.spearman$pval), breaks=100)
hist(cath.cors.spearman$p.adjusted, breaks=100)
hist(cath.cors.spearman$n_obs, breaks=100)
# cath.cors <- cath.cors.spearman
sum(cath.cors.spearman$p.adjusted < 1e-5, na.rm=T)

check.ds <- function(r) {
  # print(c(id2dataset[[r$id_1[1]]], id2dataset[[r$id_2[1]]]))
  if (id2dataset[[ r$id_1[1] ]] == id2dataset[[ r$id_2[1] ]]) {
    return(data.frame())
  } else {
    return(r)
  }
}
nrow(cath.cors.spearman)
cath.cors.filtered <- ddply(cath.cors.spearman, .variables=c("id_1", "id_2"), check.ds)
nrow(cath.cors.filtered)
cath.cors <- cath.cors.filtered

cath.cors$pval <- as.numeric(cath.cors$pval)
cath.cors$cor <- as.numeric(cath.cors$cor)
cath.cors$p.adjusted <- p.adjust(cath.cors$pval, method="BH")
write.table(cath.cors, file="data/cath_cors.tab", sep="\t", quote=F)

hist(cath.cors$cor, breaks=100)
sum(cath.cors.filtered$cor > 0, na.rm=T)
neg.cor <- subset(cath.cors.filtered, cor < 0)
hist(neg.cor$pval, breaks=100)
sum(neg.cor$p.adjusted < 0.05)

cath.fam.sizes <- ddply(cath_filtered, c("cath_id"), function(df) {
  data.frame(N=length(unique(df$stable_id)))
})
cath.fam.sizes[order(cath.fam.sizes$N, decreasing=T)[21:30], ]

cath.corr.means <- ddply(cath.cors, .variables="cath_id", function(df) { data.frame(mean=mean(df$cor, na.rm=T), pval=sum(df$pval < 0.001),
                                                                                    pval_frac=(sum(df$pval < 0.001)/sqrt(nrow(df)))) })
hist(cath.corr.means$mean, breaks=100)
cath.corr.means[order(cath.corr.means$pval_frac, decreasing=T)[1:25], ]
hi_corr_fams <- subset(cath.corr.means, V1 > 0.2)$cath_id
nrow(subset(cath.cors, cath_id %in% hi_corr_fams))

# cath.cors.filtered <- cath.cors.filtered[!is.na(cath.cors.filtered)]
# cath.cors.filtered <- data.frame(matrix(unlist(cath.cors.filtered), ncol=6, byrow=T))
# colnames(cath.cors.filtered) <- colnames(cath.cors)
# hist(as.numeric(cath.cors.filtered$pval), breaks=100, freq=F)
pdf("domain_correlations.pdf", height=8, width=9)
hist(as.numeric(cath.cors.filtered$pval), breaks=100, cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5,
     main=expression(paste("Correlation of ", omega, " profiles")), xlab="P-value")
dev.off()
hist(as.numeric(cath.cors.filtered$corr), breaks=100)

cath.cors.filtered$corr <- as.numeric(cath.cors.filtered$corr)
cath.cors.filtered$p.adjusted <- p.adjust(as.numeric(cath.cors.filtered$pval), method="BH")
sum(cath.cors.filtered$p.adjusted < 0.05, na.rm=T)

cath.cors.filtered[order(cath.cors.filtered$corr, decreasing=T)[1:50], ]

cath.cors.means <- ddply(cath.cors.filtered, .variables=c("cath_id"), function(df) { data.frame(n=sum(complete.cases(df)), mean.cor=mean(df$cor, na.rm=T)) })
hist(cath.cors.means$mean.cor, breaks=100)
cath.cors.means[order(cath.cors.means$mean.cor)[1:10], ]

## Could also calc 'joint' (multiplied) p-values -- although not clear what to do with neg. correlations

ddply(cath_all, .variables="cath_id", .fun=function(df) { print(c(length(unique(df$dataset)), length(unique(df$stable_id))))})

# Looks like there's signal in there, should now focus on extracting interesting case studies

# Mean omega from cath_filtered?
mean_omega <- ddply(cath_filtered, .variables="cath_id", function(df) { data.frame(mean_omega=mean(df$omega, na.rm=T), n=length(unique(df$stable_id))) })
mean_omega_filtered <- subset(mean_omega, n > 10)
hist(mean_omega_filtered$mean_omega, breaks=100)
top_omega <- mean_omega_filtered[order(mean_omega_filtered$mean_omega), ]

pdb_test <- subset(cath_filtered, cath_id == "2.40.30.110")
ddply(pdb_test, .variables="stable_id", function(df) { mean(df$omega, na.rm=T) })
# Normalise?
pdb_test_norm <- ddply(pdb_test, .variables="stable_id", function(df) { data.frame(pos=df$pos, omega=df$omega/mean(df$omega, na.rm=T)) })
ggplot(pdb_test, aes(x=pos, y=omega, colour=stable_id)) + geom_line()
ggplot(pdb_test_norm, aes(x=pos, y=omega, colour=stable_id)) + geom_line()
ddply(pdb_test, .variables="stable_id", function(df) {nrow(df)})

##
## Instead of using mean correlation, we can count the number of sign. adj. p
##
sign_fraction <- function(df) {
  n <- nrow(df)
  n_sign <- sum(df$p.adjusted < 0.05, na.rm=T)
  data.frame(n=n, n_sign=n_sign, frac=n_sign/n)
}

cath.cors.frac <- ddply(cath.cors.filtered, .variables="cath_id", sign_fraction)
cath.cors.frac[order(cath.cors.frac$frac, decreasing=T)[1:50], ]
cath.cors.frac[order(cath.cors.frac$n_sign, decreasing=T)[1:50], ]

##
## Check how much 'alignment coverage' sites in interesting domains have
##

##
## Extract information to construct a domain alignment
##
pep_dir <- "data/ens/73/seqsets_pep/Eutheria"

get_domain <- function(st_id, c_id) {
  domain_ann <- subset(cath_map, stable_id == st_id & cath_id == c_id)
  coords <- as.numeric(strsplit(domain_ann$coords[1], ":")[[1]])
  slr_subset <- subset(slr_raw, stable_id == st_id)
  slr_subset <- slr_subset[coords[1]:coords[2], ]
  slr_subset$Omega
  # Need the sequence as well
  # seqs <- read.fasta(paste(c(pep_dir, "/", substr(ds, 1, 2), "/", ds, ".fa"), collapse=""))
  # seq <- seqs[[st_id]][coords[1]:coords[2]]
}

get_omega <- function(st_id, c_id, seq) {
  omegas <- rep(NA, length(seq))
  omegas[is.letter(seq)] <- get_domain(st_id, c_id)
  omegas[is.na(omegas)] <- 0.0
  omegas[!is.finite(omegas)] <- 0.0
  # omegas[is.lower(omegas)] <- NA
  
  omegas
}

domain_id <- "1.10.287.600"
seqs <- read.fasta(paste(c(domain_id, ".fa"), collapse=""), forceDNAtolower=F)

# ens_ids <- unique(subset(cath_filtered, cath_id == domain_id)$stable_id)
ens_ids <- unique(c(subset(cath.cors.filtered, cath_id == domain_id & corr > 0.7)$id_1,
                    subset(cath.cors.filtered, cath_id == domain_id & corr > 0.7)$id_2))

aln <- read.fasta(paste(c(domain_id, ".fa"), collapse=""),
                  forceDNAtolower=F, set.attributes=F)

aln <- t(apply(as.matrix(aln), 1, function(s) {unlist(s)}))
aln <- aln[ens_ids, ]
cols <- apply(aln, 2, function(x) { sum(x == "-" | is.lower(x)) != nrow(aln) })
aln <- aln[, cols]


aln_ds <- list()
tracks <- list()
for (ens_id in ens_ids) {
  if (id2dataset[[ens_id]] %in% aln_ds) {
    next
  }
  aln_ds[[ens_id]] = id2dataset[[ens_id]]
  
  omegas <- get_omega(ens_id, domain_id, seqs[[ens_id]])
  omegas <- omegas[cols]
  omegas <- omegas / mean(omegas[omegas != 0.0], na.rm=T)
  slr_track <- data.frame(pos=1:length(omegas),
                          score=pmin(omegas, 1.0),
                        y_hi=pmin(omegas, 1.0),
                        id=ens_id,
                        layout="below",
                        height=3,
                        color.gradient="white,red")

  tracks[[ens_id]] <- slr_track
}  
sim <- PhyloSim()
sim$.alignment <- aln[names(aln_ds), ]

pdf("helix_hairpin_bin.pdf", width=11, height=11)
plot.PhyloSim(sim,
              plot.chars=F,
              # plot.labels=T,
              # aln.char.text.size=1.8,
              num.pages=1,
              # layout.ancestors=F,
              # aln.to.tree.size=5,
              axis.text.size=8,
              # tree.xlim=xlims[[symbol]],
              plot.legend=T,
              color.scheme='protein',
              tracks=tracks,
              tracks.ylim=c(0, 1)
)
dev.off()

cath.cors.domain <- subset(cath.cors.filtered, cath_id == domain_id & corr > 0.7)
# cath.cors.domain <- cath.cors.domain[order(cath.cors.domain$corr, decreasing=T)[1:100], ]
unique(c(cath.cors.domain$id_1, cath.cors.domain$id_2))

ddply(cath.cors.domain, .variables=c("id_1", "id_2"), function(df) {
  id_1 <- df$id_1[1]
  id_2 <- df$id_2[1]
  omega_1 <- subset(cath_filtered, stable_id == id_1)$omega
  omega_1[omega_1 == 0] <- NA
  omega_2 <- subset(cath_filtered, stable_id == id_2)$omega
  omega_2[omega_2 == 0] <- NA
  print(cor.test(omega_1, omega_2, method="spearman"))
  
  plot(omega_1, omega_2, main=paste(c(id_1, id_2)), xlab=id_1, ylab=id_2, main=)
})

unique(subset(cath_filtered, cath_id == "1.20.1070.10")$pdb_id)

cath.cors.domain <- subset(cath.cors.filtered, cath_id == domain_id & id_1 == "ENSP00000248072" & id_2 == "ENSP0000444134")
