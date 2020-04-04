library(Hmisc)

library(grid)
library(gtable)
library(cowplot)

library(stringr)

setwd("~/Documents/projects/slr_pipeline/")

sd_binom <- function(n, p) { sqrt(p * (1 - p)/(n - 1)) }

make.fractions <- function(master.df, variables, thr) {
  ddply(master.df, .variables=variables, function(df) {
    sum_pos <- sum(df$omega > 1.0 & df$Adj.Pval < thr)
    fraction <- sum_pos / nrow(df)
    data.frame(number=nrow(df), pos=sum_pos, Fraction=fraction)
  })
}

make.fractions.neg <- function(master.df, variables, thr) {
  ddply(master.df, .variables=variables, function(df) {
    sum_neg <- sum(df$omega < 1.0 & df$Adj.Pval < thr)
    fraction <- sum_neg / nrow(df)
    data.frame(number=nrow(df), neg=sum_neg, Fraction=fraction)
  })
}

slr_structure_globular <- subset(slr_structure_all, tm == FALSE)
nrow(slr_structure_globular)
# make.fractions.lower <- function(master.df, variables) {
#   ddply(master.df, .variables=variables, function(df) {
#     fraction <- sum(df$lower > 1.0 & df$Adj.Pval < 0.2) / nrow(df)
#     data.frame(number=nrow(df), Fraction=fraction)
#   })
# }

## Overall sample sizes
## How to best generate tables? Directly to Latex?
library(stargazer)
fractions.all <- make.fractions(slr_structure_all, c("tm_full", "sec_simple"), thr=0.05)
fractions.all$Fraction <- NULL
stargazer(fractions.all, summary=FALSE)


## Mapped sequence length
abs_lengths <- ddply(slr_structure_raw, c("stable_id", "pdb_id", "pdb_chain"), nrow)

ggplot(subset(abs_lengths, V1 < 3000), aes(x=V1)) +
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank(), 
        axis.text.x=element_text(margin=margin(l=0.0, r=0.0, b=0.0, t=-0.3, unit="cm")), 
        axis.text.y=element_text(margin=margin(l=0.0, r=-0.3, b=0.0, t=0.0, unit="cm"))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Length of mapped region") +
  ylab("Count") +  
  geom_histogram(bins=60)
ggsave(file=paste(results, "1a_mapping_hist_abs.pdf", sep="/"), width=7, height=5.5, unit="cm")


## Remember the outlier!
max(abs_lengths$V1)

# For the "relative" mapping lengths, we need sequence lengths
seq_lengths <- dlply(slr_all, "stable_id", function(df) { max(df$human_idx) })

rel_lengths <- ddply(abs_lengths, c("stable_id", "pdb_id", "pdb_chain"), function(df) {
  data.frame(rel_length=(df$V[1] / seq_lengths[[df$stable_id[1]]]))
})

subset(rel_lengths, rel_length > 1)

ggplot(subset(rel_lengths, rel_length < 3000), aes(x=rel_length)) +
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank(), 
        axis.text.x=element_text(margin=margin(l=0.0, r=0.0, b=0.0, t=-0.3, unit="cm")), 
        axis.text.y=element_text(margin=margin(l=0.0, r=-0.3, b=0.0, t=0.0, unit="cm"))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Coverage of mapped region") +
  ylab("Count") +  
  geom_histogram(bins=60)
ggsave(file=paste(results, "1a_mapping_hist_rel.pdf", sep="/"), width=7, height=5.5, unit="cm")

## Counts of sites in bins of 0.1 OR deciles
## Need a table breakdown of results
# pdb_master_nona <- subset(pdb_master, !is.na(rsa))
## Deciles or equal intervals?
# pdb_master_nona$rsa_breaks <- cut(pdb_master_nona$rsa, breaks=seq(0.0, 1.0, 0.1), include.lowest=T)
slr_structure_all$rsa_breaks <- cut2(slr_structure_all$rsa, g=10)

bin_labels <- str_replace(as.character(levels(slr_structure_all$rsa_breaks)), ",", ",\n")
ggplot(subset(slr_structure_all, !is.na(rsa)), aes(y=omega, x=rsa_breaks)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete(labels=bin_labels) + 
  xlab("Relative solvent accessibility") +
  ylab(expression(omega)) +
  coord_cartesian(ylim=c(0, 5)) +
  geom_boxplot()# +
  # facet_grid(. ~ rsa_breaks)

ggsave(file=paste(results, "1_omega_RSA_bins.pdf", sep="/"), width=12, height=8, unit="cm")

# Ratio of means
ddply(slr_structure_all, "rsa_breaks", function(df) {mean(df$omega, na.rm=T)})

fractions <- make.fractions(subset(slr_structure_all, !is.na(rsa)), variables=c("rsa_breaks"), 0.05)
err_limits <- ddply(fractions, c("rsa_breaks"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=rsa_breaks, y=Fraction)) +
  theme_minimal(base_size=8) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Relative solvent accessibility") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("dodgerblue2", "dodgerblue4"), guide=guide_legend(title="Solvent\naccessibility")) +
  scale_x_discrete(labels=bin_labels) + 
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "1_omega_RSA_fractions.pdf", sep="/"), width=14, height=8, unit="cm")


## Just solvent acc.
rsa.fractions <- make.fractions(slr_structure_all, c("buried"), 0.05)
err_limits <- ddply(rsa.fractions, c("buried"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(rsa.fractions, aes(x=buried, y=Fraction, fill=buried)) +
  theme_minimal(base_size=8) +
  theme(legend.position="none",
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c("lightsteelblue4", "lightsteelblue3")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  xlab("RSA") +
  ylab("Fraction under positive selection") +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "2_buried_fractions.pdf", sep="/"), width=5, height=8, unit="cm")

## Overall fractions under purifying selection
fractions <- make.fractions.neg(slr_structure_globular, variables=c("buried"), 0.05)
err_limits <- ddply(fractions, c("buried"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=buried, y=Fraction, fill=buried)) +
  theme_minimal(base_size=8) +
  theme(legend.position="none",
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Solvent\naccessibility") +
  ylab("Fraction under purifying selection") +
  scale_y_continuous(limits=c(0, 1)) +
  scale_fill_manual(values=c("lightsteelblue4", "lightsteelblue3")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)

ggsave(file=paste(results, "2a_buried_fractions_neg.pdf", sep="/"), width=5, height=8, unit="cm")


fractions <- make.fractions.neg(slr_structure_globular, variables=c("buried", "disorder", "sec_simple"), 0.05)
err_limits <- ddply(fractions, c("buried", "disorder", "sec_simple"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=disorder)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Disorder") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3"), guide=guide_legend(title="Disorder")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25) +
  facet_grid(. ~ buried)
ggsave(file=paste(results, "4b_omega_fractions_all.pdf", sep="/"), width=14, height=10, unit="cm")



hist(tm_lengths)

ggplot(data.frame(tm_lengths=tm_lengths), aes(x=tm_lengths)) +
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  xlab("Length of TM region") +
  ylab("Count") +
  geom_histogram(binwidth=1)
ggsave(file=paste(results, "3_TM_length_hist.pdf", sep="/"), width=10, height=8, unit="cm")

summary(tm_lengths)

# RSA and \omega ==============================================================
fractions <- make.fractions(slr_structure_all, variables=c("buried", "tm_full"), 0.05)
err_limits <- ddply(fractions, c("buried", "tm_full"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=buried, y=Fraction, fill=tm_full)) +
  theme_minimal(base_size=8) +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(guide=guide_legend(title="Protein type"), values=c("indianred", "coral1", "firebrick", "darkred")) +
  xlab("Relative solvent accessibility") +
  ylab("Fraction under positive selection") +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
  
ggsave(paste(results, "3_tm_buried_fractions.pdf", sep="/"), width=10, height=8, unit="cm")

# Which are the 
sort(table(subset(slr_structure_all, tm_full == "Transmembrane (protein)" & omega > 1 & Adj.Pval < 0.05)$stable_id))

subset(slr_structure_all, tm_full == "Transmembrane (protein)" & omega > 1 & Adj.Pval < 0.05 & stable_id == "ENSP00000357036") 
## How many unique TM proteins are there
## The sample size is very small for the 'actually' transmembrane
length(unique(subset(slr_structure_all, tm == TRUE)$stable_id))
## Secondary structure
## Separate and/or joint with solvent accessibility
fractions <- make.fractions(slr_structure_all, c("buried", "sec_simple"), thr=0.05)
err_limits <- ddply(fractions, c("buried", "sec_simple"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=buried)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Secondary structure") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("dodgerblue2", "dodgerblue4"), guide=guide_legend(title="Solvent\naccessibility")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
  
ggsave(file=paste(results, "4_buried_ss_fractions.pdf", sep="/"), width=10, height=8, unit="cm")

## Domain stuff
## Hard to do a further division but we should be able to do this by facets
fractions <- make.fractions(slr_structure_globular, variables=c("buried", "sec_simple", "is_domain"), 0.05)
err_limits <- ddply(fractions, c("buried", "sec_simple", "is_domain"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=is_domain)) +
  theme_minimal(base_size=8) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_fill_manual(guide=guide_legend(title="Region"), values=c("dodgerblue2", "dodgerblue4")) +
  xlab("") +
  ylab("Fraction under positive selection") +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25) +
  facet_grid(. ~ buried)
ggsave(file=paste(results, "5_omega_domains.pdf", sep="/"), width=14, height=8, unit="cm")

## NEW (corrections): Only non-gapped region as a sanity check
fractions <- make.fractions(subset(slr_structure_globular, fr_aligned == 1), variables=c("buried", "sec_simple", "is_domain"), 0.05)
err_limits <- ddply(fractions, c("buried", "sec_simple", "is_domain"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=is_domain)) +
  theme_minimal(base_size=8) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_fill_manual(guide=guide_legend(title="Region"), values=c("dodgerblue2", "dodgerblue4")) +
  xlab("Solvent exposure") +
  ylab("Fraction under positive selection") +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25) +
  facet_grid(. ~ buried)
ggsave(file=paste(results, "5_omega_domains_nogaps.pdf", sep="/"), width=14, height=8, unit="cm")


## Histograms of the distribution of solvent accessibility
# 
# ggplot(slr_structure_globular, aes(x=rsa, fill=is_domain)) +
#   theme_minimal(base_size=8) +
#   theme_minimal() +
#   xlab("Relative solvent accessibility") +
#   theme(panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()) +
#   scale_fill_discrete(guide=guide_legend(title="Region")) +
# geom_density(alpha=0.4)
# ggsave(file=paste(results, "5a_RSA_domains.pdf", sep="/"), width=13, height=11, unit="cm")


ggplot(slr_structure_globular, aes(x=rsa, fill=is_domain)) +
  theme_minimal(base_size=8) +
  xlab("Relative solvent accessibility") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_discrete(guide=guide_legend(title="Region")) +
  geom_density(alpha=0.4) +
  facet_grid(. ~ buried, scales='free_x')
ggsave(file=paste(results, "5a_RSA_domains.pdf", sep="/"), width=14, height=7, unit="cm")


# 2D histogram with marginals
main_hist <- ggplot(slr_structure_all, aes(x=rsa, y=depth)) +
  theme_minimal(base_size=8) +
  theme(legend.position="none") +
  scale_alpha_continuous(limits=c(0,2.5),breaks=seq(0,2.5,by=0.00002)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  stat_density_2d(aes(alpha=..level..), geom = "polygon", h=c(0.03, 0.015))
  # geom_density_2d()

marg_rsa <- ggplot(slr_structure_all, aes(x=rsa, y=..density..)) +
  theme_minimal(base_size=8) +  
  theme(plot.margin=margin(t=0, r=0.5, l=0, b=0, unit="cm"),
        panel.margin=margin(0, 0, 0, 0),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  geom_histogram()

marg_rd <- ggplot(slr_structure_all, aes(x=depth, y=..density..)) +
  coord_flip() +
  theme_minimal(base_size=8) +  
  theme(plot.margin=margin(t=0, r=0, l=0, b=0, unit="cm"),
        panel.margin=margin(0, 0, 0, 0),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0.0, 0.5), breaks=c(0.0, 0.25, 0.5)) +
  geom_histogram()

empty <- ggplot()+geom_point(aes(1,1), colour="white") +
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())


plot_grid(main_hist, marg_rd, marg_rsa, empty, nrow=2, ncol=2, align="hv", rel_heights=c(1, 0.25), rel_widths=c(1, 0.25)) 
ggsave(file=paste(results, "6_rsa_rd_2dhist.pdf", sep="/"), width=16, height=16, unit="cm")


## Depth vs. omega
ggplot(subset(slr_structure_all, omega <= 5), aes(x=depth, y=omega)) +
  theme_minimal(base_size=8) +
  theme(legend.position="none") +
  # scale_alpha_continuous(limits=c(0,1.0),breaks=seq(0,1.0,by=0.002)) +
  # scale_x_continuous(expand=c(0, 0)) +
  # scale_y_continuous(expand=c(0, 0)) +
  stat_density_2d(aes(alpha=..level..), geom = "polygon")

slr_structure_buried <- subset(slr_structure_all, rsa < 0.25 & !is.na(depth))

##
## Residue depth
##
# max_depth <- max(slr_structure_buried$depth)
# slr_structure_buried$depth_breaks <- cut(slr_structure_buried$depth, breaks=seq(0.0, max_depth, max_depth/10), include.lowest=T)
slr_structure_all$depth_breaks <- cut2(slr_structure_all$depth, g=10)
bin_labels <- str_replace(as.character(levels(slr_structure_all$depth_breaks)), ",", ",\n")
ggplot(subset(slr_structure_all, !is.na(depth)), aes(y=omega, x=depth_breaks)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete(labels=bin_labels) + 
  xlab("Residue depth") +
  ylab(expression(omega)) +
  coord_cartesian(ylim=c(0, 5)) +
  geom_boxplot()
ggsave(file=paste(results, "1c_omega_depth_hist.pdf", sep="/"), width=14, height=10, unit="cm")


## NEW Sequence entropy as a function of depth
ggplot(subset(slr_structure_all, !is.na(depth)), aes(y=entropy, x=depth_breaks)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete(labels=bin_labels) + 
  xlab("Residue depth") +
  ylab("Sequence entropy") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()
ggsave(file=paste(results, "entropy_depth.pdf", sep="/"), width=14, height=10, unit="cm")


## NEW Fr. aligned as a function of depth
ggplot(subset(slr_structure_all, !is.na(depth)), aes(y=fr_aligned, x=depth_breaks)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete(labels=bin_labels) + 
  xlab("Residue depth") +
  ylab("Sequence entropy") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()
ggsave(file=paste(results, "aligned_depth.pdf", sep="/"), width=14, height=10, unit="cm")


fractions <- make.fractions(subset(slr_structure_all, !is.na(depth)), variables=c("depth_breaks"), 0.05)
err_limits <- ddply(fractions, c("depth_breaks"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=depth_breaks, y=Fraction)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    scale_x_discrete(labels=bin_labels) + 
    xlab("Residue depth") +
  ylab("Fraction under positive selection") +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "1d_omega_depth_fractions_all.pdf", sep="/"), width=14, height=10, unit="cm")

## Depth fractions, just buried?
slr_structure_buried <- subset(slr_structure_all, !is.na(depth) & rsa < 0.1)
slr_structure_buried$depth_breaks <- cut2(slr_structure_buried$depth, g=10)
bin_labels <- str_replace(as.character(levels(slr_structure_buried$depth_breaks)), ",", ",\n")
fractions <- make.fractions(slr_structure_buried, variables=c("depth_breaks"), 0.05)
err_limits <- ddply(fractions, c("depth_breaks"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=depth_breaks, y=Fraction)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_x_discrete(labels=bin_labels) + 
  xlab("Residue depth") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("dodgerblue2", "dodgerblue4"), guide=guide_legend(title="Solvent\naccessibility")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "1d_omega_depth_fractions_buried.pdf", sep="/"), width=14, height=10, unit="cm")


## Specific examples
sort(table(subset(slr_structure_all, omega > 1.0 & Adj.Pval < 0.05 & depth_breaks == "[8.41,19.10]")$stable_id))

## DEXD/H-box helicase 58
subset(slr_structure_all, omega > 1.0 & Adj.Pval < 0.05 & depth_breaks == "[8.41,19.10]" & stable_id == "ENSP00000369213")

## glyceraldehyde-3-phosphate dehydrogenase 
subset(slr_structure_all, omega > 1.0 & Adj.Pval < 0.05 & depth_breaks == "[8.41,19.10]" & stable_id == "ENSP00000229239")
subset(slr_structure_all, omega > 1.0 & Adj.Pval < 0.05 & stable_id == "ENSP00000229239")

## Aldo-keto reductase family 1 member C4
subset(slr_structure_all, omega > 1.0 & Adj.Pval < 0.05 & depth_breaks == "[8.41,19.10]" & stable_id == "ENSP00000369814")
subset(slr_structure_all, omega > 1.0 & Adj.Pval < 0.05 & stable_id == "ENSP00000369814")


## Disorder
## How many sites compared to total
slr_structure_disorder <- subset(slr_structure_globular, disorder == TRUE)
nrow(slr_structure_disorder)/nrow(slr_structure_globular)
table(slr_structure_disorder$secondary)

sum(slr_structure_disorder$omega > 1 & slr_structure_disorder$Adj.Pval < 0.05 & slr_structure_disorder$buried == "exposed") / sum(slr_structure_disorder$buried == "exposed")

## The overall distribution
ggplot(subset(slr_structure_globular, omega < 5), aes(x=omega, fill=disorder)) +
  theme_minimal(base_size=8) +
  xlab("Relative solvent accessibility") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(0, 1)) +
  geom_density(alpha=0.4)

## Need to cobble together
# - all disordered sites (regardless of structure and of solvent accessibility)
# - all not-disordered sites (with structure but not disordered), divided into
# -- buried
# -- exposed
disorder_combined <- subset(slr_structure_globular, disorder == FALSE)[, c("stable_id", "ens_pos", "omega", "buried", "disorder")]
disorder_tmp <- subset(slr_all, disorder==TRUE)[, c("stable_id", "ens_pos", "Omega", "disorder")]
disorder_tmp$omega <- disorder_tmp$Omega
disorder_tmp$Omega <- NULL
disorder_tmp$buried <- NA

disorder_combined <- rbind(disorder_combined, disorder_tmp)
ggplot(disorder_combined, aes(y=omega, x=buried)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  xlab("") +
  ylab("Fraction aligned") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()

fractions <- make.fractions(slr_structure_globular, variables=c("buried", "disorder"), 0.05)
err_limits <- ddply(fractions, c("buried", "disorder"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=buried, y=Fraction, fill=disorder)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Solvent accessibility") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3"), guide=guide_legend(title="Disorder")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "6a_omega_disorder_fractions.pdf", sep="/"), width=7, height=8, unit="cm")


fractions <- make.fractions(slr_structure_globular, variables=c("buried", "disorder", "sec_simple"), 0.05)
err_limits <- ddply(fractions, c("buried", "disorder", "sec_simple"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=disorder)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Secondary structure") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3"), guide=guide_legend(title="Disorder")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25) +
  facet_grid(. ~ buried)
ggsave(file=paste(results, "6b_omega_disorder_fractions_all.pdf", sep="/"), width=14, height=8, unit="cm")


## Disorder -- fractions under purifying selection
fractions <- make.fractions.neg(slr_structure_globular, variables=c("buried", "disorder"), 0.05)
err_limits <- ddply(fractions, c("buried", "disorder"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=buried, y=Fraction, fill=disorder)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Solvent accessibility") +
  ylab("Fraction under purifying selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3"), guide=guide_legend(title="Disorder")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "6a_omega_disorder_fractions_neg.pdf", sep="/"), width=7, height=8, unit="cm")

##
## Plots for corrections
##

## Arguably a simpler thing to look at is the difference between buried and exposed regions
ggplot(slr_structure_all, aes(y=fr_aligned, x=buried)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  xlab("") +
  ylab("Fraction aligned") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()
ggsave(file=paste(results, "aligned_rsa.pdf", sep="/"), width=10, height=10, unit="cm")

ggplot(slr_structure_all, aes(y=entropy, x=buried)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  xlab("") +
  ylab("Sequence entropy") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()
ggsave(file=paste(results, "entropy_rsa.pdf", sep="/"), width=10, height=10, unit="cm")

# NEW -- Fraction aligned / domain
ggplot(slr_structure_all, aes(y=fr_aligned, x=is_domain, fill=is_domain)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        legend.position="none") +
  scale_x_discrete() + 
  xlab("") +
  ylab("Fraction aligned") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_violin(bw=0.05, kernel="biweight") +
  geom_boxplot(width=0.1, fill="white", outlier.shape=NA, size=0.25, position=position_dodge(0.9)) +
  facet_grid(. ~ buried)
ggsave(file=paste(results, "domain_fr_aligned.pdf", sep="/"), width=8, height=10, unit="cm")


## make.fractions(subset(slr_structure_all, fr_aligned==1), c("buried", "is_domain"), 0.05)


# NEW -- Entropy / domain
ggplot(slr_structure_all, aes(y=entropy, x=buried, fill=is_domain)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  xlab("") +
  ylab("Sequence entropy") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()
ggsave(file=paste(results, "entropy_domain.pdf", sep="/"), width=10, height=10, unit="cm")

## First calculate the fraction disordered
nrow(subset(slr_all, disorder == TRUE & stable_id %in% slr_structure_globular$stable_id))/nrow(subset(slr_all, stable_id %in% slr_structure_globular$stable_id))

disorder <- subset(slr_all, disorder == TRUE & stable_id %in% slr_structure_globular$stable_id)[, c("stable_id", "ens_pos", "Omega", "Pval", "fr_aligned", "entropy")]
colnames(disorder)[3] <- "omega"
disorder$status <- "disordered"

structured <- subset(slr_structure_all, disorder == FALSE & buried == "exposed")[, c("stable_id", "ens_pos", "omega", "Pval", "fr_aligned", "entropy")]
structured$status <- "exposed"
disorder <- rbind(disorder, structured)

structured <- subset(slr_structure_all, disorder == FALSE & buried == "buried")[, c("stable_id", "ens_pos", "omega", "Pval", "fr_aligned", "entropy")]
structured$status <- "buried"
disorder <- rbind(disorder, structured)

disorder$Adj.Pval <- p.adjust(disorder$Pval, "BH")

ggplot(disorder, aes(y=omega, x=status, fill=status)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  scale_fill_manual(values=c("darkorange4", "darkorange3", "orange"), guide=guide_legend(title="")) +
  theme(legend.position="none") +
  xlab("Type of structure") +
  ylab("omega") +
  coord_cartesian(ylim=c(0, 5)) +
  geom_violin(bw=0.01, kernel="biweight") +
  geom_boxplot(width=0.1, outlier.shape=NA, fill="white", size=0.25)
ggsave(file=paste(results, "omega_dist_disorder.pdf", sep="/"), width=7, height=9, unit="cm")

ddply(disorder, "status", function(df) mean(df$omega))

##
## Fraction
##
fractions <- make.fractions(disorder, c("status"), 0.05)
err_limits <- ddply(fractions, c("status"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(err_limits, aes(x=status, y=Fraction, fill=status)) +
  theme_minimal(base_size=10) +
  theme(legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3", "orange"), guide=guide_legend(title="Type of structure")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "6b_disorder_fractions_fixed.pdf", sep="/"), width=8, height=10, unit="cm")

## Only fully aligned positions
fractions <- make.fractions(subset(disorder, fr_aligned==1), c("status"), 0.05)
err_limits <- ddply(fractions, c("status"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(err_limits, aes(x=status, y=Fraction, fill=status)) +
  theme_minimal(base_size=10) +
  theme(legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Type of structure") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3", "orange"), guide=guide_legend(title="Type of structure")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "6b_disorder_fractions_fixed_nogaps.pdf", sep="/"), width=8, height=10, unit="cm")

##
## Now the fraction under purifying selection
##
fractions <- make.fractions.neg(disorder, variables=c("status"), 0.05)
err_limits <- ddply(fractions, c("status"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=status, y=Fraction, fill=status)) +
  theme_minimal(base_size=10) +
  theme(legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Type of structure") +
  ylab("Fraction under purifying selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3", "orange")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "6b_disorder_fractions_neg_fixed.pdf", sep="/"), width=8, height=10, unit="cm")


## Fraction purifying selection disorder, fully aligned
fractions <- make.fractions.neg(subset(disorder, fr_aligned==1), variables=c("status"), 0.05)
err_limits <- ddply(fractions, c("status"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=status, y=Fraction, fill=status)) +
  theme_minimal(base_size=10) +
  theme(legend.position="none",
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Type of structure") +
  ylab("Fraction under purifying selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3", "orange")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
ggsave(file=paste(results, "6b_disorder_fractions_neg_fixed_nogaps.pdf", sep="/"), width=8, height=10, unit="cm")

##
##
##
ggplot(disorder, aes(y=fr_aligned, x=status, fill=status)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  scale_fill_manual(values=c("darkorange4", "darkorange3", "orange")) +
  theme(legend.position="none") +
  xlab("Type of structure") +
  ylab("Fraction aligned") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_violin(bw=0.05, kernel="biweight") +
  geom_boxplot(width=0.1, outlier.shape=NA, fill="white", size=0.25)
ggsave(file=paste(results, "disorder_fr_aligned.pdf", sep="/"), width=7, height=9, unit="cm")


ggplot(disorder, aes(y=entropy, x=status)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  scale_x_discrete() + 
  xlab("Sequence entropy") +
  ylab("omega") +
  coord_cartesian(ylim=c(0, 1)) +
  geom_boxplot()

## Ugo's stuff
# omega vs. sequence entropy
with(slr_structure_all, plot(omega, entropy, xlim=c(0, 10)))

ggplot(subset(slr_structure_all, omega < 5 & omega != 0), aes(x=omega, y=entropy)) +
  theme_minimal(base_size=8) +
  theme(legend.position="none") +
  scale_x_continuous(limits=c(0, 5), expand=c(0, 0)) +
  scale_y_continuous(limit=c(0, 3), expand=c(0, 0)) +
  xlab(expression(omega)) +
  ylab("Sequence entropy") +
  geom_point(alpha=0.02) +
  geom_smooth(method = "gam", size = 1.5, method.args=list(bw=0.01))

ggsave(file=paste(results, "entropy_omega.pdf", sep="/"), width=10, height=10, unit="cm")

##
## Ugo's specific comments about entropy UB-12 and UB-13
##
ggplot(slr_structure_globular, aes(x=sec_simple, y=entropy, fill=buried)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Secondary structure") +
  ylab("Sequency entropy") +
  scale_fill_manual(values=c("dodgerblue2", "dodgerblue4"), guide=guide_legend(title="Solvent\naccessibility")) +
  geom_boxplot(position=position_dodge(0.9)) # +
  # geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)

ggplot(slr_structure_all, aes(x=buried, y=entropy, fill=tm_full)) +
  theme_minimal(base_size=8) +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(guide=guide_legend(title="Protein type"), values=c("indianred", "coral1", "firebrick", "darkred")) +
  xlab("Relative solvent accessibility") +
  ylab("Sequence entropy") +
  geom_boxplot(position=position_dodge())

##
## Do a gam fit
##
ggplot(subset(slr_structure_all, !is.na(rsa)), aes(y=omega, x=rsa)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  # scale_x_discrete(labels=bin_labels) + 
  xlab("Relative solvent accessibility") +
  ylab(expression(omega)) +
  coord_cartesian(ylim=c(0, 5)) +
  geom_point(size=0.2, alpha=0.05, colour="black") +
  geom_smooth(method = "gam", size = 1.5, method.args=list(bw=0.01))
ggsave(file=paste(results, "omega_rsa_new.pdf", sep="/"), width=12, height=8, unit="cm")

ggplot(subset(slr_structure_all, !is.na(depth)), aes(y=omega, x=depth)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  # scale_x_discrete(labels=bin_labels) + 
  xlab("Residue depth") +
  ylab(expression(omega)) +
  coord_cartesian(ylim=c(0, 5)) +
  geom_point(size=0.2, alpha=0.05, colour="black") +
  geom_smooth(size = 1.5, method.args=list(bw=0.01))
ggsave(file=paste(results, "omega_rd_new.pdf", sep="/"), width=12, height=8, unit="cm")

nrow(subset(slr_structure_all, buried == "buried" & omega > 0.9 & Adj.Pval < 0.05) )/nrow(subset(slr_structure_all, buried == "buried" & omega > 0.9))
nrow(subset(slr_structure_all, buried == "exposed" & omega > 0.9) )/nrow(subset(slr_structure_all, buried == "exposed" & omega > 0.9))


fractions <- make.fractions.neg(slr_structure_globular, variables=c("buried", "disorder", "sec_simple"), 0.05)
err_limits <- ddply(fractions, c("buried", "disorder", "sec_simple"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=disorder)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Disorder") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("darkorange4", "darkorange3"), guide=guide_legend(title="Disorder")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25) +
  facet_grid(. ~ buried)
ggsave(file=paste(results, "6b_omega_disorder_fractions_neg_all.pdf", sep="/"), width=14, height=8, unit="cm")

## What is up with omega when entropy != 0
## There are multiple substitutions which means that the codon model would infer 
## intermediate substitutions through codons coding for other amino acids.
dodgy_omega <- subset(slr_structure_all, omega != 0 & entropy == 0)
dodgy_omega[, c("stable_id", "Site", "dataset", "omega", "entropy")]

## 
## Simon's comment about unaligned regions
##
ggplot(slr_structure_all, aes(x=fr_aligned)) +
  theme_minimal(base_size=10, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Fraction of aligned residues") +
  ylab("Count") +  
  geom_histogram(bins=60)
ggsave(file=paste(results, "fr_aligned_hist.pdf", sep="/"), width=12, height=8, unit="cm")

sum(slr_structure_all$fr_aligned < 0.5)/nrow(slr_structure_all)
