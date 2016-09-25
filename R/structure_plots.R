sd_binom <- function(n, p) { sqrt(p * (1 - p)/(n - 1)) }

make.fractions <- function(master.df, variables, thr) {
  ddply(master.df, .variables=variables, function(df) {
    sum_pos <- sum(df$omega > 1.0 & df$Adj.Pval < thr)
    fraction <- sum_pos / nrow(df)
    data.frame(number=nrow(df), pos=sum_pos, Fraction=fraction)
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

## Counts of sites in bins of 0.1
## Need a table breakdown of results
pdb_master_nona <- subset(pdb_master, !is.na(rsa))
pdb_master_nona$rsa_breaks <- cut(pdb_master_nona$rsa, breaks=seq(0.0, 1.0, 0.1), include.lowest=T)
ggplot(pdb_master_nona, aes(x=rsa_breaks, y=omega)) + 
  theme_minimal(base_size=10, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  xlab("RSA range") +
  ylab(expression(omega)) +
  coord_cartesian(ylim=c(0, 5)) +
  geom_boxplot(outlier.size=1, size=0.25)

ggsave(file=paste(results, "omega_RSA_bins.pdf", sep="/"), width=14, height=12, unit="cm")


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
ggsave(file=paste(results, "TM_length_hist.pdf", sep="/"), width=10, height=8, unit="cm")

## Just solvent acc.
summary(tm_lengths)

rsa.fractions <- make.fractions(slr_structure_all, c("buried"), 0.05)
ggplot(rsa.fractions, aes(x=buried, y=Fraction, fill=buried)) +
  theme_minimal(base_size=10) +
  theme(legend.position="none",
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values=c("grey41", "grey25")) +
  xlab("Solvent exposure") +
  ylab("Fraction under positive selection") +
  geom_bar(stat="identity", position=position_dodge())
ggsave(file=paste(results, "buried_fractions.pdf", sep="/"), width=7, height=8, unit="cm")

# RSA and \omega ==============================================================
fractions <- make.fractions(slr_structure_all, variables=c("buried", "tm_full"), 0.05)
err_limits <- ddply(fractions, c("buried", "tm_full"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=buried, y=Fraction, fill=tm_full)) +
  theme_minimal(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(guide=guide_legend(title="Protein type"), values=c("indianred", "firebrick", "darkred")) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)

  
ggsave(paste(results, "tm_buried_fractions.pdf", sep="/"))

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
        axis.line.y = element_line(color = "black")) +
  xlab("Solvent exposure") +
  ylab("Fraction under positive selection") +
  scale_fill_manual(values=c("dodgerblue2", "dodgerblue4"), guide=guide_legend(title="Solvent\naccessibility")) +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25)
  
ggsave(file=paste(results, "buried_ss_fractions.pdf", sep="/"), width=16, height=12, unit="cm")

## Domain stuff
## Hard to do a further division but we should be able to do this by facets
fractions <- make.fractions(slr_structure_globular, variables=c("buried", "sec_simple", "is_domain"), 0.5)
err_limits <- ddply(fractions, c("buried", "sec_simple", "is_domain"), function(df) {
  sd_sample <- sd_binom(df$number[1], df$Fraction[1])
  data.frame(Fraction=df$Fraction[1], sd=sd_sample, ymin=(df$Fraction[1]-sd_sample), ymax=(df$Fraction[1]+sd_sample))
})

ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=buried)) +
  theme_minimal(base_size=10) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")) +
  scale_fill_manual(values=c("dodgerblue2", "dodgerblue4")) +
  xlab("Solvent exposure") +
  ylab("Fraction under positive selection") +
  geom_bar(stat="identity", position=position_dodge(0.9)) +
  geom_errorbar(data=err_limits, aes(ymin=ymin, ymax=ymax), position = position_dodge(0.9), width = 0.25) +
  facet_grid(. ~ is_domain)



## Histograms of the distribution of solvent accessibility
# 
ggplot(subset(slr_structure_globular, rsa < 0.25), aes(x=rsa, fill=is_domain)) +
  theme_minimal() +
  geom_density(alpha=0.4)

ggplot(subset(slr_structure_globular, rsa >= 0.25), aes(x=rsa, fill=is_domain)) +
  theme_minimal() +
  geom_density(alpha=0.4)

