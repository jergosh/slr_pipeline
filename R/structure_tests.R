# To do the tests exhaustively (within) 

ks.test(subset(pdb_master_globular, sec_simple == "Helix" & buried == "exposed" & is_domain == FALSE)$omega,
        subset(pdb_master_globular, sec_simple == "Beta sheet" & buried == "exposed" & is_domain == FALSE)$omega)

hist(subset(pdb_master_globular, sec_simple == "Helix" & buried == "exposed" & is_domain == TRUE)$omega)
hist(subset(pdb_master_globular, sec_simple == "Helix" & buried == "exposed" & is_domain == FALSE)$omega)

# Plots -----------------------------------------------------------------------

# Histograms of the overall distribution --------------------------------------
ggplot(pdb_master_globular, aes(x=omega, y=..density.., fill=sec_simple)) +
  coord_cartesian(xlim=c(0, 5)) + 
  geom_histogram(binwidth=0.1, position="dodge")

ggplot(subset(pdb_master, omega > 0), aes(x=omega)) +
  stat_bin(binwidth=0.1, geom="line") +
  theme_bw()
# geom_histogram(binwidth=0.1)

ggplot(subset(pdb_master, omega > 0), aes(x=omega)) +
  stat_bin(binwidth=0.1, geom="line") +
  theme_bw()

# Is the distribution on sites with structures the same as the overall distribution?

# ! Contemplative interlude ===================================================
## Not sure how to make stuff visible on the plots despite the unusual distribution of omega
## One idea is to remove the omega == 0 from the distribution and plot them separately on the negative 
## axis.
## Another way would be to pre-filter the alignments in a way that decreases the fraction of sites with omega == 0
##

pdb_master_nonzero <- subset(pdb_master, omega != 0)
# Both secondary structure and RSA
ggplot(pdb_master_nonzero, aes(sec_simple, omega, colour=buried)) +
  scale_y_continuous(lim=c(0, 2)) +
  geom_boxplot()

ggplot(pdb_master_nonzero, aes(buried, omega, colour=sec_simple)) +
  scale_y_continuous(lim=c(0, 2)) +
  geom_boxplot()

# Density plots of omega ======================================================
# Uncropped
ggplot(pdb_master, aes(x=omega, linetype=buried, colour=sec_simple)) +
  geom_density()

# Cropped
ggplot(pdb_master, aes(x=omega, linetype=buried, colour=sec_simple)) +
  coord_cartesian(xlim=c(0, 5), ylim=c(0, 10)) +
  geom_density()

# Overall RSA distribution ====================================================
## Distributions of RSA independent of omega
## (presumably for supplementary figures)
ggplot(pdb_master, aes(x=rsa, colour=sec_simple)) +
  theme_bw() +
  geom_density()

plot.rsa.dist <- function(pdb_master) {
  ggplot(pdb_master, aes(x=rsa, colour=sec_simple)) +
    theme_bw() +
    geom_density()
}

plot.rsa.dist(subset(pdb_master, tm_full == "Transmembrane (sites)"))
plot.rsa.dist(subset(pdb_master, tm_full == "Transmembrane (protein)"))
plot.rsa.dist(subset(pdb_master, tm_full == "Globular"))

## Master plot
pdf(paste(results,"RSA_master.pdf", sep="/"), height=10, width=14)
ggplot(pdb_master, aes(x=rsa, colour=sec_simple, linetype=tm_full)) +
  theme_bw() +
  geom_density()
dev.off()

# Focus on the globular subset and break down on domains
pdf(paste(results,"RSA_globular.pdf", sep="/"), height=10, width=14)
ggplot(pdb_master_globular, aes(x=rsa, colour=sec_simple, linetype=is_domain)) +
  theme_bw() +
  geom_density()
dev.off()




# Main plot \omega breakdown ==================================================

plot.fractions.everything <- function(fractions) {
  ggplot(fractions, aes(x=buried, y=Fraction, fill=sec_simple)) +
    #     scale_fill_manual(values=c("lightblue", "blue", "darkblue",
    #                                "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
    #                                "brown1", "brown3", "brown4")) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle("Positive selection and structure") +
    theme_bw(base_size=18)
}

plot.fractions.everything.se <- function(fractions) {
  fractions$i <- 1:nrow(fractions)
  ggplot(fractions, aes(x=buried, y=Fraction, fill=structure)) +
    #     scale_fill_manual(values=c("lightblue", "blue", "darkblue",
    #                                "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
    #                                "brown1", "brown3", "brown4")) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_segment(x=i) +
    ggtitle("Positive selection and structure") +
    theme_bw(base_size=18)
}

sebinom <- function(p, n) { sqrt(p*(1-p)/n) }
fractions.everything.globular$se <- sebinom(fractions.everything.globular$Fraction, fractions.everything.globular$number)

slr_structure_all_globular <- subset(slr_structure_all_strict, tm_full == "Globular")
fractions.everything <- make.fractions(slr_structure_all_globular, c("buried", "sec_simple"))
plot.fractions.everything(fractions.everything)
ggsave(paste(results, "omega_main_everything.pdf", sep="/"), width=11, height=8)

fractions.everything.globular <- make.fractions(subset(pdb_master, tm_full == "Globular"), c("buried", "structure"))
plot.fractions.everything(fractions.everything.globular)
ggsave(paste(results, "omega_main_globular.pdf", sep="/"), width=11, height=8)
fractions.everything.tm <- make.fractions(subset(pdb_master, tm_full != "Globular"), c("buried", "structure"))
plot.fractions.everything(fractions.everything.tm)
ggsave(paste(results, "omega_main_tm.pdf", sep="/"), width=11, height=8)

fractions.everything.globular <- make.fractions.lower(subset(slr_structure, tm_full == "Globular"), c("buried", "structure"))
plot.fractions.everything(fractions.everything.globular)
ggsave(paste(results, "omega_lower_main_globular.pdf", sep="/"), width=11, height=8)


# TODO: Have a look at the coefficient of variation

# Look at whether enzyme classification makes a difference --------------------
plot.is.enzyme <- function(fractions, main) {
  ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=is_enzyme)) +
    scale_fill_manual(values=c("lightblue", "blue", "darkblue",
                               "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
                               "brown1", "brown3", "brown4")) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(main) +
    theme_bw(base_size=18)
}

fractions.enzyme <- make.fractions(pdb_master_globular, c("sec_simple", "is_enzyme"))
plot.is.enzyme(fractions.enzyme, "Positive selection and enzyme classification (all)")
ggsave(paste(results, "omega_enzyme_all.pdf", sep="/"), width=11, height=8)
fractions.enzyme.exposed <- make.fractions(subset(pdb_master_globular, buried=="exposed"), c("sec_simple", "is_enzyme"))
plot.is.enzyme(fractions.enzyme.exposed, "Positive selection and enzyme classification (exposed)")
ggsave(paste(results, "omega_enzyme_exposed.pdf", sep="/"), width=11, height=8)
fractions.enzyme.buried <- make.fractions(subset(pdb_master_globular, buried=="buried"), c("sec_simple", "is_enzyme"))
plot.is.enzyme(fractions.enzyme.buried, "Positive selection and enzyme classification (buried)")
ggsave(paste(results, "omega_enzyme_buried.pdf", sep="/"), width=11, height=8)


fractions.everything <- make.fractions(pdb_master, c("buried", "is_enzyme"))
plot.fractions.everything(fractions.everything)

secondary.rsa.fractions <- make.fractions(pdb_master, variables=c("secondary", "buried"))
secondary.rsa.fractions <- subset(secondary.rsa.fractions, secondary != "I")

# omega_sec_buried <- function()
# Bars grouped by buriedness
ggplot(secondary.rsa.fractions, aes(x=buried, y=Fraction, fill=secondary)) +
  # scale_fill_manual(values=c("lightblue", "blue", "darkblue", "grey")) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Positive selection and secondary structure") +
  theme_bw(base_size=18)


# Now just the globular proteins ==============================================
pdb_master_globular <- subset(pdb_master, tm_full == "Globular")
head(pdb_master_globular)

fractions.globular <- make.fractions(pdb_master_globular, c("sec_simple", "buried"))

ggplot(fractions.globular, aes(x=buried, y=Fraction, fill=sec_simple)) +
  # scale_fill_manual(values=c("lightblue", "blue", "darkblue")) +
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Positive selection and secondary structure") +
  theme_bw(base_size=18)

# Density
ggplot(pdb_master, aes(x=omega, linetype=buried, colour=sec_simple)) +
  geom_density()
ggsave(paste(results, "omega_density_structure.pdf", sep="/"), width=11, height=8)

ggplot(pdb_master, aes(x=omega, linetype=buried, colour=sec_simple)) +
  coord_cartesian(xlim=c(0, 5), ylim=c(0, 10)) +
  geom_density()
ggsave(paste(results, "omega_density_structure_cropped.pdf", sep="/"), width=11, height=8)


# Or, better, histograms
ggplot(pdb_master_globular, aes(x=omega, linetype=buried, colour=sec_simple)) +
  stat_bin(aes(y=..density..), binwidth=0.1, geom="line") +
  theme_bw()

ggplot(pdb_master_globular, aes(x=omega, linetype=buried, colour=sec_simple)) +
  coord_cartesian(xlim=c(0, 5), ylim=c(0, 2.5)) +
  stat_bin(aes(y=..density..), binwidth=0.1, geom="line") +
  theme_bw()

# TODO The area under omega > 1 looks suspicious, let's look at the sample sizes

# Scatterplots of omega vs. RSA ----
ggplot(subset(pdb_master, omega > 1 & omega < 5), aes(x=rsa, y=omega, colour=sec_simple)) +
  geom_point(alpha=0.1)
ggsave(paste(results, "omega_scatterplot_cropped.pdf", sep="/"))

ggplot(pdb_master, aes(x=rsa, y=omega)) +
  geom_point() + 
  geom_density2d()

ggplot(subset(pdb_master, omega > 0), aes(x=rsa, y=omega)) +
  stat_density2d(aes(fill = ..level..), geom="polygon")

# 
##
## Disorder predictions?
##
## ...

## Need to make sure these aren't just proportional to sample sizes 

