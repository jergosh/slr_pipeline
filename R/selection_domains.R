library(plyr)
library(ggplot2)

results_dir <- "results/2015-05-30"


head(feature_data)
domain_data <- subset(feature_data, type == "pfam" & domain_id != "")
nrow(domain_data)

domain_data_structure <- subset(domain_data, stable_id %in% slr_structure_all$stable_id)
nrow(domain_data_structure)

domain_type <- rep(NA, nrow(slr_structure_all))
ddply(domain_data_structure, c("stable_id", "start", "end", "domain_id"), function(df) {
  domain_type[slr_structure_all$stable_id == df$stable_id[1] &
                   slr_structure_all$ens_pos >= (df$start[1]) &
                   slr_structure_all$ens_pos <= (df$end[1])] <<- df$domain_id[1]
})
sum(!is.na(domain_type))

slr_structure_all$domain_type <- as.factor(domain_type)
slr_structure_all$is_domain <- factor("domain", levels=c("domain", "non-domain"))
slr_structure_all$is_domain[is.na(slr_structure_all$domain_type)] <- "non-domain"

# Overview plot of the fraction under positive sel. ---------------------------

plot.domain.fractions <- function(fractions, main) {
  p <- ggplot(fractions, aes(x=sec_simple, y=Fraction, fill=is_domain)) +
    #   scale_fill_manual(values=c("lightblue", "blue", "darkblue",
    #                              "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
    #                              "brown1", "brown3", "brown4")) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(main) +
    theme_bw(base_size=18)
  
  
  print(p)
}

fractions.domains <- make.fractions(slr_structure_all, c("sec_simple", "is_domain"), thr=0.05)
plot.domain.fractions(fractions.domains, "Positive selection and domains")
ggsave(paste(results, "Domains (all).pdf", sep="/"), height=8, width=11)
fractions.domains.exposed <- make.fractions(subset(slr_structure_all, buried=="exposed"), c("sec_simple", "is_domain"))
plot.domain.fractions(fractions.domains.exposed, "Positive selection and domains (exposed)")
ggsave(paste(results, "Domains (exposed).pdf", sep="/"), height=8, width=11)
fractions.domains.buried <- make.fractions(subset(slr_structure_all, buried=="buried"), c("sec_simple", "is_domain"))
plot.domain.fractions(fractions.domains.buried, "Positive selection and domains (buried)")
ggsave(paste(results, "Domains (buried).pdf", sep="/"), height=8, width=11)

# Now with lower bounds
fractions.domains <- make.fractions(slr_structure_all, c("sec_simple", "is_domain"))
plot.domain.fractions(fractions.domains, "Positive selection and domains")
ggsave(paste(results, "Domains lower (all).pdf", sep="/"), height=8, width=11)
fractions.domains.exposed <- make.fractions(subset(slr_structure_all, buried=="exposed"), c("sec_simple", "is_domain"))
plot.domain.fractions(fractions.domains.exposed, "Positive selection and domains")
ggsave(paste(results, "Domains lower (exposed).pdf", sep="/"), height=8, width=11)
fractions.domains.buried <- make.fractions(subset(slr_structure_all, buried=="buried"), c("sec_simple", "is_domain"))
plot.domain.fractions(fractions.domains.buried, "Positive selection and domains")
ggsave(paste(results, "Domains lower (buried).pdf", sep="/"), height=8, width=11)



plot.domain.fractions <- function(fractions, main) {
  ggplot(fractions, aes(x=secondary, y=Fraction, fill=is_domain)) +
    #   scale_fill_manual(values=c("lightblue", "blue", "darkblue",
    #                              "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
    #                              "brown1", "brown3", "brown4")) +
    geom_bar(stat="identity", position=position_dodge()) +
    ggtitle(main) +
    theme_bw(base_size=18)
}

fractions.domains <- make.fractions(slr_structure_all, c("secondary", "is_domain"))
plot.domain.fractions(fractions.domains, "Positive selection and domains")
ggsave(paste(results, "Domains (all).pdf", sep="/"), height=8, width=11)
fractions.domains.exposed <- make.fractions(subset(slr_structure_all, buried=="exposed"), c("secondary", "is_domain"))
plot.domain.fractions(fractions.domains.exposed, "Positive selection and domains")
ggsave(paste(results, "Domains (exposed).pdf", sep="/"), height=8, width=11)
fractions.domains.buried <- make.fractions(subset(slr_structure_all, buried=="buried"), c("secondary", "is_domain"))
plot.domain.fractions(fractions.domains.buried, "Positive selection and domains")
ggsave(paste(results, "Domains (buried).pdf", sep="/"), height=8, width=11)

# Density of \omega depending on being in a domain



# Find structures with domain and non-domain regions
domain.fractions <- ddply(slr_structure_all, c("pdb_id", "pdb_chain"), function(df) {
  data.frame(domain_frac=sum(df$is_domain)/nrow(df), rest_frac=sum(!df$is_domain)/nrow(df))
})

fraction.structured <- rep(NA, nrow(slr_structure_all))
# Break down domains 
domain.fraction.structured <- ddply(slr_structure_all, c("stable_id", "domain_type"), function(df) {
  fraction.structured[slr_structure_all$stable_id == df$stable_id[1] &
                        slr_structure_all$domain_type == df$domain_type[1]] <<- sum(df$sec_simple %in% c("Helix", "Beta sheet"))/nrow(df)
})

ggplot(domain.fraction.structured, aes(x=V1)) +
  geom_density()

slr_structure_all$fraction.str <- fraction.structured
slr_structure_all$fraction.str.bin <- as.factor(as.integer(slr_structure_all$fraction.str*10))


ggplot(slr_structure_all, aes(x=omega, linetype=sec_simple, colour=fraction.str.bin)) +
  coord_cartesian(xlim=c(0, 5), ylim=c(0, 10)) +
  geom_density()

fractions.structured <- make.fractions(subset(slr_structure_all, buried == "exposed"), c("fraction.str.bin", "sec_simple"))
ggplot(fractions.structured, aes(x=fraction.str.bin, y=Fraction, fill=sec_simple)) +
  geom_bar(stat="identity", position=position_dodge())
  

for(i in 0:10) {
  print(i)
  fractions.structured.subset <- make.fractions(subset(slr_structure_all, fraction.str.bin == 7), c("buried", "sec_simple"))
  ggplot(fractions.structured.subset, aes(x=buried, y=Fraction, fill=sec_simple)) +
    geom_bar(stat="identity", position=position_dodge())
}


mean_omegas <- ddply(slr_structure_all, c("buried", "sec_simple"), function(df) { data.frame(omega=mean(df$omega, na.rm=T)) })
omega_norm <- ddply(slr_structure_all, c("buried", "sec_simple"), function(df) {
  data.frame(omega_norm=(df$omega / subset(mean_omegas, buried == df$buried[1] & sec_simple == df$sec_simple[1])$omega))
})
slr_structure_all$omega_norm <- omega_norm$omega_norm
pdb_master_domains <- subset(slr_structure_all, is_domain == TRUE)
domain_breakdown <- ddply(pdb_master_domains, c("domain_type"), function(df) {
  data.frame(domain_omega=mean(df$omega_norm), n=length(unique(df$stable_id)), frac=(sum(df$omega > 1)/length(df$omega)))
})
head(domain_breakdown)
nrow(domain_breakdown)
domain_breakdown_5 <- (subset(domain_breakdown, n >= 5))
domain_breakdown_10 <- (subset(domain_breakdown, n >= 10))

pdf(paste(results_dir, "pfamilies_5_scatter.pdf", sep="/"))
plot(domain_breakdown_5$domain_omega, domain_breakdown_5$frac,
     main="Protein families (5 and more members)", xlab=expression(paste("mean ", omega)),
     ylab=expression(paste("fraction ", omega, " > 1")))
dev.off()

pdf(paste(results_dir, "pfamilies_10_scatter.pdf", sep="/"))
plot(domain_breakdown_10$domain_omega, domain_breakdown_10$frac,
     main="Protein families (10 and more members)", xlab=expression(paste("mean ", omega)),
     ylab=expression(paste("fraction ", omega, " > 1")))
dev.off()

gene_breakdown <- ddply(pdb_master_domains, c("stable_id"), function(df) {
  data.frame(domain_type=df$domain_type[1], mean_omega=mean(df$omega_norm), n=length(unique(df$stable_id)), frac=(sum(df$omega > 1)/length(df$omega)))
})

domain_breakdown <- ddply(gene_breakdown, c("domain_type"), function(df) {
  # FIXME Could also take into account the site counts when calculating mean(s)
  data.frame(mean_omega=mean(df$mean_omega), sd_omega=sd(df$mean_omega), mean_frac=mean(df$frac), sd_frac=sd(df$frac),
             size=nrow(df))
})

ggplot(subset(domain_breakdown, size >= 5), aes(x=mean_omega, y=mean_frac, size=size)) +
  geom_point(fill="grey", alpha=0.5) +
  geom_segment(mapping=aes(x=mean_omega-sd_omega/10, xend=mean_omega+sd_omega/10, y=mean_frac, yend=mean_frac), size=0.1, col="red") +
  geom_segment(mapping=aes(x=mean_omega, xend=mean_omega, y=mean_frac-sd_frac/10, yend=mean_frac+sd_frac/10), size=0.1, col="red") +  
  theme_bw()

pdf(paste(results_dir, "pfamilies_genewise_scatter.pdf", sep="/"))
plot(gene_breakdown$mean_omega, gene_breakdown$frac,
     main="Protein families (single genes)", xlab=expression(paste("mean ", omega)),
     ylab=expression(paste("fraction ", omega, " > 1")))
dev.off()

# TODO How to include the pairwise correlations in the plot?

# TODO Next step could be to plot the fractions of each sec. structure in a piechart

# identify(domain_breakdown_5$domain_omega, domain_breakdown_5$frac, labels=domain_breakdown_5$domain_type)
# identify(domain_breakdown_10$domain_omega, domain_breakdown_10$frac, labels=domain_breakdown_10$domain_type)

domain_breakdown_5[order(domain_breakdown_5$domain_omega, decreasing=T)[1:10], ]
# The top domain is ubiquitin
domain_ubiquitin <- subset(slr_structure_all, domain_type == "PF00240")
sum(domain_ubiquitin$omega > 1)
subset(slr_all, stable_id == "ENSP00000405965")

# Look at the fraction of \omega
domain_breakdown_5[order(domain_breakdown_5$frac, decreasing=T)[1:10], ]
domain_breakdown_5[order(domain_breakdown_5$domain_omega, decreasing=F)[1:10], ]


# S-100
domain_s100 <- subset(slr_structure_all, domain_type == "PF01023")
nrow(domain_s100)
unique(domain_s100$stable_i)

domain_breakdown