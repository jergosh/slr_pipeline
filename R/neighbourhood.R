library(Hmisc)

# csa_sites <- subset(pdb_master_csa, pdb_id == "2b6h")
csa_sites <- pdb_master_csa
csa_sites$source <- csa_sites$csa == "Catalytic site"
write.table(file="data/catalytic_sites.tab", x=csa_sites, sep="\t", quote=F)

ggplot(pdb_master_csa, aes(x=omega, fill=csa)) +
  theme_minimal(base_size=10, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(0, 2)) +
  scale_y_continuous(limits=c(0, 5)) +
  xlab("omega") +
  ylab("Count") +  
  geom_density(alpha=0.4)

# Now that the Python script has run
csa_nbhood <- read.table("data/catalytic_site_dists.tab", sep="\t", stringsAsFactors=F, header=T)
csa_nbhood <- subset(csa_nbhood, buried == "exposed" & dist != 0)
csa_nbhood$dist_bin <- cut2(csa_nbhood$dist, cuts=seq(0, max(csa_nbhood$dist[!is.infinite(csa_nbhood$dist)]), 1))

csa_nbhood$dist_bin[as.numeric(csa_nbhood$dist_bin) > 5] <- "[ 40, 48)"
# Catalytic sites under pos. sel.


# How many very distant sites?
nrow(csa_nbhood)
nrow(subset(csa_nbhood, dist > 40))
nrow(subset(csa_nbhood, dist < 40))

# Fisher's exact test
nrow(subset(csa_nbhood, dist > 0 & dist < 4))
nrow(subset(csa_nbhood, dist > 0 & dist < 4 & omega > 1 & Adj.Pval < 0.05))

nrow(subset(csa_nbhood, omega > 1 & Adj.Pval < 0.05))
nrow(csa_nbhood)

fisher.test(matrix(c(), nrow=2))
# nrow(subset(csa_nbhood, omega > 1 & Adj.Pval < 0.05))

ggplot(csa_nbhood, aes(y=omega, x=dist)) + 
  theme_minimal(base_size=8, base_family="Helvetica") +
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  # scale_x_discrete(labels=bin_labels) + 
  # xlab("Relative solvent accessibility") +
  # ylab(expression(omega)) +
  coord_cartesian(ylim=c(0, 5)) +
  geom_point(size=0.2, alpha=0.05, colour="black") +
  geom_smooth(method = "gam", size = 1.5, method.args=list(bw=0.01))

dist_fractions <- make.fractions(csa_nbhood, c("dist_bin"), thr=0.05)

ggplot(dist_fractions, aes(x=dist_bin, y=Fraction)) +
  theme_minimal(base_size=8) +
  theme( # legend.position="none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab("Distance from catalytic site") +
  ylab("Fraction under positive selection") +
  # scale_fill_manual(values=c("dodgerblue2", "dodgerblue4"), guide=guide_legend(title="Solvent\naccessibility")) +
  # scale_x_discrete(labels=bin_labels) + 
  geom_bar(stat="identity")

plot(density(csa_nbhood$dist))
