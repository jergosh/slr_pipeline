setwd("~/Documents/projects/slr_pipeline/results/")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

joint.levels <- c("Overall", )
## First, all of globular
# Collect all data and make sure there is a common scale on all plots
# The overall proportion
# Buried vs. exposed
fractions.buried <- make.fractions(pdb_master_globular, "buried")
# Secondary structure
fractions.ss <- make.fractions(pdb_master_globular, "sec_simple")
# 'Everything'
fractions.all <- make.fractions(pdb_master_globular, c("buried", "sec_simple"))

ylim <- c(0, max(c(fraction.overall$Fraction, fractions.buried$Fraction, fractions.ss$Fraction, fractions.all$Fraction)))

fraction.overall <- data.frame(group=factor(c("All", "Buried", "Exposed"), levels=c("All", "Buried", "Exposed")),
                               Fraction=c(sum(pdb_master_globular$omega > 1.0) / nrow(pdb_master_globular),
                               fractions.buried$Fraction[1],
                               fractions.buried$Fraction[2]))
)
p.overall <- ggplot(fraction.overall, aes(x=group, y=Fraction, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) +
  scale_fill_manual(values=c("black", "indianred", "indianred2"), guide=F) + 
  scale_y_continuous(lim=ylim) +
  xlab("") + 
  theme_bw(base_size=18)

print(p.overall)
ggsave("TAC_fig2a.pdf", height=5, width=6)

fraction.ss <- data.frame(group=c(rep(c("All", "Buried", "Exposed"), each=4)),
                          secondary=factor(c(
                            as.character(fractions.ss$sec_simple)
                            ), levels=c("Helix", "Beta sheet", "Coil", "Turn")),
                          Fraction=c(
                            fractions.ss$Fraction,
                            fractions.all$Fraction
                          )
)

p.ss <- ggplot(fraction.ss, aes(x=group, y=Fraction, fill=secondary)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.8) +
  # scale_fill_manual(values=c("black", "indianred", "indianred2"), guide=F) + 
  scale_y_continuous(lim=ylim) +
  xlab("") + 
  theme_bw(base_size=18)

print(p.ss)

multiplot(p.overall, p.ss, cols=2)
ggsave("TAC_fig2b.pdf", height=5, width=10)

## Domains
fraction.domains <- make.fractions(pdb_master_globular, "is_domain")

fraction.domain.overall <- data.frame(group=factor(c("All", "Domain", "Non-domain"), levels=c("All", "Domain", "Non-domain")),
                               Fraction=c(sum(pdb_master_globular$omega > 1.0) / nrow(pdb_master_globular),
                                          fraction.domains$Fraction[1],
                                          fraction.domains$Fraction[2]))

p.domain.overall <- ggplot(fraction.domain.overall, aes(x=group, y=Fraction, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.5) +
  scale_fill_manual(values=c("black", "dodgerblue4", "dodgerblue3"), guide=F) + 
  scale_y_continuous(lim=ylim) +
  xlab("") + 
  ggtitle("Overall") +
  theme_bw(base_size=18)

print(p.domain.overall)
ggsave("TAC_fig3a.pdf", height=5, width=3.5)



fractions.domains <- make.fractions(pdb_master_globular, c("sec_simple", "is_domain"))
fractions.domains.exposed <- make.fractions(subset(pdb_master_globular, buried=="exposed"), c("sec_simple", "is_domain"))
fractions.domains.buried <- make.fractions(subset(pdb_master_globular, buried=="buried"), c("sec_simple", "is_domain"))

ylim.domains <- c(0, max(fraction.domain.overall$Fraction,
                         fractions.domains$Fraction,
                         fractions.domains.buried$Fraction,
                         fractions.domains.exposed$Fraction))

p <- ggplot(fractions.domains, aes(x=sec_simple, y=Fraction, fill=is_domain)) +
  #   scale_fill_manual(values=c("lightblue", "blue", "darkblue",
  #                              "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
  #                              "brown1", "brown3", "brown4")) +
  scale_fill_discrete(guide=F) +
  scale_y_continuous(lim=ylim.domains) + 
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Domains (all)") +
  xlab("") +
  theme_bw(base_size=18)

print(p)
ggsave("TAC_fig3b.pdf", height=5, width=6)

p <- ggplot(fractions.domains.buried, aes(x=sec_simple, y=Fraction, fill=is_domain)) +
  #   scale_fill_manual(values=c("lightblue", "blue", "darkblue",
  #                              "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
  #                              "brown1", "brown3", "brown4")) +
  scale_fill_discrete(guide=F) +
  scale_y_continuous(lim=ylim.domains) + 
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Domains (buried)") +
  xlab("") +
  ylab("") +
  theme_bw(base_size=18)

print(p)
ggsave("TAC_fig3c.pdf", height=5, width=6)

p <- ggplot(fractions.domains.exposed, aes(x=sec_simple, y=Fraction, fill=is_domain)) +
  #   scale_fill_manual(values=c("lightblue", "blue", "darkblue",
  #                              "darkolivegreen2", "darkolivegreen3", "darkolivegreen4",
  #                              "brown1", "brown3", "brown4")) +
  scale_fill_discrete() +
  scale_y_continuous(lim=ylim.domains) + 
  geom_bar(stat="identity", position=position_dodge()) +
  ggtitle("Domains (exposed)") +
  xlab("") +
  theme_bw(base_size=18)

print(p)
ggsave("TAC_fig3d.pdf", height=5, width=8)
