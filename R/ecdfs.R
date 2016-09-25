head(pdb_master_globular)

ggplot(subset(pdb_master_globular, omega < 4), aes(x=omega, colour=sec_simple, linetype=buried)) + 
  stat_ecdf() +
  coord_cartesian(ylim=c(0.5, 1)) + 
  # geom_rect(aes(ymin=0, ymax=1, xmin=1, xmax=4), colour=rgb(1, 0, 0, 0.15)) + 
  ggtitle("Empirical CDF for different structural regions") + 
  ylab("eCDF") +
  xlab(expression(omega)) +
  theme_bw() + 
  theme(axis.text=element_text(size=16),
        plot.title=element_text(size=18,face="bold"),
        legend.text=element_text(size=16),
        legend.title=element_blank())
ggsave(paste(results, "eCDFs_structure.pdf", sep="/"), width=11, height=8)


