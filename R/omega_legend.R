colfunc <- colorRampPalette(c("red", "white"))
legend_image <- as.raster(matrix(colfunc(50), ncol=1))

pdf("omega_legend.pdf", height=8, width=4)
plot(c(-0.3,1),c(0, 1),type = 'n', axes = F,xlab = '', ylab = '', main=expression(omega), cex.main=2)
text(y=seq(0,1,l=5), x=-0.2, labels=seq(0, 1, l=5), cex=1.8)

rasterImage(legend_image, 0, 0, 1,1)
dev.off()

