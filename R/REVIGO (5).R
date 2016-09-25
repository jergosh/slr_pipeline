

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0004872","receptor activity", 2.259, 0.930,-0.802, 5.649,-9.6356,0.982,0.000),
c("GO:0004888","transmembrane signaling receptor activity", 0.685,-1.348, 6.947, 5.131,-7.9842,0.885,0.000),
c("GO:0005201","extracellular matrix structural constituent", 0.015,-2.256,-1.144, 3.481,-7.2909,0.976,1.000),
c("GO:0008519","ammonium transmembrane transporter activity", 0.046,-4.923, 3.504, 3.957,-2.3032,0.951,1.000),
c("GO:0009055","electron carrier activity", 3.974, 0.393,-1.443, 5.894,-5.3426,0.982,1.000),
c("GO:0042605","peptide antigen binding", 0.001, 0.618,-1.679, 2.412,-7.9724,0.970,1.000),
c("GO:0061135","endopeptidase regulator activity", 0.123,-0.713, 0.291, 4.386,-2.9242,0.921,1.000),
c("GO:0070330","aromatase activity", 0.006, 7.238, 1.590, 3.106,-5.5433,0.780,0.000),
c("GO:0004364","glutathione transferase activity", 0.028,-2.887,-0.126, 3.746,-3.0131,0.958,0.014),
c("GO:0017171","serine hydrolase activity", 1.091,-1.425, 1.511, 5.332,-1.7062,0.954,1.020),
c("GO:0020037","heme binding", 2.465, 1.210,-1.625, 5.687,-5.9774,0.962,0.034),
c("GO:0048407","platelet-derived growth factor binding", 0.001,-1.582,-6.686, 2.465,-4.1093,0.913,0.034),
c("GO:0003823","antigen binding", 0.004, 0.437,-2.924, 2.885,-1.3159,0.970,0.037),
c("GO:0019825","oxygen binding", 0.069, 1.233,-3.127, 4.134,-3.6574,0.968,1.046),
c("GO:0033218","amide binding", 0.281, 0.426, 0.092, 4.744,-3.5164,0.967,1.053),
c("GO:0030246","carbohydrate binding", 0.811, 0.183,-0.315, 5.204,-4.4835,0.966,1.060),
c("GO:0005506","iron ion binding", 3.115,-0.310, 0.578, 5.788,-5.3788,0.959,0.071),
c("GO:0032052","bile acid binding", 0.000,-4.652, 1.862, 1.908,-1.9942,0.946,1.093),
c("GO:0036042","long-chain fatty acyl-CoA binding", 0.000,-1.496,-1.639, 1.279,-1.6410,0.964,1.099),
c("GO:0043199","sulfate binding", 0.000,-2.145,-1.901, 1.568,-1.6410,0.964,0.101),
c("GO:0047962","glycine N-benzoyltransferase activity", 0.000,-1.468, 0.755, 0.778,-1.6410,0.956,1.106),
c("GO:0047704","bile-salt sulfotransferase activity", 0.000,-0.121,-0.494, 0.845,-1.6410,0.963,1.107),
c("GO:0035605","peptidyl-cysteine S-nitrosylase activity", 0.000,-0.040, 0.151, 1.398,-1.6410,0.962,1.114),
c("GO:0035438","cyclic-di-GMP binding", 0.047,-2.774, 1.431, 3.967,-1.6410,0.946,1.133),
c("GO:0034875","caffeine oxidase activity", 0.000, 4.150,-0.546, 1.519,-4.3326,0.868,1.141),
c("GO:0047115","trans-1,2-dihydrobenzene-1,2-diol dehydrogenase activity", 0.001, 5.170,-0.549, 2.090,-4.3326,0.843,1.149),
c("GO:0001758","retinal dehydrogenase activity", 0.001, 4.997,-0.577, 2.155,-2.1337,0.855,0.150),
c("GO:0045550","geranylgeranyl reductase activity", 0.002, 5.839,-0.236, 2.628,-1.3449,0.855,0.159),
c("GO:0052692","raffinose alpha-galactosidase activity", 0.007,-3.083, 2.432, 3.166,-1.6410,0.954,0.166),
c("GO:0004033","aldo-keto reductase (NADP) activity", 0.006, 7.222,-2.433, 3.044,-3.5241,0.797,0.166),
c("GO:0004497","monooxygenase activity", 0.882, 8.114,-0.489, 5.240,-4.2517,0.808,0.225),
c("GO:0008526","phosphatidylinositol transporter activity", 0.000,-4.803, 3.448, 1.778,-1.6410,0.964,0.230),
c("GO:0004252","serine-type endopeptidase activity", 0.713,-3.011, 2.362, 5.148,-1.5596,0.955,0.245),
c("GO:0030197","extracellular matrix constituent, lubricant activity", 0.000,-1.524,-0.365, 1.000,-1.3449,0.976,0.254),
c("GO:0070975","FHA domain binding", 0.000,-2.940,-7.022, 1.362,-1.6410,0.927,0.272),
c("GO:0005139","interleukin-7 receptor binding", 0.000,-2.775,-7.130, 1.813,-1.6410,0.920,0.284),
c("GO:0033814","propanoyl-CoA C-acyltransferase activity", 0.001,-3.757, 1.057, 2.407,-1.6410,0.953,0.286),
c("GO:0019864","IgG binding", 0.000,-2.438,-6.702, 1.949,-3.0475,0.911,0.288),
c("GO:0019959","interleukin-8 binding", 0.001,-2.782,-7.002, 2.117,-1.3449,0.908,0.293),
c("GO:0016725","oxidoreductase activity", 0.214, 7.712,-0.497, 4.625,-2.8999,0.822,0.0),
c("GO:0045295","gamma-catenin binding", 0.001,-2.710,-6.876, 2.387,-1.5295,0.921,0.301),
c("GO:0047743","chlordecone reductase activity", 0.000, 6.062,-4.359, 0.602,-1.6410,0.835,0.311),
c("GO:0001848","complement binding", 0.003,-2.356,-6.593, 2.753,-1.8758,0.919,0.312),
c("GO:0047020","15-hydroxyprostaglandin-D dehydrogenase (NADP+) activity", 0.000, 6.497,-4.181, 0.845,-1.6410,0.832,0.320),
c("GO:0047718","indanol dehydrogenase activity", 0.000, 6.015,-3.834, 1.000,-2.8129,0.829,0.325),
c("GO:0036131","prostaglandin D2 11-ketoreductase activity", 0.000, 6.172,-3.444, 1.230,-1.6410,0.826,0.333),
c("GO:0047961","glycine N-acyltransferase activity", 0.001,-1.903, 0.797, 2.364,-1.3449,0.953,0.334),
c("GO:0036132","13-prostaglandin reductase activity", 0.000, 6.255,-3.352, 1.279,-1.3449,0.825,0.335),
c("GO:0047017","prostaglandin-F synthase activity", 0.000, 7.188,-3.765, 1.279,-1.6410,0.825,0.335),
c("GO:0015297","antiporter activity", 0.373,-4.341, 3.499, 4.866,-1.4864,0.952,0.338),
c("GO:0016627","oxidoreductase activity, acting on the CH-CH group of donors", 0.871, 8.281,-0.522, 5.235,-1.5400,0.808,0.341),
c("GO:0016705","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen", 0.904, 7.386,-0.557, 5.251,-1.9187,0.808,0.343),
c("GO:0016903","oxidoreductase activity, acting on the aldehyde or oxo group of donors", 0.887, 8.057,-0.524, 5.242,-1.3144,0.808,0.343),
c("GO:0016655","oxidoreductase activity, acting on NAD(P)H, quinone or similar compound as acceptor", 1.173, 7.942,-0.576, 5.364,-1.7187,0.805,0.353),
c("GO:0047042","androsterone dehydrogenase (B-specific) activity", 0.000, 7.699,-3.523, 1.857,-3.2834,0.776,0.355),
c("GO:0045703","ketoreductase activity", 0.001, 7.083,-2.900, 2.243,-1.6410,0.813,0.0),
c("GO:0004906","interferon-gamma receptor activity", 0.000,-2.007, 6.769, 0.477,-1.3449,0.916,0.373),
c("GO:0004032","alditol:NADP+ 1-oxidoreductase activity", 0.001, 6.918,-2.670, 2.438,-3.6485,0.808,0.378),
c("GO:0050649","testosterone 6-beta-hydroxylase activity", 0.000, 6.813, 3.201, 0.954,-1.6410,0.809,0.385),
c("GO:0016614","oxidoreductase activity, acting on CH-OH group of donors", 1.993, 8.765,-0.617, 5.594,-1.3681,0.799,0.387),
c("GO:0047086","ketosteroid monooxygenase activity", 0.000, 7.004, 3.179, 1.079,-4.9272,0.794,0.390),
c("GO:0019767","IgE receptor activity", 0.000,-2.116, 6.970, 1.519,-1.3449,0.909,0.430),
c("GO:0033038","bitter taste receptor activity", 0.000,-2.253, 7.045, 1.519,-2.8129,0.908,0.430),
c("GO:0001618","viral receptor activity", 0.001,-2.495, 6.900, 2.021,-2.6353,0.926,0.432),
c("GO:0047638","albendazole monooxygenase activity", 0.000, 5.192, 2.986, 0.301,-1.6410,0.808,0.446),
c("GO:0008395","steroid hydroxylase activity", 0.002, 7.306, 2.059, 2.486,-1.4687,0.792,0.453),
c("GO:0052869","arachidonic acid omega-hydroxylase activity", 0.002, 6.742, 1.712, 2.535,-2.8129,0.782,0.455),
c("GO:0004958","prostaglandin F receptor activity", 0.001,-2.587, 7.206, 2.000,-2.5185,0.896,0.458),
c("GO:0033780","taurochenodeoxycholate 6alpha-hydroxylase activity", 0.000, 5.679, 3.198, 0.602,-1.6410,0.802,0.462),
c("GO:0033783","25-hydroxycholesterol 7alpha-hydroxylase activity", 0.000, 6.194, 3.273, 0.699,-1.6410,0.800,0.466),
c("GO:0070538","oleic acid binding", 0.000,-1.462, 1.239, 1.342,-1.6410,0.949,0.467),
c("GO:0018676","(S)-limonene 7-monooxygenase activity", 0.000, 5.656, 2.703, 0.778,-3.2834,0.784,0.470),
c("GO:0047787","delta4-3-oxosteroid 5beta-reductase activity", 0.000, 4.813,-0.552, 1.447,-1.3449,0.850,0.471),
c("GO:0047522","15-oxoprostaglandin 13-oxidase activity", 0.000, 4.882,-0.548, 1.623,-1.3449,0.848,0.479),
c("GO:0018636","phenanthrene 9,10-monooxygenase activity", 0.000, 5.886, 2.705, 1.000,-4.9272,0.795,0.479),
c("GO:0015929","hexosaminidase activity", 0.040,-2.596, 2.184, 3.898,-1.4030,0.952,0.486),
c("GO:0004365","glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity", 0.032, 8.379,-0.214, 3.794,-1.3449,0.832,0.495),
c("GO:0050591","quinine 3-monooxygenase activity", 0.000, 6.939, 2.860, 1.415,-1.6410,0.788,0.497),
c("GO:0008329","pattern recognition receptor activity", 0.003,-1.313, 6.755, 2.790,-2.4475,0.909,0.513),
c("GO:0005496","steroid binding", 0.060,-1.937, 1.794, 4.071,-1.7980,0.947,0.0),
c("GO:0016709","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, NAD(P)H as one donor, and incorporation of one atom of oxygen", 0.095, 7.295, 1.041, 4.272,-1.7001,0.753,0.565),
c("GO:0032440","2-alkenal reductase [NAD(P)] activity", 0.012, 6.684,-0.483, 3.387,-1.6410,0.824,0.573),
c("GO:0016503","pheromone receptor activity", 0.017,-2.455, 7.137, 3.537,-1.6410,0.890,0.577),
c("GO:0010576","metalloenzyme regulator activity", 0.004,-1.218, 0.289, 2.851,-1.4635,0.940,0.583),
c("GO:0015116","sulfate transmembrane transporter activity", 0.137,-1.1, 3.467, 4.431,-1.5295,0.952,0.0),
c("GO:0042608","T cell receptor binding", 0.000,-2.759,-7.370, 1.699,-2.3032,0.910,0.598),
c("GO:0035410","dihydrotestosterone 17-beta-dehydrogenase activity", 0.000, 5.998,-3.580, 0.903,-1.6410,0.797,0.624),
c("GO:0047045","testosterone 17-beta-dehydrogenase (NADP+) activity", 0.000, 7.636,-3.783, 1.568,-1.3449,0.786,0.664),
c("GO:0061507","cyclic-GMP-AMP binding", 0.000,-1.144, 1.214, 0.845,-1.6410,0.959,0.668),
c("GO:0047035","testosterone dehydrogenase (NAD+) activity", 0.000, 7.662,-3.709, 1.643,-2.3032,0.784,0.668),
c("GO:0043120","tumor necrosis factor binding", 0.000,-2.564,-6.637, 1.380,-1.3449,0.912,0.670),
c("GO:0015299","solute:hydrogen antiporter activity", 0.169,-4.428, 3.522, 4.522,-1.3473,0.949,0.679),
c("GO:0047006","17-alpha,20-alpha-dihydroxypregn-4-en-3-one dehydrogenase activity", 0.000, 7.446,-3.453, 1.833,-3.2834,0.781,0.680),
c("GO:0005178","integrin binding", 0.006,-2.731,-7.408, 3.085,-2.0041,0.902,0.685));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 6 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
# p1 <- p1 + theme(legend.key = theme(text = element_text(size=20))) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen
pdf("revigo-plot-MF.pdf", height=9, width=10)
p1;
dev.off()
