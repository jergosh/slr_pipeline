

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
revigo.data <- rbind(c("GO:0001906","cell killing", 0.040,-0.402,-1.147, 3.876,-1.5213,0.994,0.000),
c("GO:0002480","antigen processing and presentation of exogenous peptide antigen via MHC class I, TAP-independent", 0.000,-3.752, 0.890, 1.763,-9.9269,0.756,1.000),
c("GO:0022617","extracellular matrix disassembly", 0.002, 1.817,-2.062, 2.542,-8.5743,0.849,0.000),
c("GO:0051707","response to other organism", 0.745,-7.552, 2.750, 5.145,-8.2828,0.739,1.0),
c("GO:0070989","oxidative demethylation", 0.001,-0.367,-0.982, 2.217,-5.5680,0.979,1.000),
c("GO:0009822","alkaloid catabolic process", 0.000, 0.454,-2.861, 1.544,-4.9052,0.964,1.0),
c("GO:0042107","cytokine metabolic process", 0.012,-0.184, 0.040, 3.363,-4.0781,0.951,1.0),
c("GO:1901687","glutathione derivative biosynthetic process", 0.000,-0.702,-1.441, 0.301,-1.8065,0.975,1.049),
c("GO:0044598","doxorubicin metabolic process", 0.000,-0.243,-1.674, 1.799,-4.7302,0.972,0.057),
c("GO:0016098","monoterpenoid metabolic process", 0.000, 1.313,-5.076, 1.863,-5.3831,0.935,1.0),
c("GO:0032078","negative regulation of endodeoxyribonuclease activity", 0.000, 1.600, 2.957, 1.342,-1.6336,0.848,0.061),
c("GO:0006690","icosanoid metabolic process", 0.012, 1.083,-4.647, 3.341,-4.1427,0.923,1.0),
c("GO:0043603","cellular amide metabolic process", 0.824,-0.362,-1.587, 5.189,-1.3465,0.968,1.0),
c("GO:0015721","bile acid and bile salt transport", 0.001, 0.753, 0.021, 2.121,-2.9858,0.884,0.089),
c("GO:0019748","secondary metabolic process", 0.144,-0.006,-2.786, 4.432,-2.9199,0.975,0.089),
c("GO:0044259","multicellular organismal macromolecule metabolic process", 0.005, 4.174, 6.605, 2.994,-7.9057,0.734,1.0),
c("GO:0006311","meiotic gene conversion", 0.000, 2.148, 0.603, 1.813,-2.7983,0.834,0.105),
c("GO:0038001","paracrine signaling", 0.000,-1.558, 0.126, 1.813,-1.3376,0.874,1.0),
c("GO:0001694","histamine biosynthetic process", 0.000, 0.090,-1.094, 1.255,-1.3376,0.966,1.0),
c("GO:0007155","cell adhesion", 0.622, 0.448, 2.450, 5.067,-3.5830,0.835,0.155),
c("GO:0035695","mitochondrion degradation by induced vacuole formation", 0.000, 0.352,-2.433, 1.431,-1.6336,0.970,0.188),
c("GO:0070383","DNA cytosine deamination", 0.000, 0.929,-1.378, 1.643,-1.6336,0.957,0.193),
c("GO:0019835","cytolysis", 0.080, 0.414, 1.522, 4.176,-3.0748,0.848,0.203),
c("GO:0000961","negative regulation of mitochondrial RNA catabolic process", 0.000, 2.425, 2.160, 0.903,-1.3376,0.806,0.204),
c("GO:0045234","protein palmitoleylation", 0.000, 0.271,-0.408, 1.380,-1.6336,0.941,0.212),
c("GO:0032385","positive regulation of intracellular cholesterol transport", 0.000, 1.574, 4.821, 1.301,-1.6336,0.730,0.216),
c("GO:1901373","lipid hydroperoxide transport", 0.000, 0.813,-0.119, 1.301,-1.6336,0.895,0.216),
c("GO:0051604","protein maturation", 0.265,-0.066, 0.262, 4.696,-2.0582,0.944,0.237),
c("GO:0006968","cellular defense response", 0.001,-6.466, 2.090, 2.410,-3.0915,0.813,0.250),
c("GO:0035606","peptidyl-cysteine S-trans-nitrosylation", 0.000, 0.217,-0.372, 1.398,-1.6336,0.946,0.251),
c("GO:0071260","cellular response to mechanical stimulus", 0.005,-6.691, 2.447, 2.992,-2.1914,0.746,0.273),
c("GO:0070208","protein heterotrimerization", 0.001, 2.214,-2.717, 2.362,-1.6690,0.956,0.274),
c("GO:0042448","progesterone metabolic process", 0.001, 2.973,-0.603, 2.228,-4.2691,0.790,0.279),
c("GO:0006749","glutathione metabolic process", 0.128, 1.309,-4.717, 4.379,-2.6058,0.933,0.291),
c("GO:0007229","integrin-mediated signaling pathway", 0.027,-3.960, 4.746, 3.710,-2.3090,0.694,0.308),
c("GO:0043062","extracellular structure organization", 0.024, 1.858, 0.082, 3.655,-3.3146,0.842,0.316),
c("GO:0071971","extracellular vesicular exosome assembly", 0.000, 1.569,-1.352, 1.447,-1.6336,0.863,0.319),
c("GO:0006805","xenobiotic metabolic process", 0.065,-6.324, 3.200, 4.087,-7.8362,0.700,0.330),
c("GO:0033559","unsaturated fatty acid metabolic process", 0.013, 1.564,-5.221, 3.376,-4.5610,0.916,0.332),
c("GO:0043252","sodium-independent organic anion transport", 0.000, 1.507,-0.643, 1.279,-1.5159,0.893,0.343),
c("GO:0072376","protein activation cascade", 0.124,-5.786, 2.786, 4.368,-5.1142,0.728,0.348),
c("GO:0035825","reciprocal DNA recombination", 0.006, 0.561,-0.982, 3.082,-1.3269,0.955,0.348),
c("GO:0009611","response to wounding", 0.223,-6.528, 2.419, 4.622,-3.0159,0.793,0.0),
c("GO:2000504","positive regulation of blood vessel remodeling", 0.000, 2.657, 7.311, 1.462,-1.3376,0.654,0.368),
c("GO:0010477","response to sulfur dioxide", 0.000,-6.621, 1.161, 1.176,-1.6336,0.852,0.368),
c("GO:0008228","opsonization", 0.000,-2.746, 1.448, 1.708,-3.3913,0.690,0.370),
c("GO:0035694","mitochondrial protein catabolic process", 0.000, 1.391,-1.966, 1.833,-1.3376,0.832,0.371),
c("GO:0046479","glycosphingolipid catabolic process", 0.004, 1.376,-5.623, 2.849,-1.6690,0.916,0.372),
c("GO:1900225","regulation of NLRP3 inflammasome complex assembly", 0.000, 3.332, 0.692, 0.699,-1.6336,0.849,0.374),
c("GO:0051606","detection of stimulus", 0.290,-6.953, 2.049, 4.735,-1.4731,0.821,0.375),
c("GO:0045071","negative regulation of viral genome replication", 0.003, 1.374, 4.819, 2.678,-3.3072,0.729,0.378),
c("GO:0090350","negative regulation of cellular organofluorine metabolic process", 0.000, 1.613, 3.595, 1.176,-1.6336,0.814,0.378),
c("GO:0002537","nitric oxide production involved in inflammatory response", 0.000,-3.174, 7.254, 1.690,-1.3376,0.648,0.382),
c("GO:0008049","male courtship behavior", 0.000,-2.161, 7.372, 1.973,-1.6336,0.685,0.386),
c("GO:0043589","skin morphogenesis", 0.001, 4.155, 6.509, 2.322,-1.3585,0.761,0.403),
c("GO:0042738","exogenous drug catabolic process", 0.001,-6.530, 0.051, 2.009,-4.2691,0.825,0.406),
c("GO:0010529","negative regulation of transposition", 0.000, 1.287, 4.019, 1.954,-1.6336,0.775,0.409),
c("GO:0035759","mesangial cell-matrix adhesion", 0.000,-0.279,-0.250, 1.491,-1.3376,0.887,0.411),
c("GO:0006821","chloride transport", 0.078, 1.583, 1.034, 4.165,-1.6619,0.860,0.413),
c("GO:0045019","negative regulation of nitric oxide biosynthetic process", 0.001, 2.109, 4.375, 2.328,-1.5885,0.793,0.425),
c("GO:0032787","monocarboxylic acid metabolic process", 1.736, 1.443,-5.406, 5.512,-1.3162,0.931,0.427),
c("GO:0009605","response to external stimulus", 1.251,-7.527, 2.223, 5.370,-3.3866,0.803,0.434),
c("GO:2000473","positive regulation of hematopoietic stem cell migration", 0.000, 1.033, 4.684, 1.531,-1.3376,0.730,0.439),
c("GO:0007339","binding of sperm to zona pellucida", 0.003, 3.941, 6.589, 2.757,-2.3125,0.716,0.443),
c("GO:0030844","positive regulation of intermediate filament depolymerization", 0.000, 2.340, 3.408, 0.778,-1.6336,0.728,0.450),
c("GO:0032534","regulation of microvillus assembly", 0.000, 2.382, 1.615, 1.398,-1.6336,0.807,0.453),
c("GO:0071460","cellular response to cell-matrix adhesion", 0.000,-6.018, 1.635, 1.447,-1.3376,0.784,0.458),
c("GO:0006948","induction by virus of host cell-cell fusion", 0.045, 0.371, 1.289, 3.929,-1.6336,0.886,0.464),
c("GO:0050909","sensory perception of taste", 0.021, 4.087, 6.924, 3.586,-7.5204,0.737,0.467),
c("GO:0071412","cellular response to genistein", 0.000,-6.717, 1.162, 0.778,-1.6336,0.778,0.483),
c("GO:0006720","isoprenoid metabolic process", 0.425, 1.448,-5.217, 4.901,-2.0841,0.915,0.494),
c("GO:0043900","regulation of multi-organism process", 0.096, 0.825, 4.792, 4.256,-1.3701,0.779,0.494),
c("GO:0001101","response to acid", 0.018,-7.210, 1.749, 3.526,-2.4735,0.821,0.495),
c("GO:0042574","retinal metabolic process", 0.001, 1.323,-5.102, 2.348,-1.6690,0.930,0.499),
c("GO:0002253","activation of immune response", 0.218,-3.959, 5.322, 4.611,-9.0202,0.430,0.0),
c("GO:1900035","negative regulation of cellular response to heat", 0.000,-4.659, 4.674, 1.477,-1.6336,0.660,0.502),
c("GO:0002412","antigen transcytosis by M cells in mucosal-associated lymphoid tissue", 0.000,-5.801, 2.335, 1.342,-1.6336,0.661,0.508),
c("GO:0050817","coagulation", 0.026, 3.902, 6.776, 3.691,-2.5424,0.737,1.0),
c("GO:0006636","unsaturated fatty acid biosynthetic process", 0.008, 1.561,-5.184, 3.198,-1.3881,0.901,0.516),
c("GO:0032930","positive regulation of superoxide anion generation", 0.001, 1.989, 5.099, 1.982,-1.9800,0.752,0.517),
c("GO:1901224","positive regulation of NIK/NF-kappaB cascade", 0.000,-3.971, 5.338, 0.602,-1.6336,0.677,0.525),
c("GO:0046600","negative regulation of centriole replication", 0.000, 2.658, 2.450, 1.398,-1.6336,0.746,0.527),
c("GO:0006922","cleavage of lamin involved in execution phase of apoptosis", 0.000, 1.541,-1.830, 1.839,-2.7983,0.821,0.528),
c("GO:0007162","negative regulation of cell adhesion", 0.009, 1.443, 4.231, 3.248,-1.8962,0.740,0.528),
c("GO:0055093","response to hyperoxia", 0.001,-7.209, 2.151, 2.352,-2.1060,0.837,0.529),
c("GO:0010572","positive regulation of platelet activation", 0.000,-1.968, 6.903, 1.973,-3.9203,0.473,0.534),
c("GO:0050665","hydrogen peroxide biosynthetic process", 0.001, 0.192,-0.973, 2.430,-1.6690,0.963,0.534),
c("GO:0002322","B cell proliferation involved in immune response", 0.000,-5.935, 2.731, 1.041,-1.3376,0.617,0.536),
c("GO:0030277","maintenance of gastrointestinal epithelium", 0.001, 3.355, 6.891, 2.270,-1.7593,0.704,0.537),
c("GO:0019732","antifungal humoral response", 0.000,-6.593, 3.151, 1.613,-1.6336,0.648,0.537),
c("GO:0016488","farnesol catabolic process", 0.000, 1.305,-4.956, 1.000,-1.3376,0.926,0.537),
c("GO:0009410","response to xenobiotic stimulus", 0.066,-7.011, 1.913, 4.094,-7.6909,0.808,0.539),
c("GO:0042035","regulation of cytokine biosynthetic process", 0.011, 3.274, 6.850, 3.300,-3.0676,0.596,0.540),
c("GO:0045575","basophil activation", 0.000,-3.276, 1.605, 1.176,-1.3376,0.682,0.543),
c("GO:0071374","cellular response to parathyroid hormone stimulus", 0.000,-6.975, 1.488, 1.833,-1.3376,0.760,0.546),
c("GO:0001676","long-chain fatty acid metabolic process", 0.025, 1.520,-5.150, 3.669,-2.7217,0.914,0.549),
c("GO:0034694","response to prostaglandin stimulus", 0.001,-7.008, 1.432, 2.326,-1.3341,0.826,0.551),
c("GO:0032504","multicellular organism reproduction", 0.084, 4.340, 6.283, 4.198,-1.7572,0.761,0.556),
c("GO:0006952","defense response", 0.961,-7.134, 2.477, 5.256,-2.3596,0.773,0.558),
c("GO:0052547","regulation of peptidase activity", 0.108, 1.965, 4.560, 4.304,-1.3380,0.816,0.561),
c("GO:0051817","modification of morphology or physiology of other organism involved in symbiotic interaction", 0.408, 1.004, 5.743, 4.884,-2.9510,0.735,0.563),
c("GO:2000617","positive regulation of histone H3-K9 acetylation", 0.000, 2.266, 3.221, 1.771,-1.6336,0.682,0.568),
c("GO:0010939","regulation of necrotic cell death", 0.002, 1.392, 3.155, 2.494,-1.3896,0.806,0.572),
c("GO:0042554","superoxide anion generation", 0.003, 0.351,-1.351, 2.715,-1.3390,0.968,0.576),
c("GO:0030887","positive regulation of myeloid dendritic cell activation", 0.000,-1.330, 4.945, 0.301,-1.6336,0.620,0.581),
c("GO:0001580","detection of chemical stimulus involved in sensory perception of bitter taste", 0.001,-3.148, 6.630, 2.045,-1.8617,0.668,0.585),
c("GO:0042737","drug catabolic process", 0.001, 0.182,-3.277, 2.086,-3.2799,0.967,0.597),
c("GO:0003014","renal system process", 0.008, 4.235, 6.830, 3.154,-1.6573,0.749,0.598),
c("GO:2000259","positive regulation of protein activation cascade", 0.000,-3.558, 5.558, 1.663,-1.3376,0.633,0.599),
c("GO:0071395","cellular response to jasmonic acid stimulus", 0.002,-6.799, 1.792, 2.605,-6.5431,0.739,0.600),
c("GO:2001199","negative regulation of dendritic cell differentiation", 0.000, 0.988, 7.269, 1.602,-4.3107,0.526,0.606),
c("GO:0042092","type 2 immune response", 0.003,-6.958, 2.667, 2.743,-1.6619,0.670,0.608),
c("GO:0042060","wound healing", 0.041,-7.751, 2.576, 3.884,-2.1395,0.778,0.609),
c("GO:0001816","cytokine production", 0.114, 3.644, 6.762, 4.330,-1.5727,0.717,0.0),
c("GO:0071657","positive regulation of granulocyte colony-stimulating factor production", 0.000, 2.620, 7.275, 1.663,-1.3376,0.630,0.617),
c("GO:0070887","cellular response to chemical stimulus", 0.430,-5.643, 3.607, 4.906,-3.5969,0.683,0.619),
c("GO:0070434","positive regulation of nucleotide-binding oligomerization domain containing 2 signaling pathway", 0.000,-4.266, 4.954, 1.477,-1.6336,0.530,0.623),
c("GO:0016045","detection of bacterium", 0.003,-6.877, 2.537, 2.719,-3.6229,0.784,0.627),
c("GO:0071230","cellular response to amino acid stimulus", 0.005,-7.071, 1.985, 2.974,-4.0092,0.728,0.632),
c("GO:0001869","negative regulation of complement activation, lectin pathway", 0.000,-3.704, 5.093, 1.633,-1.3376,0.489,0.635),
c("GO:0034341","response to interferon-gamma", 0.006,-7.374, 2.764, 3.072,-7.7300,0.598,0.638),
c("GO:2001179","regulation of interleukin-10 secretion", 0.001, 3.462, 6.750, 1.982,-2.2888,0.655,0.639),
c("GO:2000256","positive regulation of male germ cell proliferation", 0.000, 2.856, 7.480, 0.845,-1.6336,0.628,0.642),
c("GO:0050900","leukocyte migration", 0.022,-2.929, 2.311, 3.606,-2.7803,0.635,0.648),
c("GO:0002443","leukocyte mediated immunity", 0.022,-3.745, 1.039, 3.625,-8.8443,0.704,0.0),
c("GO:0071333","cellular response to glucose stimulus", 0.005,-4.722, 4.172, 2.968,-2.1914,0.671,0.655),
c("GO:0006957","complement activation, alternative pathway", 0.000,-3.960, 4.744, 1.908,-1.6690,0.514,0.0),
c("GO:0007597","blood coagulation, intrinsic pathway", 0.000,-2.059, 6.658, 1.643,-2.1060,0.607,0.659),
c("GO:0001971","negative regulation of activation of membrane attack complex", 0.000,-3.876, 5.373, 1.934,-1.6336,0.478,0.660),
c("GO:0002643","regulation of tolerance induction", 0.001, 0.780, 7.189, 2.294,-2.8787,0.520,0.0),
c("GO:0032645","regulation of granulocyte macrophage colony-stimulating factor production", 0.001, 3.296, 6.943, 2.310,-1.7593,0.669,0.663),
c("GO:0002291","T cell activation via T cell receptor contact with antigen bound to MHC molecule on antigen presenting cell", 0.000,-5.813, 2.407, 1.613,-1.3376,0.610,0.665),
c("GO:0002883","regulation of hypersensitivity", 0.001,-5.511, 4.626, 2.236,-1.9800,0.552,0.673),
c("GO:1900168","positive regulation of glial cell line-derived neurotrophic factor secretion", 0.000, 2.829, 7.198, 1.342,-1.6336,0.621,0.674),
c("GO:1900038","negative regulation of cellular response to hypoxia", 0.000,-5.153, 4.530, 1.477,-1.6336,0.636,0.676),
c("GO:0008206","bile acid metabolic process", 0.003, 1.530,-5.168, 2.698,-3.2841,0.913,0.683),
c("GO:0050850","positive regulation of calcium-mediated signaling", 0.003,-3.783, 5.656, 2.678,-1.8065,0.609,0.684),
c("GO:0033214","iron assimilation by chelation and transport", 0.110, 1.140, 4.331, 4.314,-1.3376,0.752,0.688),
c("GO:0006953","acute-phase response", 0.005,-6.654, 2.319, 2.985,-1.9451,0.778,0.689),
c("GO:0032673","regulation of interleukin-4 production", 0.002, 3.194, 6.946, 2.660,-2.9020,0.662,0.690),
c("GO:0043436","oxoacid metabolic process", 8.401, 1.419,-4.678, 6.197,-1.4214,0.927,0.0),
c("GO:0002250","adaptive immune response", 0.021,-6.575, 2.644, 3.604,-8.0174,0.638,0.691),
c("GO:0032689","negative regulation of interferon-gamma production", 0.003, 2.873, 6.895, 2.687,-1.9161,0.620,0.692),
c("GO:0044597","daunorubicin metabolic process", 0.000,-0.044,-2.152, 1.820,-4.7302,0.971,0.694),
c("GO:0032608","interferon-beta production", 0.003, 4.094, 6.360, 2.792,-1.3913,0.724,0.695),
c("GO:0010845","positive regulation of reciprocal meiotic recombination", 0.000, 2.016, 4.791, 1.623,-1.3376,0.688,0.696),
c("GO:0032613","interleukin-10 production", 0.003, 4.074, 6.364, 2.806,-2.2604,0.724,0.697),
c("GO:0043017","positive regulation of lymphotoxin A biosynthetic process", 0.000, 2.885, 7.015, 0.301,-1.6336,0.646,0.697),
c("GO:0006699","bile acid biosynthetic process", 0.001, 1.493,-5.478, 2.435,-1.6619,0.910,0.698));

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
pdf("revigo-plot.pdf", height=9, width=10)
p1;
dev.off()

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("revigo-plot.pdf");
