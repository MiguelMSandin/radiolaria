#----
#---- Libraries ------------------------------------------------------------------------------------

library(MCMCtreeR)
library(ape)
library(treeio)
library(HDInterval)
library(data.table)
library(dplyr)
library(ggplot2)

#----
setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/molClock_filtered/mcmcTree/MC07/")
#---- Setting names --------------------------------------------------------------------------------

files = list(fasta = "../all_filtered_align-linsi_trim05.fasta",
			 tree = "../all_filtered_align-linsi_trim05_gamma_BS103.RAxML_bipartitions_subset_annotated.tre",
			 calibrations = "../fossil_calibrations_v6.tsv",
			 treeOut = "all_filtered_align-linsi_trim05_gamma_BS103.RAxML_bipartitions_dated_MC07.newick",
			 treeAnnot = "all_filtered_align-linsi_trim05_gamma_BS103.RAxML_bipartitions_dated_MC07_annot.tre",
			 timeUnit = 100,
			 minTips = 0)

#----
#---- Remove branch lengths and any other node information -----------------------------------------

# Read tree
if(class(try(read.tree(files$tree), silent=TRUE))=="try-error"){
    tree <- read.nexus(files$tree) # For NEXUS
}else{
    tree <- read.tree(files$tree)  # For NEWICK
}

# Remove branch lengths
treeOut <- tree
treeOut$edge.length <- NULL
treeOut$node.label <- NULL

# Write tree
write.tree(treeOut, files$treeOut)

#_______________________________________#
# Now add the annotations from the table
#_______________________________________#


#----
#---- Read files -----------------------------------------------------------------------------------

# tree <- read.nexus(file_tree)

# Read annotated tree file _________________________________________________________________________
treei <- read.beast(files$treeAnnot)
treed <- as_tibble(treei)
tree_raw <- as.phylo(treei)

# Remove branch lengths
tree <- tree_raw
tree$edge.length <- NULL


# plot(tree)
# nodelabels(frame="none")

# stree <- extract.clade(tree, 701, root.edge=0)
# plot(stree, cex=.6)
# nodelabels(frame = "none")
# axisPhylo()


# Extracting annotations ___________________________________________________________________________
names(treed) <- gsub("!", "", names(treed))

groups <- as.vector(treed$name[!is.na(treed$name)])

# check for duplicated group names
if(any(duplicated(groups))){cat("The following annotations are duplicated:\n"); groups[duplicated(groups)]}else{cat("  No duplicates in annotated groups\n")}

# Linking the annotations to the nodes
nodes <- list()
for(g in groups){
    nodes[g] <- treed$node[grep(paste0("^", g, "$"), treed$name)]
}; rm(g)


# Read calibration file ____________________________________________________________________________

calib_raw <- fread(files$calibrations)
calib <- calib_raw[used=="yes"]

# Check all calibrations are annotated in the tree
if(!all(unique(calib$Node) %in% names(nodes))){
    cat("The following nodes are calibrated but not found in the tree:\n", unique(calib$Node)[!unique(calib$Node) %in% names(nodes)])
}
# And check all annotated nodes are calibrated
if(!all(names(nodes) %in% unique(calib$Node))){
    cat("The following nodes are annotated but there is no calibration:\n", names(nodes)[!names(nodes) %in% unique(calib$Node)])
}



#----
#---- Read calibrations file -----------------------------------------------------------------------

# Get calibrations from file _______________________________________________________________________

monophyleticGroups <- tipDes(tree, as.numeric(nodes))

calibrations <- list(node=c(), nodeName=c(), method=c(), min=c(), max=c())
for(i in names(nodes)){
    ss <- subset(calib, Node==i)
    calibrations$node <- cbind(calibrations$node, nodes[[i]])
    calibrations$nodeName <- cbind(calibrations$nodeName, ss$Node)
    calibrations$method <- cbind(calibrations$method, ss$Method)
    calibrations$min <- cbind(calibrations$min, ss$min/files$timeUnit)
    calibrations$max <- cbind(calibrations$max, ss$max/files$timeUnit)
}; rm(i, ss)


for(i in 1:length(calibrations[[1]])){
    if(i==1){
        sep1=max(nchar(calibrations$nodeName))
        sep2=max(nchar(calibrations$min, keepNA=FALSE))
        sep3=max(nchar(calibrations$max, keepNA=FALSE))
        cat("Node#",
            "\tNode", paste(rep(" ", sep1-4), collapse=""), 
            "\tmin", paste(rep(" ", sep2-3), collapse=""), 
            "\tmax", paste(rep(" ", sep3-3), collapse=""), "\tmethod\n")}
    cat(calibrations$node[i], "\t",
        calibrations$nodeName[i], paste(rep(" ", sep1-nchar(calibrations$nodeName[i], keepNA=FALSE)), collapse=""), "\t", 
        calibrations$min[i], paste(rep(" ", sep2-nchar(calibrations$min[i], keepNA=FALSE)), collapse=""), "\t", 
        calibrations$max[i], paste(rep(" ", sep3-nchar(calibrations$max[i], keepNA=FALSE)), collapse=""), "\t", 
        calibrations$method[i], "\n")
}; rm(i, sep1, sep2, sep3)

if(any(is.na(calibrations$min)) | any(is.na(calibrations$max)) | any(is.na(calibrations$method))){
    cat("\nWarning!!\nThere are some calibrations that are actually not calibrated...\nplease check and remove if necessary.\n\n")
}

for(i in 1:length(calibrations[[1]])){
    if(is.na(calibrations$min[i]) | is.na(calibrations$max[i]) | is.na(calibrations$method[i])){
        calibrations$node <- calibrations$node[-i]
        calibrations$nodeName <- calibrations$nodeName[-i]
        calibrations$method <- calibrations$method[-i]
        calibrations$min <- calibrations$min[-i]
        calibrations$max <- calibrations$max[-i]
    }
}; rm(i)
for(i in names(monophyleticGroups)){
    if(!i %in% calibrations$node){
        monophyleticGroups[i] <- NULL
    }
};rm(i)

cat("\nThere are", as.numeric(unique(lapply(calibrations, length))), "unique calibrated nodes\n\n")

#----
#---- Calibrate tree with different calibrations ---------------------------------------------------

# Export tree for MCMCTree _________________________________________________________________________
# output <- MCMCtreePhy(phy = tree, 
#                       minAge = calibrations$min,
#                       maxAge = calibrations$max, 
#                       monoGroups = monophyleticGroups, 
#                       method = calibrations$method,
#                       # shape = 10, scale = 1.5, 
#                       # estimateMode = TRUE, estimateShape = TRUE, estimateBeta = TRUE,
#                       plot = TRUE, pdfOutput = gsub("\\.[^\\.]+$", ".pdf", files$treeOut),
#                       writeMCMCtree = TRUE, MCMCtreeName = files$treeOut)


#---- Calibrate tree with skew-Normal distribution -------------------------------------------------

# Export tree for MCMCTree _________________________________________________________________________

output <- estimateSkewNormal(phy = tree, 
                      minAge = calibrations$min,
                      maxAge = calibrations$max,
                      monoGroups = monophyleticGroups,
                      writeMCMCtree = TRUE, MCMCtreeName = files$treeOut)

pdf(gsub("\\.[^\\.]+$", ".pdf", files$treeOut), width=11.69, height=8.27, paper='special')
for (i in 1:length(calibrations$node)) {
    plotMCMCtree(output$parameters[i,], 
                 method="skewNormal", 
                 title=paste0(calibrations$nodeName[i], " -node:", calibrations$node[i], 
                              "- min:", calibrations$min[i],
                              ", max:", calibrations$max[i],
                              " (location:", output$parameters[i,1], 
                              ", scale:", output$parameters[i,2], 
                              ", shape:", output$parameters[i,1], ")"),
                 lowerTime=round(calibrations$min[i]*0.95,2),
                 upperTime=round(calibrations$max[i]+calibrations$max[i]*0.2,2))
}; rm(i)
dev.off()

# skewNormal_results <- estimateSkewNormal(minAge = minimumTimes, maxAge = maximumTimes, monoGroups = monophyleticGroups, addMode = 0.05, phy = tree, plot = FALSE)


#----
#---- Explore MCMC output files --------------------------------------------------------------------

files$mcmcs = grep("^step1.*\\.txt$", dir(recursive=TRUE), value=TRUE)

mcmcs = data.frame()
for(file in files$mcmcs){
	cat("\r  Reading (", grep(file, files$mcmcs), "/", length(files$mcmcs), ") ", file, sep="", end="")
	tmp = fread(file)
	tmp$file = file
	mcmcs = rbind(mcmcs, tmp)
}; rm(file, tmp)

summary(mcmcs$lnL)
summary(mcmcs$t_n685)

ggplot(subset(mcmcs, Gen < max(mcmcs$Gen)*0.2), aes(x=Gen, y=lnL, colour=file))+
	geom_line()+
	theme(legend.position= "none")+
	theme_bw()

ggplot(mcmcs, aes(x=file, y=lnL))+
	# geom_violin()+
	geom_boxplot()+
	theme_bw()

ggplot(mcmcs, aes(x=file, y=t_n685))+
	geom_violin()+
	# geom_boxplot()+
	theme_bw()

datas = mcmcs %>% group_by(file) %>% summarise(lnLmean=mean(lnL), lnLsd=sd(lnL), rootMean=mean(t_n685), rootSD=sd(t_n685))
datas$lnLrank = rank(-datas$lnLmean)
datas$rootRank = rank(datas$rootMean)

ggplot(datas, aes(x=file, y=lnLmean))+
	geom_segment(aes(x=file, xend=file, y=lnLmean-lnLsd/2, yend=lnLmean+lnLsd/2), 
				 color="grey90", linewidth=2, lineend="round")+
	geom_point()+
	theme_bw()

ggplot(datas, aes(x=file, y=rootMean))+
	geom_segment(aes(x=file, xend=file, y=rootMean-rootSD/2, yend=rootMean+rootSD/2), 
				 color="grey90", linewidth=2, lineend="round")+
	geom_point()+
	theme_bw()

#----
#---- Plot trees -----------------------------------------------------------------------------------

geo <- data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian", "Ectasian", "Stenian", "Tonian", "Cryogenian", "Ediacran", "Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian","Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"),
				  era=c("Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Paleoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Mesoproterozoic", "Neoproterozoic", "Neoproterozoic", "Neoproterozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Paleozoic", "Mesozoic", "Mesozoic", "Mesozoic", "Cenozoic", "Cenozoic", "Cenozoic"))
geo$mid <- apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

phy = readMCMCtree("all_filtered_align-linsi_trim05_MC7_mcmc.tre")
priorTree=readMCMCtree("step-prior/all_filtered_align-linsi_trim05_MC7_step-prior.tre")
prior=fread("step-prior/all_filtered_align-linsi_trim05_MC7_step-prior.txt")
posterior = fread("all_filtered_align-linsi_trim05_MC7_mcmcChain_combined.txt.gz")

MCMC.tree.plot(phy, cex.tips = 0.2,
			   time.correction = 100,
			   scale.res = c("Eon", "Period"),
			   plot.type = "phylogram", cex.age = 0.6, cex.labels = 0.6,
			   relative.height = 0.08, col.tree = "grey40", label.offset = 4,
			   node.method = "none", no.margin = TRUE)


# pdf("figTree_figure_distributions.pdf", width=11.69, height=8.27, paper='special')
# MCMC.tree.plot(phy, MCMC.chain = posterior, cex.tips = 0.2, 
#  			   time.correction = 100, plot.type = "distributions", cex.age = 0.4, 
#  			   cex.labels = 0.5, relative.height = 0.08, col.tree = "grey40", 
#  			   scale.res = c("Eon", "Period"), no.margin = TRUE, label.offset = 4, 
#  			   density.col = "#00000050", density.border.col = "#00000080")
# dev.off()

pdf("figTree_figure.pdf", width=11.69, height=8.27, paper='special')
MCMC.tree.plot(phy, cex.tips = 0.2,
			   time.correction = 100, plot.type = "phylogram", lwd.bar = 2,
			   scale.res = c("Eon", "Period"), node.method = "bar", col.age = "lightblue4",
			   no.margin = TRUE, label.offset = 4)
dev.off()

#---- Plot distributions in specific nodes ---------------------------------------------------------

nodes = list(rads=713, 
			 polycystines=1026, Spumellaria=1181, SpumSymb=1329, 
			 Nassellaria=1027, Collodaria=1111, 
			 AcanRads=714, RadsBC=715, RadB=734, RadC=716,
			 AcanRadA=833, AcanthariaSensuLato=886, AcanSymb=959, RadA=834)

file = "step2_25burnin/all_filtered_align-linsi_trim05_MC7_mcmcChain_combined.txt.gz"
burnin = 0
factor = 100
data = fread(file)

hist(data$lnL)
# plot(data$Gen, data$lnL)
data = data[c(round(nrow(data)*burnin/100):nrow(data)),]

datas = select(data, paste0("t_n", nodes))

colnames(datas) = names(nodes)

# pdf(sub("\\.txt|\\.txt\\.gz", "_distributions.pdf", file), width=11.69, height=8.27, paper='special')
# for(n in names(datas)){
# 	ss = select(datas, n)
# 	colnames(ss) = "x"
# 	plotDist = ggplot(ss)+
# 		geom_density(aes(x=-x), fill="grey80")+
# 		labs(title=paste0(n))+
# 		theme_minimal()
# 	plot(plotDist)
# }; rm(n, ss)
# dev.off()

datam = melt(datas)

datamc = data.frame()
cat("  clade\tmedian\t95%HPD")
for(clade in unique(datam$variable)){
	ss = subset(datam, variable==clade)
	hpd = hdi(ss$value)
	cat("  ", clade, "\t", round(quantile(ss$value, 0.5)*factor), "\t(", round(hpd[1]*factor), "-", round(hpd[2]*factor), ")\n", sep="")
	ss = subset(ss, value >= hpd[1] & value <= hpd[2])
	datamc = rbind(datamc, ss)
}; rm(clade, ss, hpd)

(plotDist = ggplot(datamc)+
		geom_density(aes(x=-value), fill="grey80")+
		facet_wrap(~variable, scales="free")+
		theme_minimal())

pdf(sub("\\.txt|\\.txt\\.gz", "_distributions_HPD.pdf", file), width=11.69, height=8.27, paper='special')
plot(plotDist)
dev.off()

#

#---- Check prior Vs Posterior distributions -------------------------------------------------------

# phy = readMCMCtree("all_filtered_align-linsi_trim05_MC7_mcmc.tre")
# priorTree=readMCMCtree("step-prior/all_filtered_align-linsi_trim05_MC7_step-prior.tre")
prior=fread("step-prior/all_filtered_align-linsi_trim05_MC7_step-prior.txt")
posterior = fread("step2_25burnin/all_filtered_align-linsi_trim05_MC7_mcmcChain_combined.txt.gz")

calib = priorPosterior(MCMCPrior=prior, 
					   MCMCPosterior=posterior, 
					   inputTree="all_filtered_align-linsi_trim05_gamma_BS103.RAxML_bipartitions_dated_MC07.newick")

pr = as.data.frame(calib$prior); pr$dist = "prior"
po = as.data.frame(calib$posterior); po$dist = "posterior"
if(all(names(pr) == names(po))){prpo = rbind(pr, po); rm(pr, po)}else{cat("  Warning! Check columns match in prior and posterior files!\n")}
prpo = melt(as.data.table(prpo), id.vars=c("dist"))

prpo$value = prpo$value * 100
{prpo$calibration = as.character(prpo$variable)
	prpo$calibration[grep("1181", prpo$calibration)] = "Spumellaria (566.5-509.9)"
	prpo$calibration[grep("1318", prpo$calibration)] = "Spongopyle (37.8-33.9)"
	prpo$calibration[grep("1313", prpo$calibration)] = "Calcaromma (1-0)"
	prpo$calibration[grep("1322", prpo$calibration)] = "Hexastyloidea-group (266.2-239.6)"
	prpo$calibration[grep("1330", prpo$calibration)] = "Hollandosphaeridae (1-0)"
	prpo$calibration[grep("1183", prpo$calibration)] = "Cladococcus (64.7-59.2)"
	prpo$calibration[grep("1262", prpo$calibration)] = "Rhizosphaeroidea (162.8-146.5)"
	prpo$calibration[grep("1278", prpo$calibration)] = "Excentroconchidae (22.7-22.6)"
	prpo$calibration[grep("1338", prpo$calibration)] = "Dictyocoryne (47.8-41.2)"
	prpo$calibration[grep("1333", prpo$calibration)] = "Spongosphaera (11.6-7.24)"
	prpo$calibration[grep("1290", prpo$calibration)] = "Clade-M (106.7-89.8)"
	prpo$calibration[grep("1292", prpo$calibration)] = "Tetrapyle (5.3-3.6)"
	prpo$calibration[grep("1340", prpo$calibration)] = "Lithocyclioidea (49.5-44.6)"
	prpo$calibration[grep("1027", prpo$calibration)] = "Nassellaria (539.5-380)"
	prpo$calibration[grep("1099", prpo$calibration)] = "Plectopyramidoidea (269.5-242.6)"
	prpo$calibration[grep("1057", prpo$calibration)] = "Plagiacanthoidea (250-172)"
	prpo$calibration[grep("1059", prpo$calibration)] = "Lipmanella (47.8-41.2)"
	prpo$calibration[grep("1079", prpo$calibration)] = "Peromelissa (41.2-37.8)"
	prpo$calibration[grep("1089", prpo$calibration)] = "Artostrobioidea (200-170)"
	prpo$calibration[grep("1082", prpo$calibration)] = "Acanthodesmioidea (72.6-60)"
	prpo$calibration[grep("1044", prpo$calibration)] = "Pterocorythoidea (65.3-56.4)"
	prpo$calibration[grep("1138", prpo$calibration)] = "Collosphaeridae (36-26.9)"
	prpo$calibration[grep("1148", prpo$calibration)] = "Siphonosphaera (13.82-11.63)"
	prpo$calibration[grep("1112", prpo$calibration)] = "Sphaerozoidae (23-18.5)"
	prpo$calibration[grep("1173", prpo$calibration)] = "Collophidium (1-0)"}

unique(prpo$calibration)

tmp = as.data.frame(prpo %>% group_by(calibration) %>% summarise(min=min(value)))
prpo$calibration = factor(prpo$calibration, levels=c(tmp$calibration[order(tmp$min, decreasing = TRUE)])); rm(tmp)

(prpo_plot = ggplot(prpo, aes(x = -value, colour = dist, fill = dist))+
		# geom_histogram(alpha=0.4)+
		geom_density(alpha=0.4)+
		facet_wrap(~calibration, scales="free")+
		scale_colour_manual(values=c("#648FFF", "#FFB000"))+
		scale_fill_manual(values=c("#648FFF", "#FFB000"))+
		theme_bw())

pdf("prior_posterior.pdf", width=11.69*1.2, height=8.27*1.2, paper='special')
plot(prpo_plot)
dev.off()

#----
#---- LTT ------------------------------------------------------------------------------------------

data = fread("all_filtered_align-linsi_trim05_MC7_mcmc_LTT.tsv")

data$time = data$time * 100
data$hpd05 = data$hpd05 * 100
data$hpd95 = data$hpd95 * 100

unique(data$tree)
datas = subset(data, tree=="main" | 
			   	tree=="RAD-C" | 
			   	tree=="RAD-B" | 
			   	tree=="RAD-A" | 
			   	tree=="AcanthareaSL" | 
			   	tree=="Acantharea" | 
			   	tree=="Spumellaria" | 
			   	tree=="Entactinaria" | 
			   	tree=="Nassellaria" | 
			   	tree=="Collodaria")
datas$tree = factor(datas$tree, levels=c("main", "RAD-A", "RAD-B", "RAD-C", "AcanthareaSL", "Acantharea", 
										 "Spumellaria", "Entactinaria",  "Nassellaria", "Collodaria"))
datas$hpd05[which(datas$tree != "main")] = NA
datas$hpd95[which(datas$tree != "main")] = NA

geos = subset(geo, time > -max(datas$hpd95, na.rm=TRUE))

(lttplot = ggplot(datas)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_segment(aes(x=-hpd95, xend=-hpd05, y=lineages, yend=lineages), 
					 color="grey90", linewidth=2, lineend="round")+
		geom_line(aes(x=-time, y=lineages, color=tree)) +
		geom_text(data=geos, aes(x=mid, y=max(datas$lineages)+10, label=era), size=2.5)+
		geom_text(data=geos, aes(x=mid, y=max(datas$lineages)+5, label=period), size=2)+
		scale_x_continuous(breaks=seq((round(min(-datas$hpd95, na.rm=TRUE), -2)), 0, 100), 
						   minor_breaks=seq((round(min(-datas$hpd95, na.rm=TRUE), -2)), 0, 50)) +
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_color_manual(values=c("black", "grey80", "purple", "grey30", "yellow4", "yellow2", 
									"steelblue3", "steelblue4", "springgreen3", "orangered3"))+
		theme_classic()+
		labs(y="Ln (lineages)", x="Time (Ma)"))

pdf("all_filtered_align-linsi_trim05_MC7_mcmc_LTT.pdf", width=11.69, height=8.27, paper='special'); plot(lttplot); dev.off()


datass = subset(data, tree=="RAD-C" | tree=="RAD-B" | tree=="RAD-A" | tree=="AcanthareaSL" | tree=="Acantharea" | 
					tree=="Spumellaria" | tree=="Entactinaria" | tree=="Nassellaria" | tree=="Collodaria")
datass$tree = factor(datas$tree, levels=c("main", "RAD-A", "RAD-B", "RAD-C", "AcanthareaSL", "Acantharea", 
										  "Spumellaria", "Entactinaria",  "Nassellaria", "Collodaria"))

# geos = subset(geo, time > -max(datass$hpd95, na.rm=TRUE))
geos = subset(geo, time > -max(datass$time, na.rm=TRUE))

(lttplots = ggplot(datass)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		# geom_ribbon(aes(xmin=-hpd95, xmax=-hpd05, y=lineages, fill=tree), alpha=0.1, colour = NA) +
		geom_line(aes(x=-time, y=lineages, color=tree)) +
		geom_text(data=geos, aes(x=mid, y=max(datass$lineages)+10, label=era), size=2.5)+
		geom_text(data=geos, aes(x=mid, y=max(datass$lineages)+5, label=period), size=2)+
		# scale_x_continuous(breaks=seq((round(min(-datass$hpd95, na.rm=TRUE), -2)), 0, 100), 
		# 				   minor_breaks=seq((round(min(-datass$hpd95, na.rm=TRUE), -2)), 0, 50)) +
		scale_x_continuous(breaks=seq((round(min(-datass$time, na.rm=TRUE), -2)), 0, 100), 
						   minor_breaks=seq((round(min(-datass$time, na.rm=TRUE), -2)), 0, 50)) +
		scale_y_log10() + annotation_logticks(sides = 'l')+
		scale_color_manual(values=c("grey80", "purple", "grey30", "yellow4", "yellow2", 
									"steelblue3", "steelblue4", "springgreen3", "orangered3"))+
		# scale_fill_manual(values=c("grey80", "purple", "grey30", "yellow4", "yellow2",
		# 							"steelblue3", "steelblue4", "springgreen3", "orangered3"))+
		theme_classic()+
		labs(y="Ln (lineages)", x="Time (Ma)"))

pdf("all_filtered_align-linsi_trim05_MC7_mcmc_LTT_groups.pdf", width=11.69, height=8.27, paper='special'); plot(lttplots); dev.off()



#----
#----
#----
#---- Explore calibrations from BEAST --------------------------------------------------------------

n <- 10000

pdf(gsub("\\.[^\\.]+$", "_BEASTcalibrations.pdf", files$treeOut), width=11.69, height=8.27, paper='special')
for(node in calib$Node){
    ss <- subset(calib, Node==node)
    mean <- as.numeric(ss$Calibration %>% sub(".*\\(", "", .) %>% sub(",.*", "", .))
    sd <- as.numeric(ss$Calibration %>% sub(".*,", "", .) %>% sub("\\)", "", .))
    
    df <- data.frame(i=c(1:n), N=rnorm(n, mean, sd))
    q <- list(q01=quantile(df$N, 0.01),
              q05=quantile(df$N, 0.05),
              q95=quantile(df$N, 0.95),
              q99=quantile(df$N, 0.99))
    q <- lapply(q, function(x) round(x, 1))
    
    plot <- ggplot(df, aes(x=N)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="steelblue")+
        geom_vline(aes(xintercept=q$q01), color="steelblue", linetype="dashed", size=1)+
        geom_vline(aes(xintercept=q$q05), color="steelblue", size=1)+
        geom_vline(aes(xintercept=q$q95), color="steelblue", size=1)+
        geom_vline(aes(xintercept=q$q99), color="steelblue", linetype="dashed", size=1)+
        labs(title= paste(node, ": ", ss$Calibration, sep=""),
             subtitle=paste("Quantiles from ", n, " iterations: (", q$q01, " - (", q$q05, " - ", q$q95, ") - ", q$q99, ") Ma", sep=""),
             y="Density", x = "Age (Ma)")+
        theme_bw()
    plot(plot)
}; rm(n, node, ss, mean, sd, df, q, plot)
dev.off()


#----
#----
#----
#
