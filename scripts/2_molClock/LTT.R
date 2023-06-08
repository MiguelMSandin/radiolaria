#---- 
#---- loading packages ----

library(ape)
library(treeio)

library(BAMMtools)

library(data.table)
library(tidyr)

library(ggplot2)
library(ggtree)

# library(devtools)
# install_github("hmorlon/PANDA", dependencies = TRUE)
library(RPANDA)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, 
                          -1600 , -1400, -1200, 
                          -1000, -720, -635, 
                          -541, -485.4, -443.8, -419, -358.9, -298.9, 
                          -251.9, -201.4, -145, 
                          -66, -23.03, -2.58),
                  period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", 
                           "Calymnian","Ectasian","Stenian",
                           "Tonian","Cryogenian","Ediacran",
                           "Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian",
                           "Triassic","Jurassic","Cretaceous",
                           "Paleogene","Neogene","Quaternary"),
                  era=c("Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Paleoproterozoic",
                        "Mesoproterozoic","Mesoproterozoic","Mesoproterozoic",
                        "Neoproterozoic","Neoproterozoic","Neoproterozoic",
                        "Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic",
                        "Mesozoic","Mesozoic","Mesozoic",
                        "Cenozoic","Cenozoic","Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----  
setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/molClock_filtered/mcmcTree/MC07/step2_25burnin/")
setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/molClock_filtered/beast2/MC07/")
setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/molClock_filtered/")
#----
#---- LTT for annotated subtrees BEAST2 output with HPD --------------------------------------------

files = list(
	tree = "all_filtered_align-linsi_trim05_beast2_mc07_mcmc_annot.tre",
	factor = 1,
	HPD=TRUE)

rm(list=ls()[!ls() %in% c("geo", "files")])

treei = read.beast(files$tree)
tree = as.phylo(treei)

# Transform files
treed = as_tibble(treei)
names(treed) = gsub("!", "", names(treed))

# Select all annotated nodes _______________________________________________________________________
(annotations = sort(as.vector(treed$name[!is.na(treed$name)])))
tmp = c("Radiolaria", "Acantharia", "Acantharia-sensuLato", "Spumellaria", "Nassellaria", "Collodaria", "RAD-A", "RAD-B", "RAD-C", "RAD-X")

all(tmp %in% annotations)

annotations = tmp; rm(tmp)

# Loop through the main tree and the different annotated nodes to get the LTT data _________________
lttData = data.frame()
for(annot in annotations){
	cat("\r  Working on ", annot, " (", grep(annot, annotations), "/", length(annotations), ")                    ", sep="")
	node = as.numeric(subset(treed, name==annot)$node)
	tmp = tree_subset(treei, node=node, levels_back=0)
	tmp = as_tibble(tmp)
	tmp = subset(tmp, !is.na(height_median))
	tmp = subset(tmp, is.na(label))
	tmp = separate(tmp, 'height_0.95_HPD', c("HPD05", "HPD95"), sep=",")
	tmp$HPD05 = gsub("c\\(", "", tmp$HPD05)
	tmp$HPD95 = gsub("\\)", "", tmp$HPD95)
	tmp = tmp[order(as.numeric(tmp$height_median), decreasing=TRUE),]
	tmp = data.frame(tree=annot,
					 time=-as.numeric(tmp$height_median),
					 lineages=1:length(tmp$height_median),
					 lnLineages=log(1:length(tmp$height_median), exp(1)),
					 hpd05=-as.numeric(tmp$HPD05),
					 hpd95=-as.numeric(tmp$HPD95))
		
	lttData = rbind(lttData, tmp)
}; rm(annot, tmp, node); cat("\r  Done                                   ")

# Transform variables ______________________________________________________________________________
lttData$time = lttData$time * files$factor
if(files$HPD){lttData$hpd05 = lttData$hpd05 * files$factor}
if(files$HPD){lttData$hpd95 = lttData$hpd95 * files$factor}
lttData = subset(lttData, !is.na(time))
# lttData$hpd_min = ifelse(lttData$tree=="Main_tree", as.numeric(lttData$hpd %>% gsub("c\\(", "", .) %>% gsub(",.*", "", .)) * files$factor, NA)
# lttData$hpd_max = ifelse(lttData$tree=="Main_tree", as.numeric(lttData$hpd %>% gsub("\\)", "", .) %>% gsub(".*, ", "", .)) * files$factor, NA)

# Get summary per annotated groups _________________________________________________________________
lineages = data.frame()
for(annot in unique(lttData$tree)){
	ss = subset(lttData, tree==annot)
	lineages = rbind(lineages, data.frame(group=annot, 
										   lineages=max(ss$lineages), 
										   origin=min(ss$time),
										   slope=lm(time ~ lineages, ss)$coefficients[2],
										   slopeLog=lm(time ~ lnLineages, ss)$coefficients[2]))
}; rm(annot, ss)

# Beautifying the dataset for the plot _____________________________________________________________
sort(unique(lttData$tree))
lttData$tree = factor(lttData$tree, 
					   levels=c("Radiolaria", 
					   		 "Acantharia", "Acantharia-sensuLato", 
					   		 "Spumellaria", "Nassellaria", "Collodaria", 
					   		 "RAD-A", "RAD-B", "RAD-C", "RAD-X"))


{lttData$colour = as.character(lttData$tree)
	lttData$colour[which(lttData$colour=="Radiolaria")]="black"
	lttData$colour[which(lttData$colour=="Acantharia")]="yellow2"
	lttData$colour[which(lttData$colour=="Acantharia-sensuLato")]="yellow4"
	lttData$colour[which(lttData$colour=="Spumellaria")]="steelblue2"
	# lttData$colour[which(lttData$colour=="Entactinaria")]="steelblue4"
	lttData$colour[which(lttData$colour=="Nassellaria")]="springgreen3"
	lttData$colour[which(lttData$colour=="Collodaria")]="orangered3"
	lttData$colour[which(lttData$colour=="RAD-A")]="grey40"
	lttData$colour[which(lttData$colour=="RAD-B")]="purple"
	lttData$colour[which(lttData$colour=="RAD-C")]="grey60"
	lttData$colour[which(lttData$colour=="RAD-X")]="grey80"}

lttData$colour = factor(lttData$colour, 
						 levels=c("black", "yellow2", "yellow4", 
						 		 "steelblue2", #"steelblue4", 
						 		 "springgreen3", "orangered3", "grey40", "purple", "grey60", "grey80"))

# Selecting the dataset to be plotted ______________________________________________________________
data = lttData
# data = subset(lttData, tree=="" | tree=="")
data$hpd05 = ifelse(data$tree == "Radiolaria", data$hpd05, NA)
data$hpd95 = ifelse(data$tree == "Radiolaria", data$hpd95, NA)

if(files$HPD){
	geos = subset(geo, time >= min(data$hpd95, na.rm=TRUE))
}else{
	geos = subset(geo, time >= min(data$time))
}

(lttplot = ggplot(data, aes(x=time, y=lnLineages))+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_segment(aes(x=hpd05, xend=hpd95, y=lnLineages, yend=lnLineages), color="grey90", linewidth=2, lineend="round")+
		geom_line(aes(color=tree)) +
		# geom_point(aes(color=tree)) +
		geom_text(data=geos, aes(x=mid, y=max(data$lnLineages)+1, label=era), size=2.5)+
		geom_text(data=geos, aes(x=mid, y=max(data$lnLineages)+0.5, label=period), size=2)+
		scale_x_continuous(breaks=seq((round(min(c(data$hpd_max, data$time), na.rm=TRUE), -2)), 0, 100), 
						   minor_breaks=seq((round(min(c(data$hpd_max, data$time), na.rm=TRUE), -2)), 0, 50)) +
		scale_y_continuous(breaks=1:(round(max(data$lnLineages), 0))) +
		scale_color_manual(values=as.character(sort(unique(data$colour))))+
		theme_classic()+
		labs(y="Ln (N)", x="Time (Ma)"))

pdf(gsub("\\.tre$", "_LTT-subtress.pdf", files$tree), width=11.69, height=8.27, paper='special')
plot(lttplot)
dev.off()


# Exporting the dataset ____________________________________________________________________________
write.table(lttData, gsub("\\.tre$", "_LTT-data.tsv", files$tree), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#

#---- LTT from LTT table ---------------------------------------------------------------------------

files = list(
	table = "all_filtered_align-linsi_trim05_MC7_mcmc_annot_LTT.tsv",
	factor = -100,
	HPD=TRUE)

rm(list=ls()[!ls() %in% c("geo", "files")])

data =fread(files$table)
data = subset(data, tree!="main")
data$time = data$time * files$factor
data$hpd05 = ifelse(data$tree=="Radiolaria", data$hpd05 * files$factor, NA)
data$hpd95 = ifelse(data$tree=="Radiolaria", data$hpd95 * files$factor, NA)

sort(unique(data$tree))
data$tree = factor(data$tree, 
                       levels=c("Radiolaria", 
                                "Acantharia", "Acantharia-sensuLato", 
                                "Spumellaria", "Nassellaria", "Collodaria", 
                                "RAD-A", "RAD-B", "RAD-C", "RAD-X"))

data$colour = as.character(data$tree)
{data$colour[which(data$colour=="Radiolaria")]="black"
    data$colour[which(data$colour=="Acantharia")]="yellow3"
    data$colour[which(data$colour=="Acantharia-sensuLato")]="yellow4"
    data$colour[which(data$colour=="Spumellaria")]="steelblue3"
    data$colour[which(data$colour=="Nassellaria")]="springgreen3"
    data$colour[which(data$colour=="Collodaria")]="orangered3"
    data$colour[which(data$colour=="RAD-A")]="grey40"
    data$colour[which(data$colour=="RAD-B")]="purple"
    data$colour[which(data$colour=="RAD-C")]="grey60"
    data$colour[which(data$colour=="RAD-X")]="grey80"}

data$colour = factor(data$colour, 
                         levels=c("black", "yellow3", "yellow4", "steelblue3", "springgreen3", "orangered3", "grey40", "purple", "grey60", "grey80"))

# Selecting the dataset to be plotted ______________________________________________________________

if(files$HPD){
    geos = subset(geo, time >= min(data$hpd95, na.rm=TRUE))
}else{
    geos = subset(geo, time >= min(data$time))
}
geos$mid = apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(lttplot = ggplot(data, aes(x=time, y=lnLineages))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        geom_segment(aes(x=hpd05, xend=hpd95, y=lnLineages, yend=lnLineages), color="grey90", size=2, lineend="round")+
        geom_line(aes(color=tree)) +
        geom_text(data=geos, aes(x=mid, y=max(data$lnLineages)+1, label=era), size=2.5)+
        geom_text(data=geos, aes(x=mid, y=max(data$lnLineages)+0.5, label=period), size=2)+
        scale_x_continuous(breaks=seq((round(min(c(data$hpd95, data$time), na.rm=TRUE), -2)), 0, 100), 
                           minor_breaks=seq((round(min(c(data$hpd95, data$time), na.rm=TRUE), -2)), 0, 50)) +
        scale_y_continuous(breaks=1:(round(max(data$lnLineages), 0))) +
        scale_color_manual(values=as.character(sort(unique(data$colour))))+
        theme_classic()+
        labs(y="Ln (N)", x="Time (Ma)"))

pdf(gsub("\\.tsv$", "_plot.pdf", files$table), width=11.69, height=8.27, paper='special')
plot(lttplot)
dev.off()


#----  
#---- LTT from two different tables ----------------------------------------------------------------

rm(list=ls()[!ls() %in% c("geo")])

files = list(table1="mcmcTree/MC07/step2_25burnin/all_filtered_align-linsi_trim05_MC7_mcmc_annot_LTT.tsv",
			 factor1 = -100,
			 table2="beast2/MC07/all_filtered_align-linsi_trim05_beast2_mc07_mcmc_annot_LTT-data.tsv",
			 factor2 = 1,
			 out = "all_filtered_align-linsi_trim05_combined-MCMCtree7-beast27.pdf")

{tmp1 = fread(files$table1)
	tmp1$time = tmp1$time * files$factor1
	tmp1$hpd05 = tmp1$hpd05 * files$factor1
	tmp1$hpd95 = tmp1$hpd95 * files$factor1}

{tmp2 = fread(files$table2)
	tmp2$time = tmp2$time * files$factor2
	tmp2$hpd05 = tmp2$hpd05 * files$factor2
	tmp2$hpd95 = tmp2$hpd95 * files$factor2
	tmp2$colour = NULL}

all(names(tmp1) == names(tmp2))

{tmp1$file="MCMCTree"
	tmp2$file="BEAST2"
	tmp1$tree = paste0(tmp1$tree, "1")
	tmp2$tree = paste0(tmp2$tree, "2")}

data = rbind(tmp1, tmp2); rm(tmp1, tmp2)
data$tree = gsub("Acantharea", "Acantharia", data$tree)

# Beutifying the table _____________________________________________________________________________
data$file = factor(data$file, levels=c("MCMCTree", "BEAST2"))

unique(data$tree)
data = subset(data, tree=="Radiolaria1" | tree=="Acantharia1" | tree=="Acantharia-sensuLato1" | 
			  	tree=="Spumellaria1" | tree=="Nassellaria1" | tree=="Collodaria1"| 
			  	tree=="RAD-A1" | tree=="RAD-B1" | tree=="RAD-C1" | 
			  	tree=="Radiolaria2" | tree=="Acantharia2" | tree=="Acantharia-sensuLato2" | 
			  	tree=="Spumellaria2" | tree=="Nassellaria2" | tree=="Collodaria2"| 
			  	tree=="RAD-A2" | tree=="RAD-B2" | tree=="RAD-C2")

data$tree = factor(data$tree, 
				   levels=c("Radiolaria1", "Radiolaria2", 
				   		 "Acantharia1", "Acantharia2", "Acantharia-sensuLato1", "Acantharia-sensuLato2", 
				   		 "Spumellaria1", "Spumellaria2", "Nassellaria1", "Nassellaria2", "Collodaria1", "Collodaria2", 
				   		 "RAD-A1", "RAD-A2", "RAD-B1", "RAD-B2", "RAD-C1", "RAD-C2"))

write.table(data, sub("\\.pdf", ".tsv", files$out), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

{data$colour = as.character(data$tree)
	data$colour[which(data$colour=="Radiolaria1")]="black"
	data$colour[which(data$colour=="Radiolaria2")]="grey20"
	data$colour[which(data$colour=="Acantharia1")]="yellow2"
	data$colour[which(data$colour=="Acantharia2")]="yellow1"
	data$colour[which(data$colour=="Acantharia-sensuLato1")]="yellow4"
	data$colour[which(data$colour=="Acantharia-sensuLato2")]="yellow3"
	data$colour[which(data$colour=="Spumellaria1")]="steelblue3"
	data$colour[which(data$colour=="Spumellaria2")]="steelblue2"
	data$colour[which(data$colour=="Nassellaria1")]="springgreen3"
	data$colour[which(data$colour=="Nassellaria2")]="springgreen2"
	data$colour[which(data$colour=="Collodaria1")]="orangered3"
	data$colour[which(data$colour=="Collodaria2")]="orangered2"
	data$colour[which(data$colour=="RAD-A1")]="grey30"
	data$colour[which(data$colour=="RAD-A2")]="grey40"
	data$colour[which(data$colour=="RAD-B1")]="purple2"
	data$colour[which(data$colour=="RAD-B2")]="purple1"
	data$colour[which(data$colour=="RAD-C1")]="grey50"
	data$colour[which(data$colour=="RAD-C2")]="grey60"}

data$colour = factor(data$colour, 
					 levels=c("black", "grey20", 
					 		 "yellow2", "yellow1", "yellow4", "yellow3", 
					 		 "steelblue3", "steelblue2", "springgreen3", "springgreen2", "orangered3", "orangered2", 
					 		 "grey30", "grey40", "purple2", "purple1", "grey50", "grey60"))

data$hpd05 = ifelse(grepl("Radiolaria", data$tree), data$hpd05, NA)
data$hpd95 = ifelse(grepl("Radiolaria", data$tree), data$hpd95, NA)

# Plotting _________________________________________________________________________________________

geos = subset(geo, time >= min(data$hpd95, na.rm=TRUE))

(lttplot = ggplot(data)+
		geom_vline(xintercept = geos$time, color="lightgrey") +
		geom_segment(aes(x=hpd05, xend=hpd95, y=lnLineages, yend=lnLineages, color=tree), alpha=0.2, linewidth=2, lineend="round")+
		geom_line(aes(x=time, y=lnLineages, color=tree)) +
		geom_text(data=geos, aes(x=mid, y=max(data$lnLineages)+1, label=era), size=2.5)+
		geom_text(data=geos, aes(x=mid, y=max(data$lnLineages)+0.5, label=period), size=2)+
		scale_x_continuous(breaks=seq((round(min(c(data$hpd95, data$time), na.rm=TRUE), -2)), 0, 100), 
						   minor_breaks=seq((round(min(c(data$hpd95, data$time), na.rm=TRUE), -2)), 0, 50)) +
		scale_y_continuous(breaks=1:(round(max(data$lnLineages), 0))) +
		scale_color_manual(values=as.character(sort(unique(data$colour))))+
		theme_classic()+
		labs(y="Ln (N)", x="Time (Ma)"))

pdf(files$out, width=11.69, height=8.27, paper='special')
plot(lttplot)
dev.off()

#

#----  
#---- Diversification analysis with BAMM tools -----------------------------------------------------

setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/molClock_filtered_subset/mcmcTree/MC06/diver/div90/")

rm(list=ls()[!ls() %in% c("geo", "files")])

files = list(tree="all_filtered_align-linsi_trim05_subset_MC6_combined_rads.tre",
             treeOut=gsub("\\.[^\\.]+$", "_Ma.tre", files$tree),
             div=0.9,
             edata="test_auto_event_data.txt", 
             factor=100)

tree = read.tree(files$tree)
tree$edge.length = tree$edge.length * files$factor
write.tree(tree, files$treeOut)

plot(tree, cex=.6)
nodelabels(frame = "none")
axisPhylo()

# for(diversity_fraction in c(0.9, 0.5, 0.1)){
#     setBAMMpriors(tree, 
#                   total.taxa=Ntip(tree)+round(Ntip(tree)*(1-diversity_fraction),0),
#                   outfile=paste0("priors_BAMM_diversity", diversity_fraction*100, ".txt"))
# }

# Generate customized control file
priors = setBAMMpriors(tree, total.taxa=Ntip(tree)+round(Ntip(tree)*(1-files$div),0),
                        outfile = NULL)
generateControlFile(file='divcontrol.txt', type = 'diversification',
                    params = list(
                        treefile = files$tree,
                        outName = "test_auto",
                        globalSamplingFraction = files$div,
                        numberOfGenerations = '1000000',
                        overwrite = '1',
                        lambdaInitPrior = as.numeric(priors['lambdaInitPrior']),
                        lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']),
                        muInitPrior = as.numeric(priors['muInitPrior']),
                        expectedNumberOfShifts = 3))

system("BAMM -c divcontrol.txt")

edata = getEventData(tree, eventdata = files$edata, burnin=0.25)

head(edata$eventData)
bamm.tree = plot.bammdata(edata, lwd=2, labels = T, cex = 0.5)
addBAMMshifts(edata, cex=2)
addBAMMlegend(bamm.tree)

plotRateThroughTime(edata, ratetype="speciation")
plotRateThroughTime(edata, ratetype="netdiv")
plotRateThroughTime(edata, ratetype="extinction")




allrates = getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

best = getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
plot(best, lwd=2)
addBAMMshifts(best, cex=2)




mcmcout = read.csv("test_auto_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart = floor(0.1 * nrow(mcmcout))
postburn = mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

post_probs = table(postburn$N_shifts) / nrow(postburn)

names(post_probs)

post_probs['X'] / post_probs['Y']

tree = read.tree("all_filtered_align-linsi_trim05_subset_MC6_combined.newick")
edata = getEventData(tree, eventdata = "test_d50_event_data.txt", burnin=0.1)
shift_probs = summary(edata)

plot.bammdata(edata, lwd=2)
plot.bammdata(edata, lwd=2, legend=T)

css = credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)

# plot.credibleshiftset(css)

# best = getBestShiftConfiguration(edata, expectedNumberOfShifts=1)
# plot.bammdata(best, lwd = 2)
# addBAMMshifts(best, cex=2.5)

allrates = getCladeRates(edata)
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))

plotRateThroughTime(edata, ratetype="speciation")
plotRateThroughTime(edata, ratetype="netdiv")
# plotRateThroughTime(edata, ratetype="extinction")


#----  
#---- Read Rdata from ClaDS function (julia) -------------------------------------------------------

files = list(dir="div90/",
             file="")

load(paste0(files$dir, files$file))

rm(list=ls()[!ls() %in% c("geo", "files", "CladsOutput")])

# Checking output from ClaDS Julia _________________________________________________________________

# Plot rates of diversification in the tree
tree = CladsOutput$tree
rates = log(CladsOutput$lambdai_map)

# plot(tree)
# axisPhylo()

pdf(paste0(files$dir, gsub("\\.[^\\.]+$", "_RTTtree.pdf", files$file)), width=11.69, height=8.27, paper='special')
plot_ClaDS_phylo(tree, rates, show.tip.label=TRUE, cex=0.2, log=FALSE)
dev.off()

ltt = as.data.frame(ltt.plot.coords(CladsOutput$tree))
ltt$logN = log(ltt$N, exp(1))

# Plotting the rate through time plot ______________________________________________________________

tmp = apply(data.frame(CladsOutput$time_points, c(CladsOutput$time_points[-1], 0)), 1, mean)
clads_rtt = data.frame(time=tmp[-length(tmp)], 
                        rate=CladsOutput$RTT_map); rm(tmp)

burn = 100
chains = data.frame(toRemove=rep(NA, length(CladsOutput$rtt_chains[[1]][[1]])))
c = 0
for(chain in CladsOutput$rtt_chains){
    c = c + 1
    chaini = as.data.frame(chain)
    colnames(chaini) = paste0("Chain", c, "_iter", 1:ncol(chaini))
    chaini = chaini[,-c(1:burn)]
    chains = cbind(chains, chaini)
}; rm(c, chaini)
chains$toRemove = NULL

library(HDInterval)
hpd = as.data.frame(t(apply(t(chains), 2, hdi)))
clads_rtt$HPD_95_low = hpd$lower
clads_rtt$HPD_95_high = hpd$upper


geos = subset(geo, time>-max(clads_rtt$time)*100)
geos$mid = apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(rttPlot =ggplot(clads_rtt, aes(x=-rev(time)*100, y=rate))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        # geom_point()+
        geom_line(size=2)+
        geom_line(aes(y=HPD_95_low), colour="lightblue")+
        geom_line(aes(y=HPD_95_high), colour="lightblue")+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
        scale_x_continuous(breaks=seq((round(max(clads_rtt$time)*-100, -2)), 0, 100),
                           minor_breaks=seq((round(min(clads_rtt$time)*-100, -2)), 0, 50)) +
        theme_classic()+
        labs(y="Diversification rate (Ln)", x="Time (Ma)"))

pdf(paste0(files$dir, gsub("\\.[^\\.]+$", "_RTT.pdf", files$file)), width=11.69, height=8.27, paper='special')
plot(rttPlot)
dev.off()

# Plotting the Lineages through time plot __________________________________________________________

tmp = apply(data.frame(CladsOutput$time_points, c(CladsOutput$time_points[-1], 0)), 1, mean)
clads_ltt = data.frame(time=tmp, 
                        lineages=CladsOutput$DTT_mean); rm(tmp)

geos = subset(geo, time>-max(clads_ltt$time)*100)
geos$mid = apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

(lttPlot =ggplot(clads_ltt, aes(x=-rev(time)*100, y=log(lineages)))+
        geom_vline(xintercept = geos$time, color="lightgrey") +
        # geom_point()+
        geom_line(size=2)+
        # scale_y_log10() + annotation_logticks(sides = 'l')+
        scale_x_continuous(breaks=seq((round(max(clads_ltt$time)*-100, -2)), 0, 100),
                           minor_breaks=seq((round(min(clads_ltt$time)*-100, -2)), 0, 50)) +
        theme_classic()+
        labs(y="Diversity Through Time (Ln)", x="Time (Ma)"))

pdf(paste0(files$dir, gsub("\\.[^\\.]+$", "_DTT.pdf", files$file)), width=11.69, height=8.27, paper='special')
plot(lttPlot)
dev.off()





#

# Checking output from two different ClaDS Julia output ____________________________________________

rm(list=ls()[!ls() %in% c("geo")])

groups = list(Radiolaria="clades/clade_10_Radiolaria.Rdata",
               Foraminifera="clades/clade_7_Foraminifera_clads.Rdata",
               Diatomea="clades/clade_126_Diatomea.Rdata",
               Dinophycea="clades/clade_174_Dinophyceae.Rdata")

rtt = data.frame()
getRTTfromCladsOutput = function(cladsOutput, verbose=TRUE){
    if(verbose){cat("  Reading file\n")}
    tmp = apply(data.frame(cladsOutput$time_points, c(cladsOutput$time_points[-1], 0)), 1, mean)
    clads_rtt = data.frame(time=tmp[-length(tmp)], 
                            rate=cladsOutput$RTT_map); rm(tmp)
    
    if(verbose){cat("  Extracting chains\n")}
    burn = 100
    chains = data.frame(toRemove=rep(NA, length(cladsOutput$rtt_chains[[1]][[1]])))
    c = 0
    for(chain in cladsOutput$rtt_chains){
        c = c + 1
        cat("\r    ", c, "/", length(cladsOutput$rtt_chains), sep="")
        chaini = as.data.frame(chain)
        colnames(chaini) = paste0("Chain", c, "_iter", 1:ncol(chaini))
        chaini = chaini[,-c(1:burn)]
        chains = cbind(chains, chaini)
    }; rm(c, chaini)
    chains$toRemove = NULL
    
    if(verbose){cat("\n  Calculating HPD intervals\n")}
    hpd = as.data.frame(t(apply(t(chains), 2, hdi)))
    clads_rtt$HPD_95_low = hpd$lower
    clads_rtt$HPD_95_high = hpd$upper
    
    clads_rtt$revTime = -rev(clads_rtt$time)
    if(verbose){cat("Done\n")}
    return(clads_rtt)
}
for(g in names(groups)){
    cat("Loading data for '", g ,"' (", grep(g, names(groups)), "/", length(groups), ")\n", sep="")
    load(groups[[g]])
    tmp = CladsOutput
    
    tmpt = getRTTfromCladsOutput(cladsOutput=tmp)
    tmpt$group = g
    
    rtt = rbind(rtt, tmpt)
}; rm(g, CladsOutput, tmp, tmpt)


geos = subset(geo, time>-max(rtt$time))
geos$mid = apply(data.frame(geos$time, c(geos$time[-1], 0)), 1, mean)

rttplot = ggplot(subset(rtt, revTime>(-320)), aes(x=revTime, y=rate, colour=group))+
    # geom_vline(xintercept = geos$time, color="lightgrey") +
    geom_vline(xintercept = subset(geos, time>(-320))$time, color="lightgrey") +
    # geom_point()+
    geom_line(size=2)+
    geom_line(aes(y=HPD_95_low, colour=group))+
    geom_line(aes(y=HPD_95_high, colour=group))+
    scale_y_log10() + annotation_logticks(sides = 'l')+
    theme_classic()+
    labs(y="Diversification rate (Ln)", x="Time (Ma)", 
         title=paste0(unique(rtt$group), collapse=" - "))
rttplot

pdf("rttPlots/####.pdf", width=11.69, height=8.27, paper='special')
plot(rttplot)
dev.off()

#
#----
#