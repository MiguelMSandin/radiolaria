#----
#---- Libraries ------------------------------------------------------------------------------------

library(ape)
library(treeio)
library(data.table)
library(dplyr)
library(HDInterval)
library(ggplot2)

# geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
# 				 period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian","Ectasian","Stenian","Tonian","Cryogenian","Ediacran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene","Quaternary"),
# 				 era=c("Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Mesoproterozoic","Mesoproterozoic","Mesoproterozoic","Neoproterozoic","Neoproterozoic","Neoproterozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic",
# 				 	  "Mesozoic","Mesozoic","Mesozoic","Cenozoic","Cenozoic","Cenozoic"))
# geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

#----
setwd("~/rads/data/molClock/beast2")
#----
#---- Set file names -------------------------------------------------------------------------------

rm(list=ls()[!ls() %in% c("files")])

files = list(trees = "",
			 tree = "",
			 log = "",
			 nodeHeights="",
			 calibrations = "",
			 plot="")

#----
#---- Explore prior and posterior ------------------------------------------------------------------

log = fread(files$log)
post = log %>% select(any_of(grep("mrca.ge", names(log), value=TRUE)))
colnames(post) = names(post) %>% sub("mrca\\.ge\\(", "", .) %>% sub("\\)", "", .)
colnames(post) = names(post) %>% sub("Cladococcoidea", "Cladococcus", .)

post = melt(as.data.table(post))
colnames(post) = c("node", "value")

# Generate prior file
calib = fread(files$calibrations)

prior = data.frame()
for(i in unique(post$node)){
	if(any(grepl(i, calib$node))){
		n = nrow(log)
		if(calib$distribution[which(calib$node == i)] == "Uniform"){
			cat("  Applying a Uniform distribution to", i, "\n")
			min = calib$min[which(calib$node == i)]
			max = calib$max[which(calib$node == i)]
			tmp = runif(n, min, max)
		}else if(calib$distribution[which(calib$node == i)] == "Normal"){
			cat("  Applying a Normal distribution to", i, "\n")
			mean = calib$mean[which(calib$node == i)]
			sd = calib$sd[which(calib$node == i)]
			tmp = rnorm(n, mean, sd)
		}else if(calib$distribution[which(calib$node == i)] == "logNormal"){
			cat("  Applying a logNormal distribution to", i, "\n")
			mean = calib$mean[which(calib$node == i)]
			sd = calib$sd[which(calib$node == i)]
			min = calib$min[which(calib$node == i)]
			tmp = rlnorm(n, mean, sd) + min
		}else{cat("  No distribution found for", i, "\n")}
		prior = rbind(prior, data.frame(node=i, value=tmp))
	}else{cat(" ", i, "not found\n")}
}; rm(i, n, min, max, mean, sd, tmp)

# clean the POSTERIOR data from the 1% lowest posterior density
postc = data.frame()
for(i in unique(post$node)){
	ss = subset(post, node==i)
	hpd = hdi(ss$value, credMass = 0.99)
	cat("  Cleaning ", i, " - HPD: ", round(hpd[1]), "-", round(hpd[2]), "\n", sep="")
	out = subset(ss, value >= hpd[1] & value <= hpd[2])
	postc = rbind(postc, out)
}; rm(i, hpd, ss, out)

# clean the PRIOR data using the maximum of the observed data
priorc = data.frame()
for(i in unique(prior$node)){
	if(i == "Excentroconchidae"){
		ss = subset(prior, node==i)
		max = quantile(ss$value, 0.95)
		cat("  Cleaning", i, "from values above", max, "\n")
		out = subset(ss, value <= max & value > 0)
		priorc = rbind(priorc, out)
	}else{
		max = max(subset(postc, node == i)$value)
		cat("  Cleaning", i, "from values above", max, "\n")
		ss = subset(prior, node==i)
		out = subset(ss, value <= max & value > 0)
		priorc = rbind(priorc, out)
	}
}; rm(i, max, ss, out)

# Merge prior and posterior
postc$dist = "posterior"
priorc$dist = "prior"
data = rbind(postc, priorc)

tmp = calib %>% select(c("node", "min", "max"))
tmp$calib = paste0(tmp$node, " (", tmp$max, "-", tmp$min, ")")
tmp = tmp %>% group_by(calib) %>% mutate(mean=mean(min, max))
data = merge(data, tmp, by="node")
data$calib = factor(data$calib, levels=tmp$calib[order(tmp$mean, decreasing=TRUE)]); rm(tmp)

data = subset(data, node != "Root")

(plotPriPos = ggplot(data) +
		geom_density(aes(x=-value, fill=dist, colour=dist), alpha=0.4)+
		facet_wrap(~calib, scales="free")+
		scale_colour_manual(values=c("#648FFF", "#FFB000"))+
		scale_fill_manual(values=c("#648FFF", "#FFB000"))+
		theme_bw())

pdf(paste0(files$plot, "_prior-posterior.pdf"), width=11.69, height=8.27, paper='special')
plot(plotPriPos)
dev.off()

#----
#---- Get posterior distributions of specific nodes ------------------------------------------------

# After running the function 'treeGetNodeLengths.py' on all the trees from the mcmc run
paste0("treeGetNodeLengths.py -t ", sub("_nodeHeights\\.tsv", ".trees", files$nodeHeights))

nodes = list(rads=3, 
			 polycystines=4, Spumellaria=5, SpumSymb=12, 
			 Nassellaria=192, Collodaria=266, 
			 Spasmaria=346, RadsBC=347, RadB=348, RadC=447,
			 # AcanRadA=465, AcanthariaSensuLato=466, AcanthariaSensuStricto=598, AcanSymb=473, RadA=606)
			 AcanRadA=465, AcanthariaSensuLato=466, AcanthariaSensuStricto=468, AcanSymb=473, RadA=606)

# data = fread(files$nodeHeights)
# datas = select(data, paste0("n", nodes))

system(paste("cut -f", paste(nodes, collapse=","), files$nodeHeights, "> tmp.tsv"))
datas = fread("tmp.tsv")
system("rm -f tmp.tsv")

colnames(datas) = names(nodes)

apply(datas, 2, summary)

pdf(paste0(files$plot, "_distributionNodes_single.pdf"), width=11.69, height=8.27, paper='special')
for(n in names(datas)){
	ss = select(datas, all_of(n))
	colnames(ss) = "x"
	plotDist = ggplot(ss)+
		geom_density(aes(x=-x), fill="grey80")+
		labs(title=paste0(n))+
		theme_minimal()
	plot(plotDist)
}; rm(n, ss)
dev.off()

datam = melt(datas)

datamc = data.frame()
for(clade in unique(datam$variable)){
	ss = subset(datam, variable==clade)
	hpd = hdi(ss$value)
	ss = subset(ss, value >= hpd[1] & value <= hpd[2])
	datamc = rbind(datamc, ss)
}; rm(clade, ss, hpd)

(plotDist = ggplot(datamc)+
	geom_density(aes(x=-value), fill="grey80")+
	facet_wrap(~variable, scales="free")+
	theme_minimal())

pdf(paste0(files$plot, "_distributionNodes.pdf"), width=11.69, height=8.27, paper='special')
plot(plotDist)
dev.off()
#

# Now plot specific node per cycle
datas$cycle = 1:nrow(datas)
pdf(paste0(files$plot, "_ageCycle.pdf"), width=11.69, height=8.27, paper='special')
for(n in grep("cycle", names(datas), value=TRUE, invert=TRUE)){
	ss = select(datas, all_of(c(n, "cycle")))
	colnames(ss) = c("x", "cycle")
	plotDist = ggplot(ss)+
		geom_line(aes(x=cycle, y=x))+
		labs(title=paste0(n))+
		theme_minimal()
	plot(plotDist)
}; rm(n, ss)
dev.off()

#----
#
