
library(data.table)
library(ggplot2)

geo = data.frame(time= c(-2500, -2300, -2050, -1800, -1600 , -1400, -1200, -1000, -720, -635, -541, -485.4, -443.8, -419, -358.9, -298.9, -251.9, -201.4, -145, -66, -23.03, -2.58),
				 period=c("Siderian", "Rhyacian", "Orosirian", "Statherian", "Calymnian","Ectasian","Stenian","Tonian","Cryogenian","Ediacran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene","Quaternary"),
				 era=c("Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Paleoproterozoic","Mesoproterozoic","Mesoproterozoic","Mesoproterozoic","Neoproterozoic","Neoproterozoic","Neoproterozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic","Paleozoic",
				 	  "Mesozoic","Mesozoic","Mesozoic","Cenozoic","Cenozoic","Cenozoic"))
geo$mid = apply(data.frame(geo$time, c(geo$time[-1], 0)), 1, mean)

setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/bbb/")

data = fread("ages.tsv")

data$clade = factor(data$clade, levels=rev(c("Spumellaria", "Cladococcoidea", "Rhizosphaeroidea", "Centrocubidae", "Excentroconchidae", "Clade-M", "Trematodiscoidea", "Haliommoidea", "Clade-L", "Hexacromyidae", "Spongosphaeroidea", "Spongodiscoidea", "Lithocyclioidea", "Nassellaria", "Eucyrtidioidea", "Carpocanioidea", "Cycladophoroidea", "Lithochytridoidea", "Lithochytridoidea+Nass4", "Pterocorythoidea", "Archipillioidea", "Theopilioidea", "Plagiacanthoidea", "Acanthodesmioidea", "Artostrobioidea", "Artostrobioidea+Env", "Plectopyramidoidea", "Oroscenoidea", "Collodaria", "Sphaerozoidae", "Collosphaeridae")))
data$method = factor(data$method, levels=c("BBB", "MCMCTree", "BEAST2"))
data$median = -as.numeric(data$median)
data$HPD05 = -as.numeric(data$HPD05)
data$HPD95 = -as.numeric(data$HPD95)

geos = subset(geo, time>min(data$HPD95))

(plot = ggplot(data)+
	geom_vline(xintercept = geos$time, color="lightgrey") +
	geom_segment(aes(x=HPD05, xend=HPD95, y=clade, yend=clade, colour=method), linewidth=1, lineend="round")+
	geom_point(aes(x=median, y=clade, colour=method))+
	scale_x_continuous(breaks=seq(min(round(data$HPD95, -2)), 0, 100),
					   minor_breaks=seq(min(round(data$HPD95, -2)), 0, 50))+
	theme_classic())

pdf("ages.pdf", width=11.69, height=8.27, paper='special'); plot(plot); dev.off()
