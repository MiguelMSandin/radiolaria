
library(data.table)
library(dplyr)
library(ggplot2)

setwd("/home/miguel/Desktop/PhD/0_Thesis/2_chapter/data/molClock_filtered/beast2/")

data = fread("calibrations_minMax.tsv")

head(data)
data = as.data.frame(data %>% group_by(Node) %>% mutate(mean=mean(c(min, max)), sd=sd(c(min, max))))
data$newMean = data$mean - data$min
data$logNewMean = abs(log10(data$newMean))
data$logsd = abs(log10(data$sd))


pdf("calibrations_mc07_logNormal.pdf", width=11.69, height=8.27, paper='special')
for(i in 1:nrow(data)){
	cat("\r", round(i/nrow(data)*100), "%", sep="", end="")
	tmp = data.frame(x=rlnorm(1000, data$logNewMean[i], data$logsd[i]) + data$min[i])
	probs = qlnorm(c(0.05, 0.5, 0.95), data$logNewMean[i], data$logsd[i])+data$min[i]
	plot = ggplot(tmp, aes(x=x)) +	
		geom_density(alpha=0.2) + 
		geom_vline(xintercept=probs, colour="red4")+
		labs(title=paste0(data$Node[i], " (", data$min[i], "-", data$max[i], ")"),
			 subtitle=paste0("logNormal(µ=", round(data$logNewMean[i], 2), ", sd=", round(data$logsd[i], 2), "): ",
			 				"p05=", round(probs[1],2), " | P50=", round(probs[2],2)," | p95=", round(probs[3],2)),
			 x=NULL)+
		theme_minimal()
	plot(plot)
}; rm(i, tmp, probs)
dev.off()


# First, find parameters for a specific group ______________________________________________________
rm(list=ls()[!ls() %in% c("data")])
(i = grep("Tetra", data$Node))
# i=12
(node = data$Node[i])
(param = list(min=data$min[i], max=data$max[i],
			 mean = 0, sd = 0, offset = 0))
# param$min = round(param$min, 1)
# param$max = round(param$max, 1)
limits = list(meanLow=0.5, meanUp=1, sdLow=1, sdUp=2, steps=0.01, offsetPct=0.02, threshold=0.001)
for(o in seq(0, param$min*limits$offsetPct, by=0.1)){
	cat("\r  ", node, " (", grep(node, data$Node), "/", length(data$Node), ") ", 
		round((o/(param$min*limits$offsetPct))*100, 2), "%                    ", sep="", end="")
	for(m in seq(limits$meanLow, limits$meanUp, by = limits$steps)){
		for(s in seq(limits$sdLow, limits$sdUp, by = limits$steps)){
			low = plnorm(o, m, s)
			up = plnorm(param$max-param$min+o, m, s)
			if(0.05-limits$threshold < low & low < 0.05+limits$threshold & 0.95-limits$threshold < up & up < 0.95+limits$threshold){
				param$mean = c(param$mean, m)
				param$sd = c(param$sd, s)
				param$offset = c(param$offset, param$min-o)
			}
		}
	}
}; cat("\rDone                                        \n"); rm(o, m, s, low, up)
param
param$mean = param$mean[-1]; param$sd = param$sd[-1]; param$offset = param$offset[-1]
hist(param$mean, main=paste0("Mean (n=", length(param$mean), ")"))
hist(param$sd, main=paste0("SD (n=", length(param$sd), ")"))
hist(param$offset, main=paste0("Offset (n=", length(param$offset), ")"))
length(param$mean) == length(param$sd) & length(param$sd) == length(param$offset)

# Now choose the best scoring parameters ___________________________________________________________
chooseParam = data.frame()
# chooseSlope = data.frame()
for(i in 1:length(param$mean)){
	cat("\r  ", round(i/length(param$mean)*100), "%", sep="")
	chooseParam = rbind(chooseParam, 
						data.frame(mean=param$mean[i],
								   sd=param$sd[i],
								   offset=param$offset[i],
								   Pmin=plnorm(param$min-param$offset[i], param$mean[i], param$sd[i]),
								   Pmax=plnorm(param$max-param$offset[i], param$mean[i], param$sd[i])))
	# tmp = rlnorm(100, param$mean[i], param$sd[i])
	# chooseSlope = rbind(chooseSlope, data.frame(Function=paste0("logNormal(", param$mean[i], ",", param$sd[i], ",", param$offset[i],")"),
												# values=tmp, dates=tmp+param$offset[i]))
}; rm(i)
{chooseParam$diff05 = abs(0.05-chooseParam$Pmin)/0.05*100
chooseParam$diff95 = abs(0.95-chooseParam$Pmax)/0.95*100
chooseParam$diff05Rank = rank(chooseParam$diff05)
chooseParam$diff95Rank = rank(chooseParam$diff95)
chooseParam$diffCombined = chooseParam$diff05 + chooseParam$diff95
chooseParam$diffCombinedRank = rank(chooseParam$diffCombined)}

# ggplot(chooseSlope, aes(x=dates, colour=Function))+	geom_density(alpha=0.2)
subset(chooseParam, diffCombinedRank <= 3 | diff05Rank <= 3 | diff95Rank <= 3)
cat("  ", node, " (", param$min, "-", param$max, ")",
	"\tMean: ", subset(chooseParam, diffCombinedRank == 1)$mean, "\tSD: ", subset(chooseParam, diffCombinedRank == 1)$sd, 
	"\tOffset:", subset(chooseParam, diffCombinedRank == 1)$offset,
	sep="", end="\n")

#
