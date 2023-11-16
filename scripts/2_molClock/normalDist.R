#----
#---- Load packages --------------------------------------------------------------------------------

library(ggplot2)
library(HDInterval)
library(data.table)

#----
#---- Get single estimates -------------------------------------------------------------------------

lim = list(min=182.7, max=247.2)
m = mean(c(lim$min, lim$max))

# from the probabilities
{error= 0.001     # percentage of similarity to the limits
topTries = 0.2   # percentage of the mean
HPD = 0.95}

{sd = data.frame(sd=c(), distMin=c(), distMax=c())
for(i in seq(0, m*topTries, 0.01)){
	cat("\r", round(i/(m*topTries)*100), "%", sep="")
	pmin = qnorm(0.5-HPD/2, m, i)
	pmax = qnorm(0.5+HPD/2, m, i)
	if(((lim$min - lim$min*error) < pmin) & (pmin < (lim$min + lim$min*error)) & ((lim$max - lim$max*error) < pmax) & (pmax < (lim$max + lim$max*error))){
		sd = rbind(sd, data.frame(sd=i, 
								  # diffMin=pmin-lim$min,
								  distMin=abs(pmin-lim$min),  
								  # diffMax=pmax-lim$max, 
								  distMax=abs(pmax-lim$max)))
	}
}; rm(i, pmin, pmax); cat("\n")
	sdmin = sd$sd[which(sd$distMin == min(sd$distMin))]
	sdmax = sd$sd[which(sd$distMax == min(sd$distMax))]
	if(sdmin == sdmax){
		cat("  Mean set at", m,
			"\n  Standard deviation chosen at", sdmin, "\n")
	}else{cat("Something went wrong... please check")}}

# Now plot it
data = data.frame(x=rnorm(10^4, m, sdmin))
ggplot(data, aes(x=-x)) + 
	geom_vline(xintercept=-qnorm(0.5-HPD/2, m, sdmin), colour="skyblue2")+
	geom_vline(xintercept=-qnorm(0.5+HPD/2, m, sdmin), colour="skyblue2")+
	geom_vline(xintercept=-m, colour="skyblue2")+
	geom_segment(x=-m+sdmin, xend=-m-sdmin, y=0.1, yend=0.1, colour="skyblue2")+
	labs(title=paste0("N(", m, ", ", sdmin, ")"),
		 subtitle=paste0("HPD: ", round(qnorm(0.5+HPD/2, m, sdmin), 1), "-", round(qnorm(0.5-HPD/2, m, sdmin), 1)))+
	scale_x_continuous(limits = c(-max(data$x), 0)) +
	geom_density()+
	theme_minimal()

# Find a SD that still falls within the range but with a given mean ________________________________
lim = list(min=0, max=1)
m = 0.5

{sd = data.frame(sd=c(), distMin=c(), distMax=c())
	for(i in seq(0, m*topTries, 0.01)){
		cat("\r", round(i/(m*topTries)*100), "%", sep="")
		pmin = qnorm(0.5-HPD/2, m, i)
		pmax = qnorm(0.5+HPD/2, m, i)
		if(((lim$min - lim$min*error) < pmin) & (pmax < (lim$max + lim$max*error))){
			sd = rbind(sd, data.frame(sd=i, 
									  # diffMin=pmin-lim$min,
									  distMin=abs(pmin-lim$min),  
									  # diffMax=pmax-lim$max, 
									  distMax=abs(pmax-lim$max)))
		}
	}; rm(i, pmin, pmax); cat("\n")
	sdmin = sd$sd[which(sd$distMin == min(sd$distMin))]
	sdmax = sd$sd[which(sd$distMax == min(sd$distMax))]
	if(sdmin == sdmax){
		cat("  Mean set at", m,
			"\n  Standard deviation chosen at", sdmin, "\n")
	}else{
		cat("  Mean set at", m,
			"\n  Standard deviation could be either", sdmin, "or", sdmax, 
			"\n    And the mean between these two sd is:", mean(c(sdmin, sdmax)))
	}}

# Now plot it
SD = mean(c(sdmin, sdmax))
data = data.frame(x=rnorm(10^4, m, SD))
ggplot(data, aes(x=-x)) + 
	geom_vline(xintercept=-qnorm(0.5-HPD/2, m, SD), colour="skyblue2")+
	geom_vline(xintercept=-qnorm(0.5+HPD/2, m, SD), colour="skyblue2")+
	geom_vline(xintercept=-m, colour="skyblue2")+
	geom_segment(x=-m+SD, xend=-m-SD, y=0.1, yend=0.1, colour="skyblue2")+
	geom_density()+
	theme_minimal()

#----
#---- Plot normal distribution from table ----------------------------------------------------------

data = fread("Desktop/PhD/0_Thesis/2_chapter/data/molClock_2311/calibrations_mc09.tsv")

HPD = 0.95

pdf("Desktop/PhD/0_Thesis/2_chapter/data/molClock_2311/beast2/calibration_distributions.pdf", width=11.69, height=8.27, paper='special')
for(i in 1:nrow(data)){
	ss = data[i,]
	if(ss$distribution == "Normal"){
		cat("\r  Plotting", ss$node, " (", i, "/", nrow(data), ")", sep="")
		tmp = data.frame(x=rnorm(10^5, ss$mean, ss$sd))
		plot = ggplot(tmp, aes(x=-x)) + 
			geom_vline(xintercept=-ss$max, colour="#FFB000")+
			geom_vline(xintercept=-ss$min, colour="#FFB000")+
			geom_vline(xintercept=-qnorm(0.5-HPD/2, ss$mean, ss$sd), colour="#648FFF")+
			geom_vline(xintercept=-qnorm(0.5+HPD/2, ss$mean, ss$sd), colour="#648FFF")+
			geom_vline(xintercept=-ss$mean, colour="#648FFF")+
			geom_segment(x=-ss$mean+ss$sd, xend=-ss$mean-ss$sd, y=0.1, yend=0.1, colour="#648FFF")+
			labs(title=paste0(ss$node, " N(", ss$mean, ", ", ss$sd, ") Ma"),
				 subtitle=paste0("Prior HPD:        [", ss$max, "-", ss$min, "] Ma",
				 				"\nPosterior HPD: [", round(qnorm(0.5+HPD/2, ss$mean, ss$sd), 1), "-", round(qnorm(0.5-HPD/2, ss$mean, ss$sd), 1), "] Ma"))+
			# scale_x_continuous(limits = c(-max(data$x), 0)) +
			geom_density()+
			theme_minimal()
	}
	plot(plot)
}; rm(i, ss, tmp)
dev.off()

#----
