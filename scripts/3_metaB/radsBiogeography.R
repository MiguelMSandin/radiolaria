#----
#---- Libraries ------------------------------------------------------------------------------------

library(data.table)
library(dplyr)

library(vegan)
library(pastecs)
library(moments)

library(ggplot2)
library(treemapify)

#----
#---- Preparing final tables -----------------------------------------------------------------------

# setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/metaB/")
# 
# # Working with the TAXONOMIC table _________________________________________________________________
# # Opening new taxonomic assignation
# taxo = fread("otus/radiolaria_otusTaxo_PR2v4.14rads_211216.tsv", sep="\t")
# colnames(taxo) = c("md5sum", "target", "id", "alnlen", "mism", "opens", "qlo", "qhi", "tlo", "thi", "evalue", "bits")
# 
# # Remove certain fields from the taxonomy
# taxo = tidyr::separate(taxo, target,
#                         c("acnu", "gene", "organelle", "strain",
#                           "kingdom", "supergroup", "division", "class", "order", "family", "genus", "species"),
#                         sep="\\|", remove=T)
# taxo$acnu = gsub("\\..*", "", taxo$acnu)
# taxo$gene = NULL
# taxo$organelle = NULL
# taxo$strain = NULL
# taxo$genus = NULL
# 
# # Generating new taxonomic columns
# taxo$group = with(taxo, fifelse(family=="Collodaria_X"|family=="Collophidiidae"|family=="Collosphaeridae"|family=="Sphaerozoidae", "Collodaria",
#                                        fifelse(order=="Nassellaria", "Nassellaria",
#                                                fifelse(order=="Spumellaria", "Spumellaria", class))))
# 
# 
# # Working with the CHARACTERISTICS table ___________________________________________________________
# # Opening characteristics of the OTUs
# char = fread("otus/radiolaria.otus.v20171106_chars.tsv")
# 
# 
# # Working with the RAW ABUNDANCE tables ____________________________________________________________
# # First estimating total eukaryotic reads
# file = fread("globaldataset.otus.v20171106.tsv")
# total = select(file, -c("md5sum", "cid", "ctotab","abundance","sequence","pid","lineage","refs","taxogroup","chloroplast","symb_small","symb_host","silicification","calcification","strontification"))
# total = as.data.frame(colSums(total))
# write.table(total, "totalReads_per_sample.tsv", quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

# # Opening otu table of Radiolaria
# abun_raw = fread("otus/radiolaria.otus.v20171106_table.tsv")
# md5sum_raw = abun_raw$md5sum; abun_raw$md5sum = NULL
# 
# # Adding reads from duplicated columns
# uniques = gsub(";.*", "", colnames(abun_raw))
# uniques = uniques[!duplicated(uniques)]
# abuns = data.frame(md5sum=md5sum_raw)
# for(column in uniques){
#     tmp = rowSums(select(abun_raw, grep(column, names(abun_raw))))
#     abuns = cbind(abuns, tmp)
# };rm(column, tmp)
# 
# colnames(abuns) = c("md5sum", uniques)
# abun_raw = abuns
# rm(uniques, abuns)
# 
# # Get total reads and presence per OTU
# total_raw = data.frame(md5sum=md5sum_raw,
#                     tabun=rowSums(abun_raw[,-1]),
#                     tpres=apply(abun_raw[,-1], 1, function(x) sum(x>0)))
# 
# # merging all tables into a unique RAW table 
# out_raw = merge(taxo, char, by="md5sum", all=FALSE)
# out_raw = merge(out_raw, total_raw, by="md5sum", all=FALSE)
# out_raw = merge(out_raw, abun_raw, by="md5sum", all=FALSE)
# 
# # And exporting
# write.table(out_raw, "data/rad_otus.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
# 
# # Generating fasta files
# out_raw$taxo = with(out_raw, paste0(md5sum, "|", group, "|", family, "|", species, "|id", id, "|abun", tabun, "|pres", tpres))
# seqinr::write.fasta(sequences=as.list(out_raw$sequence), names=out_raw$taxo, nbchar=80, file.out="data/rad_otus.fasta")
# 
# 
# # Working with the MUMU ABUNDANCE tables ___________________________________________________________
# # Opening mumu post-clustered otu table
# abun_mumu = fread("otus/radiolaria_otus_mumu_abun.tsv")
# md5sum_mumu = abun_mumu$md5sum; abun_mumu$md5sum = NULL
# 
# # Adding reads from duplicated columns
# uniques = gsub(";.*", "", colnames(abun_mumu))
# uniques = uniques[!duplicated(uniques)]
# abuns = data.frame(md5sum=md5sum_mumu)
# for(column in uniques){
#     tmp = rowSums(select(abun_mumu, grep(column, names(abun_mumu))))
#     abuns = cbind(abuns, tmp)
# };rm(column, tmp)
# 
# colnames(abuns) = c("md5sum", uniques)
# abun_mumu = abuns
# rm(uniques, abuns)
# 
# # Get total reads and presence per OTU
# total = data.frame(md5sum=md5sum_mumu,
#                     tabun=rowSums(abun_mumu[,-1]),
#                     tpres=apply(abun_mumu[,-1], 1, function(x) sum(x>0)))
# 
# # merging all tables into a unique table 
# out_mumu = merge(taxo, char, by="md5sum", all=FALSE)
# out_mumu = merge(out_mumu, total, by="md5sum", all=FALSE)
# out_mumu = merge(out_mumu, abun_mumu, by="md5sum", all=FALSE)
# 
# # And exporting
# write.table(out_mumu, "data/rad_otus_mumu.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
# 
# # Generating fasta files
# out_mumu$taxo = with(out_mumu, paste0(md5sum, "|", group, "|", family, "|", species, "|id", id, "|abun", tabun, "|pres", tpres))
# seqinr::write.fasta(sequences=as.list(out_mumu$sequence), names=out_mumu$taxo, nbchar=80, file.out="data/rad_otus_mumu.fasta")
# 
# rm(list=ls()[!ls() %in% c()])

#----
#---- Setting names and variables ------------------------------------------------------------------

rm(list=ls()[!ls() %in% c()])
setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/metaB/")

files = list(otusTaxo="data/rad_otus.tsv",
			 env="raw/metadata_assembled_nonRedundant.tsv")

filters = list(id=90,
                tabun=10,
                tpres=2)

# OTU and taxonomic table
data = fread(files$otusTaxo)

if(TRUE){
    cat("Subsetting OTU table\n  - similarity id >=", filters$id,
        "\n  - total abundance >=", filters$tabun,
        "\n  - total presence >= ", filters$tpres,"\n")
    data = subset(data, id >= filters$id & tabun >= filters$tabun & tpres >= filters$tpres)
}

# contextual data
env = fread(files$env)

# Split table into taxonomy and abundances
abun = select(data, grep("TARA", names(data), value=TRUE))
taxo = select(data, grep("TARA", names(data), value=TRUE, invert=TRUE))


#----
#---- Filtering table based on environmental data --------------------------------------------------

envs = env

# select only samples coming from net and pump

# Select only samples coming from SRF, DCM and MES
sort(table(envs$depth))
envs = subset(envs, depth=="SRF" | depth=="DCM" | depth=="MES")

# select only samples of specified size fraction
sort(table(envs$size_fraction))
envs = subset(envs, size_fraction=="0.8-3" | size_fraction=="0.8-5" | size_fraction=="5-20" | size_fraction=="20-180" | size_fraction=="180-2000")
# envs = subset(envs, size_fraction=="0.8-5" | size_fraction=="0.8-3")

sort(table(envs$ocean_region))
sort(table(envs$depth))

# Make depth as a numerical variable
envs$depth_nominal = as.numeric(envs$depth_nominal)

# Subset the table
abuns = select(abun, envs$sample_ID)

#----
#---- Group together OTUs belonging to the same species --------------------------------------------

# bind taxonomic and abundances tables
taxos = select(taxo, c("group", "species"))
taxos$name = paste(taxos$group, taxos$species, sep="|")
taxos$group = NULL
taxos$species = NULL
datas = cbind(taxos, abuns)
 
# Summarise table
datas = datas %>% group_by(name) %>% summarise_all(sum)
 
# And split again
abuns = select(datas, grep("TARA", names(datas), value=TRUE))
taxos = select(datas, grep("TARA", names(datas), value=TRUE, invert=TRUE))

#----
#---- PERMANOVA ------------------------------------------------------------------------------------

file_abun = abuns
file_env = envs

# Transform abundances file
file_abun = as.data.frame(t(file_abun))
# colnames(file_abun) = taxos$name
file_abun = decostand(file_abun, "hellinger")
file_abun = file_abun[, unlist(lapply(file_abun, function(x) sum(x)>0))]

# Check if samples matches
all(rownames(file_abun) == file_env$sample_ID)

# PERMANOVA
permanova = adonis2(file_abun ~ size_fraction, file_env, permutations=1000, method="jaccard")
permanova
cat("Percentage of variability explained by the selected variable is:  ", round(permanova$R2[1]*100, 2), "%\n", sep="")
#
#----
#---- Size fractions -------------------------------------------------------------------------------

file_abun = abun # with OTUs and not taxonomic clusters
file_taxo = taxo
file_env = envs

# Prepare dataset __________________________________________________________________________________
# Select only interesting columns
file_taxo = select(file_taxo, c("md5sum", "group", "species", "id", "tabun", "tpres", "symb_host"))
file_env = select(file_env, c("sample_ID", "station", "depth", "size_fraction", 
                               "biome", "ocean_region", "latitude", "longitude"))

# bind taxonomic and abundances files
file = cbind(file_taxo, file_abun)

# melt the table and remove rows with 0 reads
file = melt(file, id.vars=c("md5sum", "group", "species", "id", "tabun", "tpres", "symb_host"))
colnames(file) = c("md5sum", "group", "species", "id", "tabun", "tpres", "symb_host", "sample_ID", "abundance")
file = file[abundance>0]

# merge with environmental table
file = merge(file, file_env, by="sample_ID")

# Explore size fractions Vs number of OTUs and reads _______________________________________________
# Summarize table for plotting
plot_treeMap_data = file |> group_by(size_fraction, group) |> summarise(OTUs=length(unique(md5sum)), reads=sum(abundance))
plot_treeMap_data = melt(as.data.table(plot_treeMap_data), id.vars=c("size_fraction", "group"))

# Create factors for representation
plot_treeMap_data$group = factor(plot_treeMap_data$group, levels=c("Acantharea", "Spumellaria", "Nassellaria", "Collodaria", "RAD-A", "RAD-B", "RAD-C", "Radiolaria_X"))
plot_treeMap_data$size_fraction = factor(plot_treeMap_data$size_fraction, levels=c("0.8-3", "0.8-5", "5-20", "20-180", "180-2000"))


# And plot
(plot_treeMap = ggplot(plot_treeMap_data, aes(area=value, fill=group))+
    geom_treemap()+
    facet_grid(variable~size_fraction)+
    geom_treemap_text(aes(label= paste(group, value, sep="\n")), place = "centre")+
    theme(legend.position= "none")+
    scale_fill_manual(values=c("yellow3", "steelblue3", "springgreen3", "orangered3", "grey30", "purple", "grey50", "grey70")))

# Exporting
pdf("treeMap_sizeFractions_vs_OTUs&reads.pdf", width=11.69, height=8.27/2, paper='special')
plot(plot_treeMap)
dev.off()

# Get descriptive stats per size fractions _________________________________________________________
# Regarding the number of reads
tmp = plot_treeMap_data
tmp$colonial = ifelse(tmp$group=="Collodaria", "colonial", "solitary")
tmp = subset(tmp, variable=="reads") %>% group_by(size_fraction, colonial) %>% summarise(total=sum(value))
tmp = dcast(as.data.table(tmp), size_fraction ~ colonial)
tmp$total = tmp$colonial + tmp$solitary
tmp$colonialFraction = tmp$colonial/tmp$total
tmp
tmp = subset(plot_treeMap_data, (size_fraction=="0.8-3" | size_fraction=="0.8-5") & variable=="reads")
tmp %>% group_by(size_fraction) %>% mutate(value/sum(value))

# Regarding the number of OTUs
data %>% group_by(group) %>% summarise(length(unique(sequence)))
tmp = subset(plot_treeMap_data, variable=="OTUs")
as.data.frame(tmp %>% group_by(size_fraction) %>% mutate(value/sum(value)*100))

# Explore size fractions Vs ocean region ___________________________________________________________
# Summarize table for plotting
plot_treeMap_data = file |> group_by(size_fraction, ocean_region, group) |> summarise(OTUs=length(unique(md5sum)), reads=sum(abundance))

# Create factors for representation
plot_treeMap_data$group = factor(plot_treeMap_data$group, levels=c("Acantharea", "Spumellaria", "Nassellaria", "Collodaria", "RAD-A", "RAD-B", "RAD-C", "Radiolaria_X"))
plot_treeMap_data$size_fraction = factor(plot_treeMap_data$size_fraction, levels=c("0.8-3", "0.8-5", "5-20", "20-180", "180-2000"))
plot_treeMap_data$ocean_region = factor(plot_treeMap_data$ocean_region, levels=c("MS", "RS", "IO", "SO", "SPO", "NPO", "SAO", "NAO", "AO"))

# Plot number of OTUs
(plot_treeMap = ggplot(plot_treeMap_data, aes(area=OTUs, fill=group))+
        geom_treemap()+
        facet_grid(size_fraction ~ ocean_region)+
        geom_treemap_text(aes(label= paste(group, OTUs, sep="\n")), place = "centre")+
        theme(legend.position= "none")+
        scale_fill_manual(values=c("yellow3", "steelblue3", "springgreen3", "orangered3", "grey30", "purple", "grey50", "grey70")))

# Exporting
pdf("treeMap_sizeFractions_vs_oceanRegion_OTUnumber.pdf", width=11.69, height=8.27/2, paper='special')
plot(plot_treeMap)
dev.off()

# Plot number of READS
(plot_treeMap = ggplot(plot_treeMap_data, aes(area=reads, fill=group))+
        geom_treemap()+
        facet_grid(size_fraction ~ ocean_region)+
        geom_treemap_text(aes(label= paste(group, OTUs, sep="\n")), place = "centre")+
        theme(legend.position= "none")+
        scale_fill_manual(values=c("yellow3", "steelblue3", "springgreen3", "orangered3", "grey30", "purple", "grey50", "grey70")))

# Exporting
pdf("treeMap_sizeFractions_vs_oceanRegion_OTUreads.pdf", width=11.69, height=8.27/2, paper='special')
plot(plot_treeMap)
dev.off()

#----
#---- Estimating contribution of Radiolaria to the total eukaryotic community ----------------------

file_abun = abun # with OTUs and not taxonomic clusters
file_taxo = taxo
file_env = envs

# Total Eukaryotic reads ___________________________________________________________________________
total = fread("raw/totalReads_per_sample.tsv", sep="\t", header=FALSE)
colnames(total) = c("sample_ID", "reads")
total$sample_ID = gsub(";.*", "", total$sample_ID)

# Adding reads from duplicated samples
uniques = unique(total$sample_ID)
tmp = data.frame(sample=c(), reads=c())
for(row in uniques){
    tmp1 = sum(as.numeric(subset(total, sample_ID==row)$reads))
    tmp = rbind(tmp, data.frame(row, tmp1))
};rm(row, tmp1, uniques)
total = tmp; rm(tmp)
colnames(total) = c("sample_ID", "totalReads")

# Calculate total radiolarian reads ________________________________________________________________
rads = as.data.frame(colSums(file_abun))
rads$sample_ID = rownames(rads)
colnames(rads) = c("rads", "sample_ID")

# Merge tables
file = merge(total, rads, by="sample_ID", all.x=FALSE)

# Total proportion _________________________________________________________________________________
file$contribution = file$rads/file$totalReads
# Subset based on environmental data
file = file[file$sample_ID %in% file_env$sample_ID,]

cat("Radiolaria contribute an average ", round(mean(file$contribution*100), 2), 
    "% (+-",round(sd(file$contribution*100), 2), "%) to the total eukaryotic community.\n  With a minimum of ",
	round(min(file$contribution*100), 4), " and a maximum of ", round(max(file$contribution*100), 3), sep="")

summary(file$contribution)
ggplot(file, aes(x=contribution)) + 
    geom_density(fill="steelblue", alpha=0.4) + 
    geom_vline(xintercept = mean(file$contribution), color="red4")


# Plot per station _________________________________________________________________________________
# Prepare dataset
file_taxo = select(file_taxo, c("md5sum", "group", "species", "id", "tabun", "tpres", "symb_host"))
file_env = select(file_env, c("sample_ID", "station", "depth", "size_fraction", 
                               "biome", "ocean_region", "latitude", "longitude"))

# bind taxonomic and abundances files
tmp = cbind(file_taxo, file_abun)

# melt the table and remove rows with 0 reads
tmp = melt(tmp, id.vars=c("md5sum", "group", "species", "id", "tabun", "tpres", "symb_host"))
colnames(tmp) = c("md5sum", "group", "species", "id", "tabun", "tpres", "symb_host", "sample_ID", "abundance")
tmp = tmp[abundance>0]

# merge with environmental table
tmp = merge(tmp, file_env, by="sample_ID")

# And add the contribution
file = merge(file, tmp, by="sample_ID"); rm(tmp)

# Now replace "TARA_148b" by "TARA_148"
file$station = gsub("TARA_148b", "TARA_148", file$station)

# Summaries table for stations, ocean regions, size fractions and depths
plot_data = file |> 
    group_by(ocean_region, station, depth, size_fraction, group) |> 
    summarise(totalReads=sum(totalReads), rads=sum(rads), abundance=sum(abundance))
plot_data$contribution = plot_data$rads/plot_data$totalReads

# Create factors
plot_data$ocean_region = factor(plot_data$ocean_region, levels=c("MS", "RS", "IO", "SO", "SPO", "NPO", "SAO", "NAO", "AO"))
plot_data$station = as.character(gsub("TARA_", "", plot_data$station))
plot_data$depth = factor(plot_data$depth, levels=c("SRF", "DCM", "MES"))
plot_data$size_fraction = factor(plot_data$size_fraction, levels=c("0.8-3", "0.8-5", "5-20", "20-180", "180-2000"))
plot_data$group = factor(plot_data$group, levels=c("Acantharea", "Spumellaria", "Nassellaria", "Collodaria", "RAD-A", "RAD-B", "RAD-C", "Radiolaria_X"))

# And finally the plot
(plot_abun = ggplot(plot_data, aes(x=station, y=abundance, fill=group))+
    geom_bar(position="fill", stat="identity")+
    geom_point(aes(x=station, y=contribution), colour="black")+
    facet_grid(depth+size_fraction~ocean_region, scales="free", space="free")+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))+
    scale_fill_manual(values=c("yellow3", "steelblue3", "springgreen3", "orangered3", "grey30", "purple", "grey50", "grey70")))

# width=11.69, height=8.27
pdf(paste("plot_abundancePerStation_depth_and_sizeFraction.pdf", sep=""), width=11.69*2, height=12, paper='special')
plot(plot_abun)
dev.off()




#----
#---- Prepare ABUNDANCE file for Redundancy Analysis (RDA) -----------------------------------------

# Rename files so we can quickly change the input data without modifying the script ________________
file_abun = abuns

# Transform abundances _____________________________________________________________________________
file_abun = as.data.frame(t(file_abun))
colnames(file_abun) = taxos$name

file_abun = decostand(file_abun, "hellinger")

# Removing 0 columns
file_abun = file_abun[, unlist(lapply(file_abun, function(x) sum(x)>0))]

# Select species that explain more variability by Escoufier's selection criteria ___________________
file_abun_selection = escouf(file_abun, level=0.95)
plot(file_abun_selection)
file_abun_selection$level = 0.90 # We choose a treshold of 90%
lines(file_abun_selection)
file_abun = extract(file_abun_selection)

# species_contribution = c()
# for(i in 1:length(file_abun_selection$RV)){
#     if(i == 1){species_contribution[i] = file_abun_selection$RV[[i]]
#     }else{species_contribution[i] = file_abun_selection$RV[[i]] - file_abun_selection$RV[[i-1]]}
# }; rm(i)
# species_contribution = setNames(species_contribution, names(file_abun_selection$RV))

# Remove empty rows
file_abun = file_abun[rowSums(file_abun)!=0,]

#---- Prepare ENVIRONMENTAL file for Redundancy Analysis (RDA) -------------------------------------

# Rename files so we can quickly change the input data without modifying the script ________________
file_env = envs

# Match order with abundance file
file_env = file_env[file_env$sample_ID %in% row.names(file_abun)]
file_env = as.data.frame(file_env[order(row.names(file_abun)),])
rownames(file_env) = file_env$sample_ID

if(all(rownames(file_abun) == rownames(file_env))){cat("All good! :)")}else{cat("Warning! abundances and environmental files have different row ordering")}

# Safe qualitative variables _______________________________________________________________________
file_env_qual = file_env[, !unlist(lapply(file_env, is.numeric))]

# Extracting only numerical variables ______________________________________________________________
file_env = file_env[, unlist(lapply(file_env, is.numeric))]

# Have a look at the distribution of the environmental variables ___________________________________
file_env_long = data.frame(variable=colnames(file_env), skewness=round(apply(file_env, 2, skewness), 2))
file_env_long = merge(melt(as.data.table(file_env), measure=1:ncol(file_env)), file_env_long, by="variable", sort=FALSE)

# (skew_plot = ggplot(file_env_long) +
#         geom_histogram(aes(x=value), fill="steelblue2", alpha=0.4)+
#         facet_wrap(~variable, scales = "free")+
#         theme_bw())

# Transform environmental variables ________________________________________________________________
# # Apply a log transformation
file_env = decostand(file_env, "log", 2)

# # Apply a standardizarion
# file_env = decostand(file_env, "standardize", 2)

# # Or instead, transform only parameters highly skewed ____________________________________________
# skew_threshold = 3
# for(i in 1:length(unique(skew$variable))){
#     if(abs(skew$skewness[i]) > skew_threshold){
#         file_env[,i] = decostand(file_env[,i], "log")
#     }
# };rm(i)

# Replacing NAs for mean values ____________________________________________________________________
for (i in 1:ncol(file_env)) {
    file_env[is.na(file_env[,i]),i]=mean(file_env[,i], na.rm=TRUE)
}; rm(i)

if(!any(apply(file_env, 2, function(x) any(is.na(x))))){cat("All good! :)")}else{cat("Warning! There are NAs in your contextual data file")}

# Adding the filter as centroids ___________________________________________________________________
# file_env$FilterA = file_env_qual$size_fraction
file_env$FilterB = file_env_qual$depth


#----
#---- Redundancy Analysis (RDA) --------------------------------------------------------------------

rda_data = rda(file_abun~., file_env)

# summary(rda_data)
RsquareAdj(rda_data)

# We use the ordistep function to select variables through permutation tests :
(seed = sample(1:10^9, 1))
TeachingDemos::char2seed(seed, set=TRUE) # Ordistep is a function involving randomness so we set a seed for reproductability
rda_both = ordistep(rda(file_abun~1, data=file_env), scope=formula(rda_data), direction="both", 
                     steps=10^8, # steps=10^8, pstep=5000
					permutations = 10^4) # permutations = 10^4

formula(rda_both)

rda_data=rda(formula=formula(rda_both), data=file_env, na.action=na.omit)
# summary(rda_data)
RsquareAdj(rda_data)

rda_anova = anova.cca(rda_data, by="axis", model="direct", step=1000, cutoff=1, parallel = getOption("mc.cores"))

# Get stats ________________________________________________________________________________________
cat("The total ", nrow(abun), " metabarcodes were clustered into ", nrow(abuns), " unique lineages,\n",
	"  of which ", ncol(file_abun), " were selected by Escoufier's method of equivalent vectors as the most significant lineages.", sep="")
tmp = as.data.frame(scores(rda_data, display="sp", scaling=2))
tmp$contribution = sqrt(tmp$RDA1^2 + tmp$RDA2^2)
tmp[order(tmp$contribution),]

cat("RDA R²:", round(RsquareAdj(rda_data)$adj.r.squared*100, 2), "% (unadjusted R²: ",round(RsquareAdj(rda_data)$r.squared*100, 2),"%)\n", 
	"  Variability explained by the first two axis: ", round(summary(rda_data)$cont$importance[2,1]*100+summary(rda_data)$cont$importance[2,2]*100, 2), "%\n",
	"    Axis 1 (RDA1): ", round(summary(rda_data)$cont$importance[2,1]*100, 2),"%\n",
	"    Axis 2 (RDA2): ", round(summary(rda_data)$cont$importance[2,2]*100, 2),"%\n", sep="")
tmp = as.data.frame(scores(rda_data, choices=1:2, display="bp", scaling=2))
tmp[order(tmp[,1]),]
tmp[order(tmp[,2]),]



#

#---- Prepare files for plotting Redundancy Analysis (RDA) -----------------------------------------

# Create samples coordinates table _________________________________________________________________
plot_samples = as.data.frame(scores(rda_data, scaling=2)$sites)
plot_samples$sample_ID = rownames(plot_samples)
plot_samples = merge(plot_samples, 
                 select(file_env_qual, c("sample_ID", "station", "depth", "size_fraction", "biome", "ocean_region")), 
                 by="sample_ID", all.x=T, all.y=F)

# Adding environmental variables ___________________________________________________________________
plot_env = as.data.frame(scores(rda_data, choices = 1:2, display="bp", scaling=2))
# we get rid of all qualitative variables (the centroids we previously added)
plot_env = subset(plot_env, !grepl("Filter", rownames(plot_env)))

# # Get only variables that explains the most (optional)
# plot_env = subset(plot_env, abs(RDA1)>=0.01 & abs(RDA2)>=0.01)
plot_env = subset(plot_env, sqrt(RDA1^2 + RDA2^2) >= 0.01)

# Create species coordinates table _________________________________________________________________
plot_species = as.data.frame(scores(rda_data, display="sp", scaling=2))
plot_species$group = gsub("\\|.*","", rownames(plot_species))
plot_species$species = gsub(".*\\|","", rownames(plot_species))

# # Keep only properly represented species in the triplot (optional)
# plot_species = subset(plot_species, abs(RDA1)>=0.01 & abs(RDA2)>=0.01)
plot_species = subset(plot_species, sqrt(RDA1^2 + RDA2^2) >= 0.01)

# add colors for the species
plot_species$colour = plot_species$group
{plot_species$colour[which(plot_species$colour=="Nassellaria")]="springgreen3"
    plot_species$colour[which(plot_species$colour=="Acantharea")]="yellow3"
    plot_species$colour[which(plot_species$colour=="Spumellaria")]="steelblue3"
    plot_species$colour[which(plot_species$colour=="Collodaria")]="orangered3"
    plot_species$colour[which(plot_species$colour=="RAD-A")]="grey40"
    plot_species$colour[which(plot_species$colour=="RAD-B")]="purple"
    plot_species$colour[which(plot_species$colour=="RAD-C")]="grey60"
    plot_species$colour[which(plot_species$colour=="Radiolaria_X")]="grey80"}

# Adding centroids for the filter and depth factors ________________________________________________
plot_filter = as.data.frame(rda_data$CCA$centroids[,1:2])
plot_filter$filter = gsub("Filter.","",rownames(plot_filter))
# plot_filter$filter = factor(plot_filter$filter, levels=c("0.8-5", "5-20", "20-180", "180-2000"))
plot_filter$filter = factor(plot_filter$filter, levels=c("SRF", "DCM", "MES"))


#---- Plotting Redundancy Analysis (RDA) -----------------------------------------------------------

# First we plot the samples, variability explained by represented variables and some general aesthetics
(rda_plot = ggplot() + 
        geom_point(data=plot_samples ,aes(RDA1, RDA2), color="grey80")+
        xlab(paste0("RDA1 (", round(summary(rda_data)$cont$importance[2,1]*100, 2),"%)")) +  
        ylab(paste0("RDA2 (", round(summary(rda_data)$cont$importance[2,2]*100, 2),"%)")) +
        labs(title=paste0("RDA R²: ", round(RsquareAdj(rda_data)$adj.r.squared*100, 2), 
                          " (unadjusted R²: ",round(RsquareAdj(rda_data)$r.squared*100, 2),")"))+
        geom_hline(yintercept = 0, linetype='dotted') +
        geom_vline(xintercept = 0, linetype='dotted') +
        theme(plot.title=element_text(hjust=0.5))+
        theme_minimal())

# Then the environmental variables
factor = 1 # 1.5 we multiply the arrows by a factor for a better visualization (optional)
(rda_plot = rda_plot + 
        geom_segment(data=plot_env, aes(xend=RDA1*factor, yend= RDA2*factor),
                     x=0, y=0, linewidth= 0.4, color= 'steelblue4', arrow= arrow(length= unit(0.2,"cm"))) +
        geom_text(data= plot_env, aes(RDA1*factor, RDA2*factor, 
                      label=rownames(plot_env)), color="steelblue4", size=2.5))

# Now we plot the species
(rda_plot = rda_plot  + geom_segment(data=plot_species, aes(xend=RDA1, yend=RDA2), x=0, y=0,
                          arrow=arrow(length=unit(0.2,"cm")), size=0.4, color=plot_species$colour)+
        geom_text(data=plot_species ,aes(RDA1, RDA2, label=species), 
                  color=plot_species$colour, size=2))

# And finally we add the centroids (if we had selected them)
(rda_plot = rda_plot + 
        geom_point(data=plot_filter, aes(x=RDA1, y=RDA2, shape=filter), color="darkorange2", size=4) +
        scale_shape_manual(values = c(21,22,24)))

#---- Export the plot into a pdf for further aesthetics modification in Inkscape -------------------
# filters
if(any(grepl("^RDA.*\\.pdf", dir()))){
	tmp = grep("^RDA.*\\.pdf", dir(), value=TRUE)
	tmp = tmp %>% sub("\\.pdf$", "", .) %>% sub(".*\\-", "", .)
	tmp = max(as.numeric(tmp))
	name = paste0("RDA_lineages_id", filters$id, "_reads", filters$tabun, "_pres", filters$tpres, 
				  "_Hellinger_envLog_pico-nano_steps108_perm104_", 
				  format(Sys.Date(), "%y%m%d"), "-", tmp+1, ".pdf")
}else{
	name = paste0("RDA_lineages_id", filters$id, "_reads", filters$tabun, "_pres", filters$tpres, 
				  "_Hellinger_envLog_pico-nano_steps108_perm104_", 
				  format(Sys.Date(), "%y%m%d"), "-1.pdf")
}
name
pdf(name, width=11.69, height=8.27, paper='special')
plot(rda_plot)
dev.off()

# save.image(file="RDA_speciesHellinger_envLog_id95_reads10_pres2.Rdata")

# Now plotting the stations by the environmental parameters ________________________________________

# (tmp = ggplot(plot_samples, aes(x=RDA1, y=RDA2, label=gsub("TARA_","",station), colour=ocean_region)) + 
#         geom_text(aes())+
#         theme_minimal())

# pdf(paste("station_names.pdf", sep=""), width=11.69, height=8.27, paper='special')
# plot(tmp)
# dev.off()


#---- 
load(file="RDA_speciesHellinger_envLog_id95_reads10_pres2.Rdata")
#----
#---- Distance-based Redundancy Analysis (dbRDA) ---------------------------------------------------
# dbrda.files = capscale(file_abuns~., file_envs, distance="jaccard")
# 
# summary(dbrda.files)
# RsquareAdj(dbrda.files)
# # Adj Rsq = 29.84% (26.43%) for all size fractions
# 
# # We use the ordistep function to select variables through permutation tests :
# TeachingDemos::char2seed('sandin', set=T) # Ordistep is a function involving randomness so we set a seed for reproductability
# dbrda.step.both = ordistep(capscale(file_abuns~1, data=file_envs), scope=formula(dbrda.files), 
#                           distance="jaccard", direction="both", steps=5000)  # pstep=5000)
# 
# formula(dbrda.step.both)
# 
# dbrda.files=rda(formula=formula(rda.step.both), data=file_envs, na.action=na.omit)
# summary(dbrda.files)
# RsquareAdj(dbrda.files)
# Adj Rsq = 26.71% (24.61%) for all size fractions

#----
