#----
#---- Libraries ------------------------------------------------------------------------------------

library(data.table)

#----
#---- Open the relevant files ----------------------------------------------------------------------

setwd("~/Desktop/PhD/0_Thesis/2_chapter/data/metaB/")

# OTU and taxonomic table
samples <- fread("data/rad_otus.tsv")

# contextual data
env <- fread("raw/TARA_registies_mesoscale.tsv")

# Contextual data to complement

env1 <- fread("raw/TARA_sample_enviro.tsv")
env2 <- fread("raw/TARA_SAMPLES_CONTEXT_METHODS_20180320.tsv")

#----
#---- Extract only samples where we have biological material ---------------------------------------

samples <- grep("TARA_", names(samples), value=TRUE)

env <- env[env$`Sample ID (TARA_barcode#)` %in% samples,]

env1 <- env1[env1$`Sample ID (TARA_barcode#)` %in% samples,]
env2 <- env2[env2$`TARA_barcode#` %in% samples,]

#----
#---- Modify names and qualitative values for simplicity -------------------------------------------

colnames(env) <- c("sample_ID",
                   "campaign",
                   "station",
                   "event",
                   "date_time",
                   "latitude",
                   "longitude",
                   "depth",
                   "depth_nominal",
                   "sample_method",
                   "biome",
                   "ocean_region",
                   "biogeographic_province",
                   "fluorescence", # fCDOM [ppb (QSE)]
                   "PAR",
                   "PAR_mean8days",
                   "PAR_mean30days",
                   "NPP_C_mean8days",
                   "NPP_C_mean30days",
                   "moon_phase_nominal",
                   "moon_phase_proportion",
                   "SSD", # Sunshine duration
                   "season",
                   "season_period",
                   "bathymetry",
                   "latitude_closestCoast",
                   "longitude_closestCoast",
                   "coast_distance",
                   "SSM", # Strain sub-mesoscale index
                   "SST",
                   "geostrophic_velocity_longitudinal", # u [cm/s]
                   "geostrophic_velocity_latitudinal", # v [cm/s]
                   "Okubo_Weiss",
                   "Lyapunov",
                   "residence_time",
                   "latitude_continentalShelf",
                   "longitude_continentalShelf",
                   "age",
                   "latitude_origin",
                   "longitude_origin")


env$event <- gsub(".*_", "", env$event)

env$depth <- env$depth %>% gsub("\\].*", "", .) %>% gsub("\\[", "", .)

env$ocean_region <- env$ocean_region %>% gsub("\\[O\\] North Atlantic Ocean \\(MRGID:1912\\)", "NAO", .) %>% 
    gsub("\\].*", "", .) %>% gsub("\\[", "", .)

env$biogeographic_province <- env$biogeographic_province %>% gsub("\\[.*\\] ", "", .) %>% gsub(" \\(.*\\)", "", .)

#----
#---- Merge with other contextual tables to grasp for other parameters -----------------------------

# Keep only columns of interest from the secondary tables and merge with main table
env1 <- select(env1, c("Sample ID (TARA_barcode#)", 
                       "Sal (in the selected environmental...)", 
                       "Tpot [°C] (in the selected environmental...)", 
                       "O2 [µmol/kg] (in the selected environmental...)", 
                       "[NO3]- [µmol/l] (in the selected environmental...)", 
                       "[NO2]- [µmol/l] (in the selected environmental...)", 
                       "[NO3]- + [NO2]- [µmol/l] (in the selected environmental...)", 
                       "[PO4]3- [µmol/l] (in the selected environmental...)", 
                       "Si(OH)4 [µmol/l] (in the selected environmental...)", 
                       "Chl a [mg/m**3] (in the selected environmental...)", 
                       "beta470 [m/sr] (in the selected environmental...)",  
                       "bb470 [1/m] (in the selected environmental...)", 
                       "bbp470 [1/m] (in the selected environmental...)", 
                       "bac660 [1/m] (in the selected environmental...)", 
                       "bacp [1/m] (in the selected environmental...)", 
                       "MLD [m] (in the selected environmental...)", 
                       "D chl m [m] (in the selected environmental...)", 
                       "Depth max Brunt Väisälä freq [m] (in the selected environmental...)", 
                       "Depth max O2 [m] (in the selected environmental...)", 
                       "Depth min O2 [m] (in the selected environmental...)", 
                       "Depth nitracline [m] (in the selected environmental...)"))
colnames(env1) <- c("sample_ID", 
                    "salinity", 
                    "temperature", 
                    "O2", 
                    "NO3", 
                    "NO2", 
                    "NO3NO2", 
                    "PO43", 
                    "SiOH4", 
                    "ChlA", 
                    "angular_scattering_coef", # beta470
                    "optical_backscattering_coef", # bb470 
                    "particle_backscattering_coef", # bbp470
                    "optical_attenuation_coef", # bac660 
                    "particle_optical_attenuation_coef", # bacp
                    "depth_mix_layer", 
                    "depth_ChlA_max", 
                    "depth_max_Brunt", 
                    "depth_max_O2", 
                    "depth_min_O2", 
                    "depth_nitracline")

env2 <- select(env2, c("TARA_barcode#", 
                       "Size fraction, lower threshold", 
                       "Size fraction, upper threshold"))
colnames(env2) <- c("sample_ID", 
                    "sf_low", 
                    "sf_big")

env2$size_fraction <- paste(env2$sf_low, env2$sf_big, sep="-")
env2$sf_low <- NULL
env2$sf_big <- NULL

# And finally merge
env <- merge(env, env1, by="sample_ID", sort=FALSE)
env <- merge(env, env2, by="sample_ID", sort=FALSE)


#----
#---- Export tables --------------------------------------------------------------------------------

write.table(env, "raw/metadata_assembled.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Keep only non-redudant columns and order them
envs <- select(env, c("sample_ID", "campaign", "station",
                      "event",
                      "date_time",
                      "latitude","longitude",
                      "depth",
                      "sample_method",
                      "biome", "ocean_region", "biogeographic_province",
                      "size_fraction",
                      "depth_nominal",
                      "temperature", 
                      "SST",
                      "O2", "NO3", "NO2", "NO3NO2", "PO43", "SiOH4", 
                      "ChlA", 
                      "fluorescence", # fCDOM [ppb (QSE)]
                      "PAR",
                      "NPP_C_mean8days",
                      "SSD", # Sunshine duration
                      "bathymetry",
                      "coast_distance",
                      "Okubo_Weiss",
                      "Lyapunov",
                      "residence_time",
                      "age", 
                      "salinity", 
                      "angular_scattering_coef", # beta470
                      "optical_backscattering_coef", # bb470 
                      "particle_backscattering_coef", # bbp470
                      "optical_attenuation_coef", # bac660 
                      "particle_optical_attenuation_coef", # bacp
                      "depth_mix_layer", "depth_ChlA_max", "depth_max_Brunt", "depth_max_O2", "depth_min_O2", "depth_nitracline"))

write.table(envs, "raw/metadata_assembled_nonRedundant.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#----
#----
