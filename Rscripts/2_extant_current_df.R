library(readr)
library(dplyr)
library(MASS)
library(car)
library(ellipse)
library(ggplot2)
library(psych)
library(GPArotation)
library(Matrix)
library(MVN)
library(ecodist)
library(ape)
library(tiff)
library(raster)
library(gdalUtilities)
library(rgdal)
library(rgeos)
library(mapview)
library(tidyverse)
library(maps)
library(phonTools)
library(phangorn)
library(phylobase)
library(adephylo)
library(oceanmap)
library(abind)
library(devtools)
library(tictoc)
library(profvis)
library(stats)
library(reshape2)

rm(list = ls())
if(!is.null(dev.list())) dev.off()
base_dir <- "C:/Users/tomos/OneDrive/PhD/ExtinctHosts_PathogenRichness/"
species <- list.files(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Current"), pattern = "tif$")
Phylacine_trait <- read.csv(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Traits/Trait_data.csv"))
Phylacine_trait <- data.frame(Phylacine_trait)
Extant_species <- Phylacine_trait %>% filter(IUCN.Status.1.2 != 'EP' & IUCN.Status.1.2 != 'EX' & IUCN.Status.1.2 != 'EW')
Extant_species <- Extant_species %>% filter(Terrestrial == 1 | Aerial == 1)
Extant_species <- Extant_species %>% filter(Marine != 1)
Extant_species <- Extant_species %>% filter(Binomial.1.2 != 'Homo_sapiens')
no_extant_sp <- nrow(Extant_species)

###############################################################
#### Generate location dfs for extinct and extant mammals #####
###############################################################

### Generate data frames of all the raster columns, rows and cell values where each species is located:

extant_current_df <- data.frame(species = character(), row = numeric(), col = numeric(), ras_cell_value = numeric())

for(jj in 1:no_extant_sp){ 
  #tic("sleeping")
  focal_name <- Extant_species[jj,1]
  print(paste0(focal_name, "; number ", jj, " of ", no_extant_sp))
  focal_species <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Current/", focal_name,".tif"))
  presence_coords <- which(values(focal_species) == 1)
  presence_rows <- rowFromCell(focal_species, presence_coords)
  presence_cols <- colFromCell(focal_species, presence_coords)
  name_rep <- rep(focal_name, length(presence_coords))
  focal_entry <- cbind(name_rep, presence_rows, presence_cols, presence_coords)
  extant_current_df <- rbind(extant_current_df, focal_entry)
  #toc()
}
names(extant_current_df) <- c("species", "row", "col", "ras_cell_value")
saveRDS(extant_current_df, file = paste0(base_dir,'intermediates/extant_current_df.rda'))

