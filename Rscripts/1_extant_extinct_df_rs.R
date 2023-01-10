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
library(tictoc)
library(profvis)
library(stats)
library(reshape2)

rm(list = ls())
if(!is.null(dev.list())) dev.off()
base_dir <- "C:/Users/tomos/OneDrive/PhD/ExHosts_PathRichness/"
species <- list.files(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural"), pattern = "tif$")
Phylacline_trait <- read.csv(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Traits/Trait_data.csv"))
Phylacline_trait <- data.frame(Phylacline_trait)
Extinct_species <- Phylacline_trait %>% filter(IUCN.Status.1.2 == 'EP' | IUCN.Status.1.2 == 'EX' | IUCN.Status.1.2 == 'EW')
Extinct_species <- Extinct_species %>% filter(Terrestrial == 1 | Aerial == 1)
Extant_species <- Phylacline_trait %>% filter(IUCN.Status.1.2 != 'EP' & IUCN.Status.1.2 != 'EX' & IUCN.Status.1.2 != 'EW')
Extant_species <- Extant_species %>% filter(Terrestrial == 1 | Aerial == 1)
Extant_species <- Extant_species %>% filter(Marine != 1)
Extant_species <- Extant_species %>% filter(Binomial.1.2 != 'Homo_sapiens')
no_extinct_sp <- nrow(Extinct_species)
no_extant_sp <- nrow(Extant_species)

# Nexus file contains 1000 trees. The simplify argument means that only the first tree is used.
# Eventually needs modifying to use all trees, or the consensus tree (simple_tree.nex in this directory)
phylo_tree <- readNexus(file = paste0(base_dir, "PHYLACINE_v1.2.1/Data/Phylogenies/Complete_phylogeny.nex"), simplify = TRUE, type = c("all"))
phylo_tree <- as(phylo_tree, "phylo")
#distNodes(phylo_tree, node = c("Bison_bison","Loxodonta_africana"), clus=0.5)
distance_matrix <- cophenetic.phylo(phylo_tree)
saveRDS(distance_matrix, file = paste0(base_dir,'intermediates/distance_matrix.rda'))
#pairwise <- distTips(phylo_tree, tips = c("Bison_bison","Loxodonta_africana"), method = c("patristic"), useC = TRUE)
#Distance_matrix <- readDist("Complete_phylogeny.nex", format = "nexus")

### Create the mean and median phylogenetic distance matrix ###
phylo_tree2 <- readNexus(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Phylogenies/Complete_phylogeny.nex"), simplify = F, type = c("all"))
#phylo_tree2 <- as(phylo_tree2, "phylo")

mean_distance_matrix <- distance_matrix
mean_distance_matrix[,] <- 0
for(i in 1:1000){
  print(i)
  test_phylo_tree <- phylo_tree2[[i]]
  test_phylo_tree <- as(test_phylo_tree, "phylo")
  test_distance_matrix <- cophenetic.phylo(test_phylo_tree)
  mean_distance_matrix <- mean_distance_matrix + test_distance_matrix
}
mean_distance_matrix <- mean_distance_matrix / 1000
saveRDS(mean_distance_matrix, file = paste0(base_dir,'intermediates/mean_distance_matrix.rda'))

### Read in phylogeny based on Phylacine github:
#phylo_tree2 <- readNexus(file = paste0(base_dir, "PHYLACINE_v1.2.1/Data/Phylogenies/Complete_phylogeny.nex"))
#names(phylo_tree2) <- NULL
#set.seed(42)
#phylo_tree2 <- phylo_tree2[sample(1:1000, 1000)]
#consensus_tree <- consensus(phylo_tree2, p = 1, check.labels = TRUE, rooted = TRUE)

###############################################################
#### Generate location dfs for extinct and extant mammals #####
###############################################################

### Generate data frames of all the raster columns, rows and cell values where each species is located, both extinct and extant:
extinct_df <- data.frame(species = character(), row = numeric(), col = numeric(), ras_cell_value = numeric())
for(jj in 1:no_extinct_sp){ 
  #jj <- 1
  #tic("sleeping")
  focal_name <- Extinct_species[jj,1]
  print(paste0(focal_name, "; number ", jj, " of ", no_extinct_sp))
  focal_species <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", focal_name,".tif"))
  presence_coords <- which(values(focal_species) == 1)
  presence_rows <- rowFromCell(focal_species, presence_coords)
  presence_cols <- colFromCell(focal_species, presence_coords)
  name_rep <- rep(focal_name, length(presence_coords))
  focal_entry <- cbind(name_rep, presence_rows, presence_cols, presence_coords)
  extinct_df <- rbind(extinct_df, focal_entry)
  #toc()
}
names(extinct_df) <- c("species", "row", "col", "ras_cell_value")
saveRDS(extinct_df, file = paste0(base_dir,'intermediates/extinct_df.rda'))

extant_df <- data.frame(species = character(), row = numeric(), col = numeric(), ras_cell_value = numeric())
for(jj in 1:no_extant_sp){ 
  #tic("sleeping")
  focal_name <- Extant_species[jj,1]
  print(paste0(focal_name, "; number ", jj, " of ", no_extant_sp))
  focal_species <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", focal_name,".tif"))
  presence_coords <- which(values(focal_species) == 1)
  presence_rows <- rowFromCell(focal_species, presence_coords)
  presence_cols <- colFromCell(focal_species, presence_coords)
  name_rep <- rep(focal_name, length(presence_coords))
  focal_entry <- cbind(name_rep, presence_rows, presence_cols, presence_coords)
  extant_df <- rbind(extant_df, focal_entry)
  #toc()
}
names(extant_df) <- c("species", "row", "col", "ras_cell_value")
saveRDS(extant_df, file = paste0(base_dir,'intermediates/extant_df.rda'))


###############################################################
#### Generate raster stacks for extinct and extant mammals ####
###############################################################

extinct_raster_stack <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", Extinct_species[1,1],".tif")) 
for(jj in 2:no_extinct_sp){ 
  #tic("sleeping")
  #print(jj)
  temp_species <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", Extinct_species[jj,1],".tif"))
  extinct_raster_stack <- addLayer(extinct_raster_stack, temp_species)
  #toc()
}
saveRDS(extinct_raster_stack, file = paste0(base_dir,'intermediates/extinct_rs.rda'))

extant_raster_stack <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", Extant_species[1,1],".tif")) 
for(jj in 2:no_extant_sp){ 
  #tic("sleeping")
  #print(jj)
  temp_species <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", Extant_species[jj,1],".tif"))
  extant_raster_stack <- addLayer(extant_raster_stack, temp_species)
  #toc()
}
saveRDS(extant_raster_stack, file = paste0(base_dir,'intermediates/extant_rs.rda'))


#Cumulative data frame showing where all extinct AND extant species are located
all_sp_df <- rbind(extant_df, extinct_df)
cumulative_sp_df <- all_sp_df %>% 
  group_by(ras_cell_value, row, col) %>%
  dplyr::summarise(cumulative_sp = n())
saveRDS(cumulative_sp_df, file = paste0(base_dir,'intermediates/cumulative_sp_df.rda'))

# Land area data frame, showing all terrestrial land (assuming that at least a single Phylacine mammal is found there)
land_area_df <- cumulative_sp_df
names(land_area_df)[names(land_area_df) == 'cumulative_sp'] <- 'land_area'
land_area_df$land_area <- 1
saveRDS(land_area_df, file = paste0(base_dir,'intermediates/landarea_df.rda'))

# Data frame showing the cumulative number of extant terrestrial mammals in each location
all_extant_sp_df <- extant_df %>%
  group_by(ras_cell_value, row, col) %>%
  dplyr::summarise(cumulative_sp = n())
saveRDS(all_extant_sp_df, file = paste0(base_dir,'intermediates/cumulative_extant_df.rda'))

# Data frame showing the cumulative number of extinct terrestrial mammals in each location
all_extinct_sp_df <- extinct_df %>%
  group_by(ras_cell_value, row, col) %>%
  dplyr::summarise(cumulative_sp = n())
saveRDS(all_extinct_sp_df, file = paste0(base_dir,'intermediates/cumulative_extinct_df.rda'))





### Generate a raster of land area and the combination of all extinct and extant species:
cumulative_sp_df <- readRDS(file = paste0(base_dir,'intermediates/cumulative_sp_df.rda')) 
cumulative_sp_raster <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", Extinct_species[1,1],".tif"))
names(cumulative_sp_raster) <- "Extinct_and_extant_sp"
cumulative_sp_raster[] <- 0

for(q in 1:length(cumulative_sp_df$ras_cell_value)){
  #q <- 1
  #print(q)
  focal_coords <- c(as.numeric(cumulative_sp_df$row[q]), as.numeric(cumulative_sp_df$col[q]))
  #print(focal_coords)
  focal_x <- xFromCol(cumulative_sp_raster, as.numeric(cumulative_sp_df$col[q]))
  focal_y <- yFromRow(cumulative_sp_raster, as.numeric(cumulative_sp_df$row[q]))
  focal_xy <- c(focal_x, focal_y)
  raster_cell_value <- cellFromXY(cumulative_sp_raster, focal_xy)
  #print(raster_cell_value)
  cumulative_sp_raster[raster_cell_value] <- cumulative_sp_df$cumulative_sp[q]
}

cumulative_sp_raster
image(cumulative_sp_raster)
saveRDS(cumulative_sp_raster , file = paste0(base_dir,'intermediates/cumulative_sp_raster.rda')) 
land_area_raster <- cumulative_sp_raster
land_area_raster[land_area_raster > 0] <- 1
image(land_area_raster)
saveRDS(land_area_raster, file = paste0(base_dir,'intermediates/land_area_raster.rda')) 


### Generate a raster for each of the extinct and extant species:
### Extinct:
cumulative_extinct_df <- readRDS(file = paste0(base_dir,'intermediates/cumulative_extinct_df.rda')) 
cumulative_extinct_raster <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", Extinct_species[1,1],".tif"))
names(cumulative_extinct_raster) <- "All_extinct_sp"
cumulative_sp_raster[] <- 0
for(q in 1:length(cumulative_extinct_df$ras_cell_value)){
  #print(q)
  focal_x <- xFromCol(cumulative_extinct_raster, as.numeric(cumulative_extinct_df$col[q]))
  focal_y <- yFromRow(cumulative_extinct_raster, as.numeric(cumulative_extinct_df$row[q]))
  focal_xy <- c(focal_x, focal_y)
  raster_cell_value <- cellFromXY(cumulative_extinct_raster, focal_xy)
  cumulative_extinct_raster[raster_cell_value] <- cumulative_extinct_df$cumulative_sp[q]
}
cumulative_extinct_raster
image(cumulative_extinct_raster)
saveRDS(cumulative_extinct_raster , file = paste0(base_dir,'intermediates/cumulative_extinct_raster.rda')) 


### Extant:
cumulative_extant_df <- readRDS(file = paste0(base_dir,'intermediates/cumulative_extant_df.rda')) 
cumulative_extant_raster <- raster(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Ranges/Present_natural/", Extant_species[1,1],".tif"))
names(cumulative_extant_raster) <- "All_extant_sp"
cumulative_extant_raster[] <- 0
for(q in 1:length(cumulative_extant_df$ras_cell_value)){
  #print(q)
  focal_x <- xFromCol(cumulative_extant_raster, as.numeric(cumulative_extant_df$col[q]))
  focal_y <- yFromRow(cumulative_extant_raster, as.numeric(cumulative_extant_df$row[q]))
  focal_xy <- c(focal_x, focal_y)
  raster_cell_value <- cellFromXY(cumulative_extant_raster, focal_xy)
  cumulative_extant_raster[raster_cell_value] <- cumulative_extant_df$cumulative_sp[q]
}
cumulative_extant_raster
image(cumulative_extant_raster)
saveRDS(cumulative_extant_raster , file = paste0(base_dir,'intermediates/cumulative_extant_raster.rda')) 

### Generate a matrix for each of the extinct and extant species using the all_extinct_sp_df and all_extant_sp_df:

### Extinct:
# Generate a matrix of zeros with the correct dimensions:
extinct_matrix <- array(rep(0, 360*142*no_extinct_sp), dim=c(360,142,no_extinct_sp))
# Populate
for(i in 1:no_extinct_sp){
  #i <- 2
  print(i)
  focal_species <- Extinct_species$Binomial.1.2[i]
  #print(focal_species)
  df_subset <- extinct_df[extinct_df$species == focal_species,]
  for(j in 1:length(df_subset$species)){
    extinct_matrix[as.numeric(df_subset$col[j]), (142 - as.numeric(df_subset$row[j])), i] <- 1
  }
}
saveRDS(extinct_matrix, file = paste0(base_dir,'intermediates/extinct_matrix2.rda'))
cumulative_matrix_extinct <- rowSums(extinct_matrix, dims = 2)

### Extant:
# Generate a matrix of zeros with the correct dimensions:
extant_matrix <- array(rep(0, 360*142*no_extant_sp), dim=c(360, 142, no_extant_sp))
# Populate
for(i in 1:no_extant_sp){
  #i <- 2
  print(i)
  focal_species <- Extant_species$Binomial.1.2[i]
  #print(focal_species)
  df_subset <- extant_df[extant_df$species == focal_species,]
  for(j in 1:length(df_subset$species)){
    extant_matrix[as.numeric(df_subset$col[j]), (142 - as.numeric(df_subset$row[j])), i] <- 1
  }
}
saveRDS(extant_matrix, file = paste0(base_dir,'intermediates/extant_matrix2.rda'))
cumulative_matrix_extant <- rowSums(extant_matrix, dims = 2)
