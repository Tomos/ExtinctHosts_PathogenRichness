library(tidyverse)
library(stringi)
library(sp)
library(rgdal)
library(raster)
library(classInt)
library(rasterVis)
library(maptools)
library(mgcv)
library(parallel)
library(rgeos)
library(gridExtra)
library(RColorBrewer)
library(latticeExtra)
library(terra)
library(tictoc)
library(sf)
library(reshape2)
library(rnaturalearth)
library(inlmisc)
library(ggplot2)
library(cartomisc)
library(scico)
library(vign)

rm(list = ls())
if(!is.null(dev.list())) dev.off()
base_dir <- "C:/Users/tomos/OneDrive/PhD/Host-parasite-paper/"
map_dir <- paste0(base_dir, "maps/")

blankPlot <- ggplot()+geom_blank(aes(1,1)) +
  cowplot::theme_nothing()

map_type <- "missing" #predicted_max
noExtinction_gmpd2_raster <- readRDS(paste0(base_dir, "intermediates/noExtinction-GMPD2_helminth_rasters_ALL.rds"))
Extinction_gmpd2_raster <- readRDS(paste0(base_dir, "intermediates/Extinction-GMPD2_helminth_rasters_ALL.rds"))


spdf_world <- ne_countries(continent = c("africa", "europe", "north america", "south america", "asia", "oceania", "arctic"))
spTransform(spdf_world, projection(noExtinction_gmpd2_raster$raster[[1]]))
country_count <- nrow(spdf_world@data)
spdf_world@data$id <- 1:country_count
country_fort <- fortify(spdf_world, by = 'id')


spdf_world <- ne_countries(continent = c('africa', 'north america', 'south america', 'europe', 'asia', 'oceania'))
#spdf_world <- ne_countries()
spdf_world_fort <- fortify(spdf_world)

unique(noExtinction_gmpd2_raster$orders)

noExtinction_gmpd2_carnivora_raster <- noExtinction_gmpd2_raster[noExtinction_gmpd2_raster$orders == "CARNIVORA" & 
                                                               noExtinction_gmpd2_raster$data_type == "predicted_max",]
noExtinction_gmpd2_cetartiodactyla_raster <- noExtinction_gmpd2_raster[noExtinction_gmpd2_raster$orders == "CETARTIODACTYLA" & 
                                                                 noExtinction_gmpd2_raster$data_type == "predicted_max",]
noExtinction_gmpd2_primates_raster <- noExtinction_gmpd2_raster[noExtinction_gmpd2_raster$orders == "PRIMATES" & 
                                                                 noExtinction_gmpd2_raster$data_type == "predicted_max",]


Extinction_gmpd2_carnivora_raster <- Extinction_gmpd2_raster[Extinction_gmpd2_raster$orders == "CARNIVORA" & 
                                                           Extinction_gmpd2_raster$data_type == "predicted_max",]
Extinction_gmpd2_cetartiodactyla_raster <- Extinction_gmpd2_raster[Extinction_gmpd2_raster$orders == "CETARTIODACTYLA" & 
                                                                 Extinction_gmpd2_raster$data_type == "predicted_max",]
Extinction_gmpd2_primates_raster <- Extinction_gmpd2_raster[Extinction_gmpd2_raster$orders == "PRIMATES" & 
                                                          Extinction_gmpd2_raster$data_type == "predicted_max",]


plot(noExtinction_gmpd2_carnivora_raster$raster[[1]])
plot(Extinction_gmpd2_carnivora_raster$raster[[1]])
plot(noExtinction_gmpd2_cetartiodactyla_raster$raster[[1]])
plot(Extinction_gmpd2_cetartiodactyla_raster$raster[[1]])
plot(noExtinction_gmpd2_primates_raster$raster[[1]])
plot(Extinction_gmpd2_primates_raster$raster[[1]])


#Extinction_shaw_virus_biasLayers <- readRDS(paste0(base_dir, "intermediates/Extinction-Shaw_virus-Reso_sixth-bias_layers.rds"))$bias_shape[1]$viruses
#noExtinction_shaw_virus_biasLayers <- readRDS(paste0(base_dir, "intermediates/noExtinction-Shaw_virus-Reso_sixth-bias_layers.rds"))$bias_shape[1]$viruses

difference_gmpd2_carnivora_map <- Extinction_gmpd2_carnivora_raster$raster[[1]] - noExtinction_gmpd2_carnivora_raster$raster[[1]]
difference_gmpd2_carnivora_map <- (difference_gmpd2_carnivora_map / noExtinction_gmpd2_carnivora_raster$raster[[1]]) * 100
difference_gmpd2_cetartiodactyla_map <- Extinction_gmpd2_cetartiodactyla_raster$raster[[1]] - noExtinction_gmpd2_cetartiodactyla_raster$raster[[1]]
difference_gmpd2_cetartiodactyla_map <- (difference_gmpd2_cetartiodactyla_map / noExtinction_gmpd2_cetartiodactyla_raster$raster[[1]]) * 100
difference_shaw_primates_map <- Extinction_gmpd2_primates_raster$raster[[1]] - noExtinction_gmpd2_primates_raster$raster[[1]]
difference_shaw_primates_map <- (difference_shaw_primates_map / noExtinction_gmpd2_primates_raster$raster[[1]]) * 100

step_size <- 20
max_for_graphing <- plyr::round_any(max(difference_gmpd2_carnivora_map@data@max, 
                                        difference_gmpd2_cetartiodactyla_map@data@max,
                                        difference_shaw_primates_map@data@max,
                                
                                        abs(difference_gmpd2_carnivora_map@data@min), 
                                        abs(difference_gmpd2_cetartiodactyla_map@data@min),
                                        abs(difference_shaw_primates_map@data@min)
                        ), step_size) 


#combined_biasLayers_sldf <- Extinction_shaw_virus_biasLayers_sldf + noExtinction_shaw_virus_biasLayers_sldf
diff_bias <- cartomisc::gplot_data(difference_gmpd2_carnivora_map, maxpixels = 2332800)
#legend_title <- "Percentage difference in\nviral richness -\nCarnivora"
 
#abs_graphing <- round(max(abs(min(diff_bias$value, na.rm = T)), abs(max(diff_bias$value, na.rm = T))))

legend_title <- "Percentage difference\nin helminth richness -\nCarnivora"
diff_plot1 <- ggplot() +
  theme_classic() +
  geom_tile(data = diff_bias, aes(x = x, y = y, fill = value)) +
  scale_fill_gradientn(legend_title,  
                       colours = topo.colors(max(diff_bias$value, na.rm = TRUE), 
                                             alpha = 1, 
                                             rev = FALSE),
                       limits = c(min(diff_bias$value, na.rm = T), max(diff_bias$value, na.rm = T)),
                       na.value="lightgrey") +
  xlab("Longditude") +
  ylab("Latitude") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(), 
        axis.line = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  #scale_fill_gradient(low = "cyan", high = "orange", name = "UVI") +
  #geom_line(data = combined_biasLayers_sldf, aes(x = long, y = lat, group = group), colour = 'white') +
  geom_polygon(data = spdf_world_fort, aes(x = long, y = lat, group = group), colour = "black", fill = NA) 

#diff_plot1 <- diff_plot1 + 
#  scico::scale_fill_scico(legend_title, palette = "oleron", na.value = 'white', # previously 'vic'
#                          limits = c(-1, 1) * max_for_graphing) +
#  scico::scale_color_scico(legend_title, palette = "oleron", na.value = 'white', # previously 'vic'
#                           limits = c(-1, 1) * max_for_graphing) 

diff_plot1 <- diff_plot1 + 
  scale_fill_stepsn(legend_title,
                    colours=c("darkslateblue", "cornflowerblue", "white", "coral1", "darkred"),
                    breaks=seq(-max_for_graphing, max_for_graphing, by=step_size),
                    limits = c(-1, 1) * max_for_graphing, 
                    na.value = "white") +
  guides(fill = guide_colorsteps(
    barwidth=unit(30,'points'),
    barheight=unit(300,'points')
  )) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title = element_text(size=20))

#plot(diff_plot1)
ggsave(paste0(base_dir,"Figures/GMPD2_", map_type, "_helminth_withVSwithout_extinction-carnivora.png"), diff_plot1, dpi = 500, width = 15, height = 8)
#dev.off()

legend_title <- "Difference in helminth richness\nwith vs without extinctions"
diff_plot1 <- diff_plot1 + 
  scale_fill_stepsn(legend_title,
                    colours=c("darkslateblue", "cornflowerblue", "white", "coral1", "darkred"),
                    breaks=seq(-max_for_graphing, max_for_graphing, by=step_size),
                    limits = c(-1, 1) * max_for_graphing, 
                    na.value = "white") +
  guides(fill = guide_colorsteps(
    barwidth=unit(40,'points'),
    barheight=unit(400,'points')
  )) +
  theme(legend.text=element_text(size=23)) +
  theme(legend.title = element_text(size=23))

#diff_plot1 <- diff_plot1 +
#  scico::scale_fill_scico(legend_title, palette = "oleron", na.value = 'white', # previously 'vic'
#                          limits = c(-1, 1) * max_for_graphing) 
#plot(diff_plot1)
leg <- cowplot::get_legend(diff_plot1)
leg <- ggplotify::as.ggplot(leg)
diff_plot1 <- diff_plot1 + theme(legend.position = "none")



#combined_biasLayers_sldf <- Extinction_shaw_virus_biasLayers_sldf + noExtinction_shaw_virus_biasLayers_sldf
diff_bias <- cartomisc::gplot_data(difference_gmpd2_cetartiodactyla_map, maxpixels = 2332800)
#legend_title <- "Percentage difference in\nviral richness -\nCetartiodactyla"
legend_title <- "Percentage difference\nin helminth richness -\nCetartiodactyla"
diff_plot2 <- ggplot() +
  theme_classic() +
  geom_tile(data = diff_bias, aes(x = x, y = y, fill = value)) +
  scale_fill_gradientn(legend_title,  
                       colours = topo.colors(max(diff_bias$value, na.rm = TRUE), 
                                             alpha = 1, 
                                             rev = FALSE),
                       limits = c(min(diff_bias$value, na.rm = T), max(diff_bias$value, na.rm = T)),
                       na.value="lightgrey") +
  xlab("Longditude") +
  ylab("Latitude") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(), 
        axis.line = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  #scale_fill_gradient(low = "cyan", high = "orange", name = "UVI") +
  #geom_line(data = combined_biasLayers_sldf, aes(x = long, y = lat, group = group), colour = 'white') +
  geom_polygon(data = spdf_world_fort, aes(x = long, y = lat, group = group), colour = "black", fill = NA) 

#diff_plot2 <- diff_plot2 + 
#  scico::scale_fill_scico(legend_title, palette = "oleron", na.value = 'white', # previously 'vic'
#                          limits = c(-1, 1) * max_for_graphing) +
#  scico::scale_color_scico(legend_title, palette = "oleron", na.value = 'white', # previously 'vic'
#                           limits = c(-1, 1) * max_for_graphing) 

#plot(diff_plot2)

diff_plot2 <- diff_plot2 + 
  scale_fill_stepsn(legend_title,
                    colours=c("darkslateblue", "cornflowerblue", "white", "coral1", "darkred"),
                    breaks=seq(-max_for_graphing, max_for_graphing, by=step_size),
                    limits = c(-1, 1) * max_for_graphing, 
                    na.value = "white") +
  guides(fill = guide_colorsteps(
    barwidth=unit(30,'points'),
    barheight=unit(300,'points')
  )) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title = element_text(size=20))

ggsave(paste0(base_dir,"Figures/GMPD2_", map_type, "_helminth_withVSwithout_extinction-cetartiodactyla.png"), diff_plot2, dpi = 500, width = 15, height = 8)
#dev.off()

diff_plot2 <- diff_plot2 + theme(legend.position = "none")



#combined_biasLayers_sldf <- Extinction_shaw_virus_biasLayers_sldf + noExtinction_shaw_virus_biasLayers_sldf
diff_bias <- cartomisc::gplot_data(difference_shaw_primates_map, maxpixels = 2332800)
#legend_title <- "Percentage difference in\nviral richness -\nPrimates"
legend_title <- "Percentage difference\nin helminth richness -\nChiroptera"
diff_plot4 <- ggplot() +
  theme_classic() +
  geom_tile(data = diff_bias, aes(x = x, y = y, fill = value)) +
  scale_fill_gradientn(legend_title,  
                       colours = topo.colors(max(diff_bias$value, na.rm = TRUE), 
                                             alpha = 1, 
                                             rev = FALSE),
                       limits = c(min(diff_bias$value, na.rm = T), max(diff_bias$value, na.rm = T)),
                       na.value="lightgrey") +
  xlab("Longditude") +
  ylab("Latitude") +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(), 
        axis.line = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  #scale_fill_gradient(low = "cyan", high = "orange", name = "UVI") +
  #geom_line(data = combined_biasLayers_sldf, aes(x = long, y = lat, group = group), colour = 'white') +
  geom_polygon(data = spdf_world_fort, aes(x = long, y = lat, group = group), colour = "black", fill = NA) 

#diff_plot4 <- diff_plot4 + 
#  scico::scale_fill_scico(legend_title, palette = "oleron", na.value = 'white', # previously 'vic'
#                          limits = c(-1, 1) * max_for_graphing) +
#  scico::scale_color_scico(legend_title, palette = "oleron", na.value = 'white', # previously 'vic'
#                           limits = c(-1, 1) * max_for_graphing) 
#plot(diff_plot4)

diff_plot4 <- diff_plot4 + 
  scale_fill_stepsn(legend_title,
                    colours=c("darkslateblue", "cornflowerblue", "white", "coral1", "darkred"),
                    breaks=seq(-max_for_graphing, max_for_graphing, by=step_size),
                    limits = c(-1, 1) * max_for_graphing, 
                    na.value = "white") +
  guides(fill = guide_colorsteps(
    barwidth=unit(30,'points'),
    barheight=unit(300,'points')
  )) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title = element_text(size=20))

ggsave(paste0(base_dir,"Figures/Shaw_", map_type, "_helminth_withVSwithout_extinction-primates.png"), diff_plot4, dpi = 500, width = 15, height = 8)
#dev.off()

diff_plot4 <- diff_plot4 + theme(legend.position = "none")


x1 <- blankPlot
x2 <- blankPlot

# Combined plot
combined_plot <- cowplot::plot_grid(cowplot::plot_grid(diff_plot1, diff_plot2, diff_plot4, x1, x1,
                                                       nrow = 5,
                                                       align = "h",
                                                       rel_widths = c(10),
                                                       rel_heights = c(5, 5, 5, 5, 5),
                                                       labels = c("a) Carnivora", "b) Cetartiodactyla", "c) Primates", "", ""),
                                                       label_size = 15,
                                                       hjust = 0),
                                    cowplot::plot_grid(x2, leg, x2, x2, x2,
                                                       nrow = 5,
                                                       align = "h",
                                                       rel_widths = c(4),
                                                       rel_heights = c(5, 5, 5, 5, 5),
                                                       labels = c("", ""),
                                                       label_size = 15,
                                                       hjust = 0),
                                    nrow = 1,
                                    ncol = 2,
                                    align = "h",
                                    rel_widths=c(10,4),
                                    labels = c("",""),
                                    hjust = 0)
#plot(combined_plot)
#ggsave(file=paste0(base_dir, "Figures/", map_type, "_combined.pdf"), combined_plot, width = 15, height = 24)
#dev.off()

ggsave(file = paste0(base_dir, "Figures/", map_type, "_helminth_perOrder_combined.png"), combined_plot, dpi = 500, width = 15, height = 40)

