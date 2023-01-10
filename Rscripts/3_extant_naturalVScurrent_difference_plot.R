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
library(arsenal)

rm(list = ls())
if(!is.null(dev.list())) dev.off()
base_dir <- "C:/Users/tomos/OneDrive/PhD/ExtinctHosts_PathogenRichness/"
Phylacine_trait <- read.csv(paste0(base_dir, "PHYLACINE_v1.2.1/Data/Traits/Trait_data.csv"))
Phylacine_trait <- data.frame(Phylacine_trait)
Extant_species <- Phylacine_trait %>% filter(IUCN.Status.1.2 != 'EP' & IUCN.Status.1.2 != 'EX' & IUCN.Status.1.2 != 'EW')
Extant_species <- Extant_species %>% filter(Terrestrial == 1 | Aerial == 1)
Extant_species <- Extant_species %>% filter(Marine != 1)
Extant_species <- Extant_species %>% filter(Binomial.1.2 != 'Homo_sapiens')
no_extant_sp <- nrow(Extant_species)

extant_current_df <- readRDS(paste0(base_dir,'intermediates/extant_current_df.rda'))
extant_df <- readRDS(paste0(base_dir,'intermediates/extant_df.rda'))
extant_current_df$row <- as.numeric(extant_current_df$row)
extant_current_df$col <- as.numeric(extant_current_df$col)
extant_df$row <- as.numeric(extant_df$row)
extant_df$col <- as.numeric(extant_df$col)

extant_matrix <- array(rep(0, 360*142*no_extant_sp), dim=c(360,142,no_extant_sp))
extant_current_matrix <- array(rep(0, 360*142*no_extant_sp), dim=c(360,142,no_extant_sp))
extant_increase_matrix <- array(rep(0, 360*142*no_extant_sp), dim=c(360,142,no_extant_sp))
extant_decrease_matrix <- array(rep(0, 360*142*no_extant_sp), dim=c(360,142,no_extant_sp))
difference_matrix <- array(rep(0, 360*142*no_extant_sp), dim=c(360,142,no_extant_sp))
counter <- 0

for(i in 1:no_extant_sp){ #no_extant_sp
  #i <- 1
  #print(i)
  focal_host <- Extant_species$Binomial.1.2[i]
  #print(focal_host)
  focal_current <- NULL
  focal_current <- as.data.frame(extant_current_df[extant_current_df$species == focal_host,])
  
  for(j in 1:nrow(focal_current)){
    extant_current_matrix[focal_current$col[j], focal_current$row[j], i] <- 1
  }
  #image(extant_current_matrix[,,i])
  
  focal_natural <- as.data.frame(extant_df[extant_df$species == focal_host,])
  for(j in 1:nrow(focal_natural)){
    extant_matrix[focal_natural$col[j], focal_natural$row[j], i] <- 1
  }

  
  focal_diff1 <- NULL
  focal_diff2 <- NULL
  focal_diff3 <- NULL
  
  focal_diff1 <- suppressMessages(anti_join(focal_natural, focal_current))
  focal_diff2 <- suppressMessages(anti_join(focal_current, focal_natural))
  focal_diff3 <- rbind(focal_diff1, focal_diff2)

  for(j in 1:nrow(focal_diff1)){
    extant_decrease_matrix[focal_diff1$col[j], focal_diff1$row[j], i] <- 1
  }
  #image(extant_decrease_matrix[,,i])
  
  for(j in 1:nrow(focal_diff2)){
    extant_increase_matrix[focal_diff2$col[j], focal_diff2$row[j], i] <- 1
  }
  #image(extant_increase_matrix[,,i])
  
  for(j in 1:nrow(focal_diff3)){
    difference_matrix[focal_diff3$col[j], focal_diff3$row[j], i] <- 1
  }
  #image(difference_matrix[,,i])
  
  if(nrow(focal_natural) > nrow(focal_current)){
    #print(paste0("Index: ", i, " of ", no_extant_sp, "; host: ", focal_host))
    print(paste0(focal_host, ": Current range smaller than past range"))
    counter <- counter + 1
    #image(extant_matrix[,,i])
    #image(extant_decrease_matrix[,,i])
    #image(extant_current_matrix[,,i])
  }else{}
  
}
print(paste0("No. hosts with larger past ranges: ", counter))

saveRDS(extant_current_matrix, paste0(base_dir, "intermediates/extant_current_matrix"))
saveRDS(difference_matrix, paste0(base_dir, "intermediates/extant_diff_matrix"))
saveRDS(extant_increase_matrix, paste0(base_dir, "intermediates/extant_increase_matrix"))
saveRDS(extant_decrease_matrix, paste0(base_dir, "intermediates/extant_decrease_marix"))

cumulative_matrix_extant <- rowSums(extant_matrix, dims = 2)
image(cumulative_matrix_extant)

cumulative_matrix_extant_current <- rowSums(extant_current_matrix, dims = 2)
image(cumulative_matrix_extant_current)
cumulative_matrix_extant_decrease <- rowSums(extant_decrease_matrix, dims = 2)
image(cumulative_matrix_extant_decrease)
cumulative_matrix_extant_increase <- rowSums(extant_increase_matrix, dims = 2)
image(cumulative_matrix_extant_increase)
cumulative_matrix_extant_diff <- rowSums(difference_matrix, dims = 2)
image(cumulative_matrix_extant_diff)

saveRDS(cumulative_matrix_extant_current, paste0(base_dir, "intermediatescumulative_matrix_extant_current"))
saveRDS(cumulative_matrix_extant_diff, paste0(base_dir, "intermediates/cumulative_matrix_extant_diff"))
saveRDS(cumulative_matrix_extant_decrease, paste0(base_dir, "intermediates/cumulative_matrix_extant_decrease"))
saveRDS(cumulative_matrix_extant_increase, paste0(base_dir, "intermediates/cumulative_matrix_extant_increase"))

#image(cumulative_matrix_extant_current)
#image(cumulative_matrix_extant_diff)

cumulative_matrix_extant_current <- readRDS(paste0(base_dir, "intermediates/cumulative_matrix_extant_current"))
cumulative_matrix_extant_diff <- readRDS(paste0(base_dir, "intermediates/cumulative_matrix_extant_diff"))
cumulative_matrix_extant_decrease <- readRDS(paste0(base_dir, "intermediates/cumulative_matrix_extant_decrease"))
cumulative_matrix_extant_increase <- readRDS(paste0(base_dir, "intermediates/cumulative_matrix_extant_increase"))

melt_cum_extant_curr <- reshape2::melt(cumulative_matrix_extant_current)
melt_cum_extant_curr$value[melt_cum_extant_curr$value == 0] <- NA
#legend_title <- expression(atop(textstyle("No of mammal species"),
#                                atop(textstyle("lost since 130,000 ybp"),
#                                     )))

legend_title <- paste0("Number of extant mammals\n(human modified ranges):")

extant_curr_plot <- ggplot() +
  geom_raster(data = melt_cum_extant_curr, aes(Var1, 142-Var2, fill = value)) +
  scale_fill_gradientn(legend_title, #"No of extinct mammal\nspecies lost during\n130,000 ybp", 
                       colours = topo.colors(no_extant_sp, 
                                             alpha = 1, 
                                             rev = FALSE),
                       limits = c(min(melt_cum_extant_curr$value, na.rm = T), max(melt_cum_extant_curr$value, na.rm = T)),
                       #limits = c(min(melt_cum_extinct_sp, na.rm = T), max(melt_cum_extinct_sp$value, na.rm = T)),  #max(melt_cum_extant_sp$value, na.rm = T)), #to have the extinct and extant plots on the same scale
                       na.value="lightgrey") +
  xlab("Longditude (resolution 96.5 km)") +
  ylab("Latitude (resolution 96.5 km)") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 1.5), text = element_text(size = 15))  

plot(extant_curr_plot)
ggsave(paste0(base_dir,"Figures/Extant_current_mammals.png"), width = 14, height = 7)

leg1 <- cowplot::get_legend(extant_curr_plot)
leg1 <- ggplotify::as.ggplot(leg1)
extant_curr_plot <- extant_curr_plot + theme(legend.position = "none")




melt_diff_extant_decrease <- reshape2::melt(cumulative_matrix_extant_decrease)
melt_diff_extant_decrease$value[melt_diff_extant_decrease$value == 0] <- NA
#legend_title <- paste0("Per-pixel count of mammals with\nreduced ranges from human acivity")
legend_title <- paste0("Number of extant mammals\nwith an decreased distribution\n(per pixel):")
extant_decrease_plot <- ggplot() +
  geom_raster(data = melt_diff_extant_decrease, aes(Var1, 142-Var2, fill = value)) +
  scale_fill_gradientn(legend_title, #"No of extinct mammal\nspecies lost during\n130,000 ybp", 
                       colours = topo.colors(no_extant_sp, 
                                             alpha = 1, 
                                             rev = FALSE),
                       limits = c(min(melt_diff_extant_decrease$value, na.rm = T), max(melt_diff_extant_decrease$value, na.rm = T)),
                       #limits = c(min(melt_cum_extinct_sp, na.rm = T), max(melt_cum_extinct_sp$value, na.rm = T)),  #max(melt_cum_extant_sp$value, na.rm = T)), #to have the extinct and extant plots on the same scale
                       na.value="lightgrey") +
  xlab("Longditude (resolution 96.5 km)") +
  ylab("Latitude (resolution 96.5 km)") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 1.5), text = element_text(size = 15))  

plot(extant_decrease_plot)
ggsave(paste0(base_dir,"Figures/Extant_decrease_mammals.png"), width = 14, height = 7)

leg2 <- cowplot::get_legend(extant_decrease_plot)
leg2 <- ggplotify::as.ggplot(leg2)
extant_decrease_plot <- extant_decrease_plot + theme(legend.position = "none")






melt_diff_extant_increase <- reshape2::melt(cumulative_matrix_extant_increase)
melt_diff_extant_increase$value[melt_diff_extant_increase$value == 0] <- NA
#legend_title <- paste0("Per-pixel count of mammals with\nincreased ranges from human acivity")
legend_title <- paste0("Number of extant mammals\nwith an increased distribution\n(per pixel):")
extant_increase_plot <- ggplot() +
  geom_raster(data = melt_diff_extant_increase, aes(Var1, 142-Var2, fill = value)) +
  scale_fill_gradientn(legend_title, #"No of extinct mammal\nspecies lost during\n130,000 ybp", 
                       colours = topo.colors(no_extant_sp, 
                                             alpha = 1, 
                                             rev = FALSE),
                       limits = c(min(melt_diff_extant_increase$value, na.rm = T), max(melt_diff_extant_increase$value, na.rm = T)),
                       #limits = c(min(melt_cum_extinct_sp, na.rm = T), max(melt_cum_extinct_sp$value, na.rm = T)),  #max(melt_cum_extant_sp$value, na.rm = T)), #to have the extinct and extant plots on the same scale
                       na.value="lightgrey") +
  xlab("Longditude (resolution 96.5 km)") +
  ylab("Latitude (resolution 96.5 km)") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 1.5), text = element_text(size = 15))  

plot(extant_increase_plot)
ggsave(paste0(base_dir,"Figures/Extant_increase_mammals.png"), width = 14, height = 7)

leg3 <- cowplot::get_legend(extant_increase_plot)
leg3 <- ggplotify::as.ggplot(leg3)
extant_increase_plot <- extant_increase_plot + theme(legend.position = "none")



combined_legend <- cowplot::plot_grid(leg1, leg2, leg3,
                                      nrow = 3,
                                      align = "h",
                                      rel_widths = c(3),
                                      labels = c("", "", ""),
                                      label_size = 15,
                                      hjust = 0)

plot(combined_legend)
combined_plot <- cowplot::plot_grid(extant_curr_plot, extant_decrease_plot, extant_increase_plot,
                                    nrow = 3,
                                    align = "h",
                                    rel_widths = c(10),
                                    labels = c("a", "b", "c"),
                                    label_size = 15,
                                    hjust = 0)

combined_plot <- cowplot::plot_grid(combined_plot, combined_legend,
                                    ncol = 2,
                                    align = "h",
                                    rel_widths = c(10, 3),
                                    labels = c("", ""),
                                    label_size = 15,
                                    hjust = 0)


ggsave(file=paste0(base_dir, "Figures/extant_curr_diff_combined.png"), combined_plot, dpi = 500, width = 14, height = 21)


