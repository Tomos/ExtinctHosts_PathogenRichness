library(tidyverse)
library(stringi)
library(sp)
library(rgdal)
library(raster)
library(classInt)
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
library(lattice)
library(rasterVis)

#n_cores = detectCores()
#n_cores_use = round(nrow(models) / ceiling(nrow(models) / (n_cores - 1)))
#options(mc.cores = n_cores_use)

rm(list = ls())
if(!is.null(dev.list())) dev.off()
base_dir <- "C:/Users/tomos/OneDrive/PhD/Host-parasite-paper/"
#P <- rprojroot::find_rstudio_root_file
#source(paste0(base_dir, "Olival-functions.R"))
# Load the Mammals' shapefile TERRESTRIAL_MAMMALS, downloaded from: http://www.iucnredlist.org/technical-documents/spatial-data#mammals
# This study uses version 2015-2, which we provide via Amazon S3 storage.
terr = shapefile(paste0("C:/Users/tomos/OneDrive/PhD/Literature/Olival_et_al_2017_dataset/Mammals_Terrestrial/Mammals_Terrestrial/Mammals_Terrestrial.shp"), verbose = T)
terr@data$BINOMIAL = stri_replace_all_fixed(terr@data$BINOMIAL, " ", "_")

# select only extant species
terr = subset(terr, PRESENCE == 1)

# Read taxonomic information.
taxa = read_csv(paste0("C:/Users/tomos/OneDrive/PhD/Literature/Olival_et_al_2017_dataset/ecohealthalliance-HP3-928327a/data/IUCN_taxonomy_23JUN2016.csv")) %>%
  dplyr::select(-c(7:23), spp_ID = `Species ID`)

# Join taxonomic information to spatial data
terr@data = left_join(terr@data, taxa, by = c("ID_NO" = 'spp_ID'))

hp3_orders = c('CARNIVORA', 'CETARTIODACTYLA', 'CHIROPTERA', 'PRIMATES', 'RODENTIA')
# Read HP3 results
# Get HP3 all virus predictions

tag <- "noExtinction"
repo_path <- "Shaw_bacteria"
repo <- "Shaw"
reso <- 1/6
reso_txt <- "Reso_sixth"
all_bacteria_gam <- readRDS(paste0(base_dir, "intermediates/noExtinction_shaw_bact_models_tw.rds"))$model[[1]]
hp3_hosts <- read_rds(paste0(base_dir, "intermediates/Shaw_bacteria_dataset.rds"))
hp3_all = as_tibble(left_join(all_bacteria_gam$model, hp3_hosts))
hp3_all_pred = within(hp3_all, {hDiseaseZACitesLn = max(hDiseaseZACitesLn)})
hp3_all = mutate(hp3_all,
                 observed = TotVirusPerHost,
                 predicted = as.vector(unname(predict(all_bacteria_gam, hp3_all, type="response"))),
                 predicted_max = as.vector(unname(predict(all_bacteria_gam, hp3_all_pred, type="response"))),
                 resid = predicted - observed,
                 missing = predicted_max - observed,
                 number = 1,
                 hHostNameFinal = as.character(hHostNameFinal))

# just zoonotic viruses
##zoo_viruses_gam <- readRDS(paste0(base_dir, "ecohealthalliance-HP3-928327a/intermediates/all_zoonoses_models.rds"))$model[[1]]
##hp3_zoo = left_join(zoo_viruses_gam$model, hp3_hosts) %>%
##  left_join(dplyr::select(hp3_all, hHostNameFinal, predicted_max)) %>%
##  dplyr::rename(vir_prediction_max = predicted_max) %>%
##  as_tibble()
##hp3_zoo_pred =  within(hp3_zoo, {hDiseaseZACites = max(hDiseaseZACites)
##LnTotNumVirus = log(vir_prediction_max)})

##hp3_zoo = mutate(hp3_zoo,
##                 observed = NSharedWithHoSa,
##                 predicted = as.vector(unname(predict(zoo_viruses_gam, hp3_zoo, type="response"))),
##                 predicted_max = as.vector(unname(predict(zoo_viruses_gam, hp3_zoo_pred, type="response"))),
##                 resid = predicted - observed,
##                 missing = predicted_max - observed,
##                 number = 1,
##                 hHostNameFinal = as.character(hHostNameFinal))


hp3_hosts = read_rds(paste0(base_dir, "intermediates/Shaw_bacteria_dataset.rds")) %>%
  filter(hWildDomFAO == 'wild', hMarOTerr == 'Terrestrial') %>%
  dplyr::select(hHostNameFinal) %>%
  mutate(hp3 = 1)


# Create spatial polygons for all_viruses & all_zoonoses
data_hp3_all <- terr
##data_hp3_zoo <- terr
data_host <- terr
#z <- data_hp3_all@data
# Join spatial polygons and hp3 data frames
data_hp3_all@data = full_join(data_hp3_all@data, hp3_all, by = c("BINOMIAL" = 'hHostNameFinal')) %>%
  mutate(order_group = if_else(Order %in% hp3_orders, Order, "OTHER"))

z <- data_hp3_all@data

##data_hp3_zoo@data = full_join(data_hp3_zoo@data, hp3_zoo, by = c("BINOMIAL" = 'hHostNameFinal')) %>%
##  mutate(order_group = if_else(Order %in% hp3_orders, Order, "OTHER"))

data_host@data <- left_join(data_host@data, hp3_hosts, by = c("BINOMIAL" = 'hHostNameFinal')) %>%
  mutate(observed = hp3) %>%
  mutate(predicted = 1) %>%
  mutate(missing = if_else(is.na(hp3) & predicted == 1, 1, 0)) %>%
  mutate(order_group = if_else(Order %in% hp3_orders, Order, "OTHER"))

zz <- data_host@data

# Remove NAs (i.e., species with no data available)
data_hp3_all = data_hp3_all[!is.na(data_hp3_all$observed),]


##data_hp3_zoo = data_hp3_zoo[!is.na(data_hp3_zoo$observed),]

# Auxilliary functions

# Create color palette. Requires rasterVis and RColorBrewer
myTheme = rasterTheme(region = rev(brewer.pal(11, 'RdYlGn')))
myTheme$fontsize$text <- 5.6
myTheme$axis.line$lwd <- 0.2

myTheme2 = RdBuTheme()
myTheme2$fontsize$text <- 5.6
myTheme2$axis.line$lwd <- 0.2

myTheme3 <- rasterTheme(region = viridis::viridis(11))
myTheme3$fontsize$text <- 5.6
myTheme3$axis.line$lwd <- 0.2

data(wrld_simpl)
wrld_simpl = subset(wrld_simpl, NAME != 'Antarctica')
world_layer <-latticeExtra::layer(sp::sp.polygons(wrld_simpl, lwd=1, col='gray50'))

make_map <- function(model, orders, data_type, raster_res) {
  if(model == 'hosts') {
    sp_data_frame = data_host
  } else if (model == 'zoonoses') {
    sp_data_frame = data_hp3_zoo
  } else if (model == 'bacteria') {
    sp_data_frame = data_hp3_all
  }
  template <- raster(resolution = raster_res)
  subsetted_shapes <- subset(sp_data_frame, (order_group %in% orders))
  my_raster = raster()
  my_raster = rasterize(subsetted_shapes, template, data_type, fun = 'sum', silent = TRUE)
  #make_png(my_raster, theme, filename, png_res)
  return(my_raster)
}


make_cutoff <- function(resids, predictions, sig_level = 0.05, nreps = 1000, vals = seq(from=1, to=100, by=1)) {
  N = length(resids)
  stopifnot(N == length(predictions))
  vals <- unique(c(vals[vals <= N], N))
  q <- matrix(NA, nrow = nreps, ncol=length(vals))
  for (i in seq_len(nreps)) {
    for (j in seq_along(vals)) {
      samp = sample.int(N, vals[j])
      q[i, j] <- sum(resids[samp])/sum(predictions[samp])
    }
  }
  cutoffs <- apply(q, 2, function(dis) {
    quantile(dis, c(sig_level/2, 1-sig_level/2))
  })
  return(list(lower=approxfun(vals, cutoffs[1,], rule=1),
              upper=approxfun(vals, cutoffs[2,], rule=1)))
}

#rasters = crossing(orders=c(hp3_orders,"OTHER"),
#                   model=c('hosts'),
#                   ###model=c('hosts', 'zoonoses', 'viruses'),
#                   data_type = c('observed', 'predicted', 'predicted_max', 'missing', 'resid', 'number')
#) %>%
#  filter(!(model == 'hosts' & data_type %in% c('resid', 'number', 'predicted_max')))

rasters = crossing(orders=c(hp3_orders,"OTHER"),
                   #model=c('hosts', 'zoonoses', 'bacteria'), #####
                   model=c('bacteria'), #####
                   data_type = c('observed', 'predicted', 'predicted_max', 'missing', 'resid', 'number')
                   #data_type = c('predicted_max')
) %>%
  filter(!(model == 'hosts' & data_type %in% c('resid', 'number', 'predicted_max')))


#rs <- make_map(
#   model = rasters$model[3],
#   orders = rasters$orders[3],
#   data_type = rasters$data_type[3],
#   raster_res = 1
#)



tic()
rasters$raster <- mcmapply(make_map,
                           model = rasters$model,
                           orders = rasters$orders,
                           data_type = rasters$data_type,
                           raster_res = reso,
                           mc.preschedule = FALSE,
                           mc.cores = 1, 
                           mc.silent= TRUE
                           )
toc()

saveRDS(rasters, paste0(base_dir, "intermediates/", tag, "-", reso_txt, "-", repo_path, "_rasters.rds"))
unique(data_hp3_all@data$Order)
unique(rasters$orders)

############################
rasters <- readRDS(paste0(base_dir, "intermediates/", tag, "-", reso_txt, "-", repo_path, "_rasters.rds"))

rasters2 <- rasters %>%
  group_by(model, data_type) %>%
  do({
    new_row = .[1,]
    new_row$orders <- "ALL"
    sum_raster = do.call("sum", c(unname(.$raster, NULL), na.rm=TRUE))

    sum_raster = calc(sum_raster, fun=function(x) ifelse(x == 0, NA, x))
    print(sum_raster) #####
    plot(sum_raster) #####
    new_row$raster <- list(sum_raster)
    #out = bind_rows(., new_row)
    saveRDS(new_row, paste0(base_dir, "maps/", tag, "-", reso_txt, "-map_intermediate_", repo, "-", new_row$model, "-", new_row$data_type))
    #return(out)
  })

bacteria_missing <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-map_intermediate_", repo, "-", "bacteria-missing"))
bacteria_number <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-map_intermediate_", repo, "-", "bacteria-number"))
bacteria_observed <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-map_intermediate_", repo, "-", "bacteria-observed"))
bacteria_predicted <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-map_intermediate_", repo, "-", "bacteria-predicted"))
bacteria_predicted_max <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-map_intermediate_", repo, "-", "bacteria-predicted_max"))
bacteria_resid <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-map_intermediate_", repo, "-", "bacteria-resid"))

rasters2 <- bind_rows(list(
    rasters, 
    bacteria_missing,
    bacteria_number,
    bacteria_observed,
    bacteria_predicted,
    bacteria_predicted_max,
    bacteria_resid
  ))

rasters3 <- rasters2 %>%
  group_by(orders, model) %>%
  do({
    if(all(.$model == "hosts")) {
      return(.)
    } else {
      new_row = .[1,]
      new_row$data_type <- "bias"
      bias_raster = .$raster[[which(.$data_type == "resid")]]/.$raster[[which(.$data_type == "predicted")]]
      plot(bias_raster)
      new_row$raster <- list(bias_raster)
      #out = bind_rows(., new_row)
      print(new_row)
      saveRDS(new_row, paste0(base_dir, "maps/", tag, "-", reso_txt, "-", "map_intermediate_", repo, "-", new_row$model, "-", new_row$orders, "-BIAS"))
    }
  }) 

carn_bacteria_bias <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-", "map_intermediate_", repo, "-bacteria-CARNIVORA-BIAS"))
cert_bacteria_bias <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-", "map_intermediate_", repo, "-bacteria-CETARTIODACTYLA-BIAS"))
chir_bacteria_bias <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-", "map_intermediate_", repo, "-bacteria-CHIROPTERA-BIAS"))
prim_bacteria_bias <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-", "map_intermediate_", repo, "-bacteria-PRIMATES-BIAS"))
rode_bacteria_bias <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-", "map_intermediate_", repo, "-bacteria-RODENTIA-BIAS"))
all_bacteria_bias <- readRDS(paste0(base_dir, "maps/", tag, "-", reso_txt, "-", "map_intermediate_", repo, "-bacteria-ALL-BIAS"))

rasters3 <- bind_rows(list(
  rasters2, 
  carn_bacteria_bias,
  cert_bacteria_bias,
  chir_bacteria_bias,
  prim_bacteria_bias,
  rode_bacteria_bias,
  all_bacteria_bias
))

rasters4 <- rasters3 %>%
  group_by() %>%
  arrange(orders, model, data_type) %>%
  filter(orders != "OTHER")

saveRDS(rasters4, paste0(base_dir, "intermediates/", tag, "-Shaw_bacteria_rasters_ALL.rds"))

bias_layers = rasters4 %>%
  filter(model != "hosts") %>%
  filter(data_type %in% c("bias", "number")) %>%
  spread("data_type", "raster")


angle_lines <- SpatialLines(lapply(seq(from=0, to=720, by=1), function(z) Lines(list(Line(matrix(c(-720 + z, z, z, z-720), ncol=2))), ID = z)),  proj4string = CRS(proj4string(bias_layers$bias[[1]])))

bias_cutshape <- function(model, orders, bias_raster, number_raster, cutoff = 0.05) {
  if(model == "bacteria") {
    mod_data <- hp3_all
  } else if (model=="zoonoses") {
    mod_data <- hp3_zoo
  }
  if(orders != "ALL") {
    mod_data = mod_data[mod_data$hOrder == orders,]
  }

  bias_funs <- make_cutoff(mod_data$resid, mod_data$predicted, nreps=100000, sig_level = 0.05, vals = seq(from=1, to=101, by=4))
  bias_matrix <- as.matrix(bias_raster)
  num_matrix <- as.matrix(number_raster)
  bias_mat2 <- ifelse((bias_matrix > bias_funs$upper(num_matrix)) | bias_matrix < bias_funs$lower(num_matrix), 1, NA)

  bias_shapes <- fasteraster::raster2vector(bias_mat2) %>%
  {Filter(function(x) nrow(x) > 3, .) } %>%
    lapply(function(x) cbind(360*x[,2]/ncol(bias_matrix) - 180, -180*x[,1]/nrow(bias_matrix) + 90))
  z <- seq_along(bias_shapes)
  bias_p <- SpatialPolygons(lapply(z, function(z) Polygons(list(Polygon(bias_shapes[[z]])), z)), proj4string = CRS(proj4string(bias_raster)))
  #saveRDS(bias_p, paste0(base_dir, "test.rds"))
  if(length(bias_p) > 0){
    bias_out <- rgeos::gIntersection(angle_lines, gBuffer(bias_p, width = .0001))
  }
}

bias_layers$bias_shape <- mcmapply(bias_cutshape, model=bias_layers$model, orders=bias_layers$orders, bias_raster = bias_layers$bias, number_raster = bias_layers$number,
                                   mc.cores = 1)

saveRDS(bias_layers, file = paste0(base_dir, "intermediates/", tag, "-", repo_path, "-", reso_txt, "-bias_layers.rds"))

##bias_layers <- readRDS(paste0(base_dir, "intermediates/", repo_path, "-", reso_txt, "bias_layers.rds"))

##make_png <- function(my_raster, orders, model, data_type, png_res) {
##  if(model == "hosts" & data_type=="missing") {
##    TheTheme <- BuRdTheme()
##  } else if ((data_type %in% c("predicted_max", "observed", "missing")) | model=="hosts") {
##    TheTheme <- rasterTheme(region = rev(brewer.pal(11, 'RdYlGn')))
##  } else {
##    TheTheme <- rasterTheme(region = viridis::viridis(11))
##  }
##  TheTheme$fontsize$text <- 30
##  TheTheme$axis.line$lwd <- 2
  

##  if(data_type %in% c("predicted_max", "missing") &
##     model != "hosts") {
##    bias_pt_layer = bias_layers$bias_shape[[which(bias_layers$orders == orders & bias_layers$model == model)]]
##  } else {
##    bias_pt_layer = angle_lines[0,]
##  }

##  filename = paste0(base_dir, "intermediates2/", repo, "_maps", orders, "_", model, "_", data_type, ".png")
##  png(filename, width = ncol(my_raster) + 1, height = nrow(my_raster) + 1, res = png_res, type = "cairo-png", antialias = "subpixel", family = "Arial")
##  print(lattice::levelplot(my_raster, 
##                           layers = 1,
##                           par.settings = TheTheme,
##                           margin= FALSE,
##                           ylab = '',
##                           xlab = '',
##                           scales = list(draw=FALSE),
##                           colorkey=list(width = 1,
##                                         axis.text = list(fontfamily="Helvetica", fontsize=5.6)),
##                           maxpixels = ncell(my_raster),
##                           xlim = c(-180, 180),
##                           ylim = c(-58, 90)
##                           ) +
          #layer(sp.points(bias_pt_layer, pch=20, cex=0.8, col="grey20"), data=list(bias_pt_layer=bias_pt_layer)) +
##          layer(sp.lines(bias_pt_layer, lwd=3, col="grey20"), data=list(bias_pt_layer=bias_pt_layer)) +
##          world_layer)
##  dev.off()
  #return(bias_pt_layer)

##}

##rasters5 = filter(rasters4, data_type %in% c('observed', 'predicted', 'predicted_max', 'missing')) %>%
##  filter(!(model=="hosts" & data_type=="predicted_max")) %>%
##  filter(!(model!="hosts" & data_type=="predicted"))

##mcmapply(make_png, my_raster=rasters5$raster, orders=rasters5$orders, model=rasters5$model, data_type=rasters5$data_type, png_res=150, mc.cores = 1)

