#### Load in libraries

library(sp)
library(sf)
library(tidyverse)
library(automap)
library(patchwork)
library(RColorBrewer)
library(gstat)
library(ggspatial)

#### Load in Data 

crs= 20137 # This is the EPSG code for our UTM zone "+proj=utm +zone=37"

filenames.da <- c("Processed Data/BaleMountainsNP_MeanAbundance_2010_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2011_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2012_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2013_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2015_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2016_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2017_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2019_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2020_WetSeason_res100m.csv",
                  "Processed Data/BaleMountainsNP_MeanAbundance_2021_WetSeason_res100m.csv")

da.wet.list <- lapply(filenames.da, read.csv)  # Create list of all dataframes for all years

da.wet.list <- lapply(da.wet.list, st_as_sf, coords = c("UTMX.1","UTMY.1"), crs = crs)

# Clean Data

hist(da.wet.list[[1]] %>% pull(Mountain.Nyala) %>% as.numeric()) # There is one Mountain Nyala record that is an order of magnitude larger than all other observations (500, presumably a typo, should be 50). We drop this.

da.wet.list[[1]] <- da.wet.list[[1]] %>% 
  filter(!(Mountain.Nyala >= "500"))

hist(da.wet.list[[1]] %>% pull(Mountain.Nyala) %>% as.numeric()) 

# Merge Livestock

livestock_names <- c("Cattle", "Sheep...goats", "Horse", "Donkey", "Mixed.Livestock") 
da.wet.list <- lapply(da.wet.list, function(df) {
  for (var in livestock_names) {
    if (!var %in% names(df)) {
      df[[var]] <- 0
    }
  }
  return(df)
})

da.wet.list <-lapply(da.wet.list, function(df) {
  df %>% mutate(Livestock=Cattle+Sheep...goats+Horse+Donkey+Mixed.Livestock)
}
)

hist(da.wet.list[[4]] %>% pull(Livestock) %>% as.numeric()) 


####  Make Grid ####


bbox <- st_bbox(da.wet.list[[3]])


# Minimum coordinates
min_x <- bbox["xmin"]
min_y <- bbox["ymin"]

# Lengths with padding
x_length <- (bbox["xmax"] - bbox["xmin"]) + 2000
y_length <- (bbox["ymax"] - bbox["ymin"]) + 2000

cellsize <- 500
ncol <- round(x_length / cellsize, 0)
nrow <- round(y_length / cellsize, 0)

# Create grid as sf object
x_seq <- seq(min_x, by = cellsize, length.out = ncol)
y_seq <- seq(min_y, by = cellsize, length.out = nrow)

grid_cells <- expand.grid(x = x_seq, y = y_seq)

grid_sf <- st_as_sf(grid_cells, coords = c("x", "y"), crs = crs)


#### Read in shapefile of Bale and crop grid to Gaysay Area ####

Bale_shapefile <- sf::st_read("Shapefiles/BaleMountainsSHP_WDPA_WDOECM_Feb2023_Public_2281_shp", "WDPA_WDOECM_Feb2023_Public_2281_shp-polygons")
Bale_shapefile %>% st_set_crs(st_crs(da.wet.list))

Bale_shapefile<-st_transform(Bale_shapefile, crs= crs)

xmin = min_x - 1000
ymin = min_y - 1000
xmax = min_x + x_length + 1000
ymax = min_y + y_length + 1000

bbox_sfc <- st_sfc(st_polygon(list(matrix(c(
  xmin, ymin,
  xmin, ymax,
  xmax, ymax,
  xmax, ymin,
  xmin, ymin
), ncol = 2, byrow = TRUE))))


Bale_shapefile<-st_crop(Bale_shapefile, bbox_sfc)

grid_sf <-st_crop(grid_sf, Bale_shapefile)



#### Read in shapefile of road network to add context to plots ####

Roads <- sf::st_read("Shapefiles/Roads", "Ethiopia_roads_network_clipped")

Roads<-st_crop(Roads, y = grid_sf)

# Interpolation for all large ungulate species 

years <- c(2010,2011,2012,2013,2015,2016,2017,2019,2020,2021)

### Livestock Plots


Livestock_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=Livestock ~ 1, newdata=grid_sf, idp=2.0)

livestock_min_value <- min(sapply(Livestock_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.05, na.rm = TRUE)
}))
livestock_max_value <- max(sapply(Livestock_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))


Livestock_plots <- lapply(1:length(Livestock_wet.idw.list), function(i) { ggplot() + 
    geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
    geom_sf(data = Livestock_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
    scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(livestock_min_value,livestock_max_value),
                          oob = scales::squish,
                          breaks = c(livestock_min_value, livestock_max_value),
                          labels = c("Lower", "Higher")) +  
    labs(color = "Relative Abundance") +
    xlab("Latitude") +
    ylab("Longitude") +
    ggtitle(paste("All Livestock")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })



#### Cattle Plots ####


Cattle_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=Cattle ~ 1, newdata=grid_sf, idp=2.0)
  
Cattle_min_value <- min(sapply(Cattle_wet.idw.list, function(sf) {
  min(sf$var1.pred, na.rm = TRUE)
  }))
Cattle_max_value <- max(sapply(Cattle_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))


Cattle_plots <- lapply(1:length(Cattle_wet.idw.list), function(i) { ggplot() + 
  geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
  geom_sf(data = Cattle_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
  scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(Cattle_min_value, Cattle_max_value),
                        oob = scales::squish,
                        breaks = c(Cattle_min_value, Cattle_max_value),
                        labels = c("Lower", "Higher")) +  
  labs(color = "Relative Abundance") +
  xlab("Latitude") +
  ylab("Longitude") +
  ggtitle(paste("Cattle")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })


# Sheep & Goats


Sheep.Goats_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=`Sheep...goats` ~ 1, newdata=grid_sf, idp=2.0)

Sheep_Goats_min_value <- min(sapply(Sheep.Goats_wet.idw.list, function(sf) {
  min(sf$var1.pred, na.rm = TRUE)
}))
Sheep_Goats_max_value <- max(sapply(Sheep.Goats_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))


Sheep_Goats_plots <- lapply(1:length(Sheep.Goats_wet.idw.list), function(i) { ggplot() + 
    geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
    geom_sf(data = Sheep.Goats_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
    scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(Sheep_Goats_min_value, Sheep_Goats_max_value),
                          oob = scales::squish,
                          breaks = c(Sheep_Goats_min_value, Sheep_Goats_max_value),
                          labels = c("Lower", "Higher")) +  
    labs(color = "Relative Abundance") +
    xlab("Latitude") +
    ylab("Longitude") +
    ggtitle(paste("Sheep & Goats")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })


# Horses


Horses_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=Horse ~ 1, newdata=grid_sf, idp=2.0)

Horse_min_value <- min(sapply(Horses_wet.idw.list, function(sf) {
  min(sf$var1.pred, na.rm = TRUE)
}))
Horse_max_value <- max(sapply(Horses_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))


Horse_plots <- lapply(1:length(Horses_wet.idw.list), function(i) { ggplot() + 
    geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
    geom_sf(data = Horses_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
    scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(Horse_min_value, Horse_max_value),
                          oob = scales::squish,
                          breaks = c(Horse_min_value, Horse_max_value),
                          labels = c("Lower", "Higher")) +  
    labs(color = "Relative Abundance") +
    xlab("Latitude") +
    ylab("Longitude") +
    ggtitle(paste("Horses")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })


# Mountain Nyala


MN_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=Mountain.Nyala ~ 1, newdata=grid_sf, idp=2.0)

MN_min_value <- min(sapply(MN_wet.idw.list, function(sf) {
  min(sf$var1.pred, na.rm = TRUE)
}))
MN_max_value <- max(sapply(MN_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))


MN_plots <- lapply(1:length(MN_wet.idw.list), function(i) { ggplot() + 
    geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
    geom_sf(data = MN_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
    scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(MN_min_value, MN_max_value),
                          oob = scales::squish,
                          breaks = c(MN_min_value, MN_max_value),
                          labels = c("Lower", "Higher")) +  
    labs(color = "Relative Abundance") +
    xlab("Latitude") +
    ylab("Longitude") +
    ggtitle(paste("Mountain Nyala")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })

# Bushbuck


BB_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=Menelik.s.Bushbuck ~ 1, newdata=grid_sf, idp=2.0)

BB_min_value <- min(sapply(BB_wet.idw.list, function(sf) {
  min(sf$var1.pred, na.rm = TRUE)
}))
BB_max_value <- max(sapply(BB_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))

BB_plots <- lapply(1:length(BB_wet.idw.list), function(i) { ggplot() + 
    geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
    geom_sf(data = BB_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
    scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(BB_min_value, BB_max_value),
                          oob = scales::squish,
                          breaks = c(BB_min_value, BB_max_value),
                          labels = c("Lower", "Higher")) +  
    labs(color = "Relative Abundance") +
    xlab("Latitude") +
    ylab("Longitude") +
    ggtitle(paste("Menelik's Bushbuck")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })

# Reedbuck


RB_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=Bohor.Reedbuck ~ 1, newdata=grid_sf, idp=2.0)

RB_min_value <- min(sapply(RB_wet.idw.list, function(sf) {
  min(sf$var1.pred, na.rm = TRUE)
}))
RB_max_value <- max(sapply(RB_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))

RB_plots <- lapply(1:length(RB_wet.idw.list), function(i) { ggplot() + 
    geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
    geom_sf(data = RB_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
    scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(RB_min_value, RB_max_value),
                          oob = scales::squish,
                          breaks = c(RB_min_value, RB_max_value),
                          labels = c("Lower", "Higher")) +  
    labs(color = "Relative Abundance") +
    xlab("Latitude") +
    ylab("Longitude") +
    ggtitle(paste("Bohor Reedbuck")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })

# Warthog 


WH_wet.idw.list <- lapply(da.wet.list, gstat::idw, formula=Warthog ~ 1, newdata=grid_sf, idp=2.0)

WH_min_value <- min(sapply(WH_wet.idw.list, function(sf) {
  min(sf$var1.pred, na.rm = TRUE)
}))
WH_max_value <- max(sapply(WH_wet.idw.list, function(sf) {
  quantile(sf$var1.pred,  probs = 0.95, na.rm = TRUE)
}))

WH_plots <- lapply(1:length(WH_wet.idw.list), function(i) { ggplot() + 
    geom_sf(data = Bale_shapefile, colour = "black", fill = NA, linewidth = 2) +
    geom_sf(data = WH_wet.idw.list[[i]], aes(color = var1.pred), size = 5) +  
    scale_color_gradientn(colours = brewer.pal(9, "YlOrRd"), limits = c(WH_min_value, WH_max_value),
                          oob = scales::squish,
                          breaks = c(WH_min_value, WH_max_value),
                          labels = c("Lower", "Higher")) +  
    labs(color = "Relative Abundance") +
    xlab("Latitude") +
    ylab("Longitude") +
    ggtitle(paste("Warthog")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90)) +
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
    #                                  pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
    #                                  style = ggspatial::north_arrow_nautical(
    #                                    fill = c("black", "white"), line_col = "black")) +
    geom_sf(data = Roads, color = "darkgrey", linetype = "dotted", linewidth = 1.5) +
    coord_sf(crs = 20137) })



##patchwork 

patchwork.plots<- lapply(1:length(Livestock_plots), function(i) {
(Livestock_plots[[i]]+MN_plots[[i]])/(WH_plots[[i]]+RB_plots[[i]]+BB_plots[[i]]) + plot_layout(guides = "collect") +
  plot_annotation(title = paste0("Heatmaps - ", years[[i]]),
                  theme = theme(
                    plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10)))) 
})


patchwork.plots[[10]]
