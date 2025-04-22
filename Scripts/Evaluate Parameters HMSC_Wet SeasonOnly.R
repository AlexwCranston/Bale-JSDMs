# ---- Load Libraries

library(Hmsc)
library(colorspace)
library(vioplot)
library(corrplot)
library(reshape2)
library(patchwork)
library(ggspatial)
library(ggplot2)

# ---- Set Directory

localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

# Load Models and Set Names

load("models/models.CA_thin_100_samples_250_chains_4.Rdata")
models_fitted.CA <- models_fitted

filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
load(filename_unfitted)
models.CA.list <- unfitted.models.list[[1]]
names(models_fitted.CA) <- names(models.CA.list)


load("models/models.presab_thin_100_samples_250_chains_4.Rdata")
models_fitted.presab <- models_fitted

filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
load(filename_unfitted)
models.presab.list <- unfitted.models.list[[2]]
names(models_fitted.presab) <- names(models.presab.list)

# Input Years


Year <- c(2010,2011,2012,2013,2015,2016,2017,2019,2020,2021)

# Evaluate Model Fit

preds.CA.list <-  lapply(models_fitted.CA, function(model) {
  model<-computePredictedValues(model)
  return(model)
})


results.CA.list <- lapply(seq_along(models_fitted.CA), function(i) {
  results.CA <- evaluateModelFit(hM=models_fitted.CA[[i]], predY=preds.CA.list[[i]])
  results.CA<-as.data.frame(results.CA)
  rownames(results.CA)<-colnames(models_fitted.CA[[i]]$Y)
  return(results.CA)
})


preds.presab.list <- lapply(models_fitted.presab, function(model) {
  model<-computePredictedValues(model)
  return(model)
})


results.presab.list <- lapply(seq_along(models_fitted.presab), function(i) {
  results.presab <- evaluateModelFit(hM=models_fitted.presab[[i]], predY=preds.presab.list[[i]])
  results.presab<-as.data.frame(results.presab)
  rownames(results.presab)<-colnames(models_fitted.presab[[i]]$Y)
  return(results.presab)
})


## Correlation Matrix

## First CA models

OmegaCor.list <- lapply(models_fitted.CA, function(model) {
  computeAssociations(model)})

supportLevel = 0.95

toPlot.list <-  lapply(seq_along(OmegaCor.list), function(i) {
  result <- OmegaCor.list[[i]][[1]]$support > supportLevel | OmegaCor.list[[i]][[1]]$support < (1 - supportLevel)
  result <- result * OmegaCor.list[[i]][[1]]$mean
  return(result)
})

cor_matrix_melted.list <- lapply(toPlot.list, function(matrix) {
 melted_matrix <- melt(matrix)
 return(melted_matrix)
}) 

corrplots.list.CA <- lapply(seq_along(cor_matrix_melted.list), function(i) {
  plot <- ggplot(cor_matrix_melted.list[[i]], aes(Var1, Var2, fill = value)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "gray", limits = c(-1,1)) + 
    scale_x_discrete(labels = c(
      "Mountain.Nyala" = "Mountain nyala",
      "Bohor.Reedbuck" = "Bohor reedbuck",
      "Menelik.s.Bushbuck" = "Menelik's bushbuck"
    )) +
    scale_y_discrete(labels = c(
      "Mountain.Nyala" = "Mountain nyala",
      "Bohor.Reedbuck" = "Bohor reedbuck",
      "Menelik.s.Bushbuck" = "Menelik's bushbuck"
    )) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, size = 7, face = "bold"),
          axis.text.y = element_text(size = 7, face = "bold") ) +
    labs(title = paste(Year[[i]]), x = "", y = "") + 
    coord_fixed() 
  return(plot)
})


(corrplots.list.CA[[1]] + corrplots.list.CA[[2]]+corrplots.list.CA[[3]])/
  (corrplots.list.CA[[4]] + corrplots.list.CA[[5]] + corrplots.list.CA[[6]])/
  (corrplots.list.CA[[7]] + corrplots.list.CA[[9]]+corrplots.list.CA[[10]]) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Correlation Matrices - Conditional Abundance Models",
                  theme = theme(
                    plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10))))

# Then Presab models 

OmegaCor.list <- lapply(models_fitted.presab, function(model) {
  computeAssociations(model)})

supportLevel = 0.95

toPlot.list <-  lapply(seq_along(OmegaCor.list), function(i) {
  result <- OmegaCor.list[[i]][[1]]$support > supportLevel | OmegaCor.list[[i]][[1]]$support < (1 - supportLevel)
  result <- result * OmegaCor.list[[i]][[1]]$mean
  return(result)
})

cor_matrix_melted.list <- lapply(toPlot.list, function(matrix) {
  melted_matrix <- melt(matrix)
  return(melted_matrix)
}) 

corrplots.list.presab <- lapply(seq_along(cor_matrix_melted.list), function(i) {
  plot <- ggplot(cor_matrix_melted.list[[i]], aes(Var1, Var2, fill = value)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "gray", limits = c(-1,1)) + 
    scale_x_discrete(labels = c(
      "Mountain.Nyala" = "Mountain nyala",
      "Bohor.Reedbuck" = "Bohor reedbuck",
      "Menelik.s.Bushbuck" = "Menelik's bushbuck"
    )) +
    scale_y_discrete(labels = c(
      "Mountain.Nyala" = "Mountain nyala",
      "Bohor.Reedbuck" = "Bohor reedbuck",
      "Menelik.s.Bushbuck" = "Menelik's bushbuck"
    )) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, size = 7, face = "bold"),
          axis.text.y = element_text(size = 7, face = "bold") ) +
    labs(title = paste(Year[[i]]), x = "", y = "") + 
    coord_fixed() 
    return(plot)
})

(corrplots.list.presab[[1]] + corrplots.list.presab[[2]]+corrplots.list.presab[[3]])/
  (corrplots.list.presab[[4]] + corrplots.list.presab[[5]] + corrplots.list.presab[[6]])/
  (corrplots.list.presab[[7]] + corrplots.list.presab[[9]]+corrplots.list.presab[[10]]) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Correlation Matrices - Presence/Absence Models",
                  theme = theme(
                    plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10))))
  



# Gamma parameters

## First Conditional Abundance

postGamma.CA.list <- lapply(models_fitted.CA, function(model) {
  getPostEstimate(model, parName = "Gamma")
  })


toPlot_Gamma.list <-  lapply(seq_along(postGamma.CA.list), function(i) {
  result <- postGamma.CA.list[[i]]$support > supportLevel | postGamma.CA.list[[i]]$support < (1 - supportLevel)
  result <- result * postGamma.CA.list[[i]]$mean
  return(result)
})


gamma_matrix_melted.list <- lapply(toPlot_Gamma.list, function(matrix) {
  melted_matrix <- melt(matrix)
  return(melted_matrix)
}) 

# Define rename variables by their original names
var_labels <- c(
  "1" = "Grassland",
  "2" = "Shrubs/Heath",
  "3" = "Wetland",
  "4" = "Woodland",
  "5" = "Elevation",
  "6" = "Distance from Road")

y_labels <- c("1" = "Wild",
              "2" = "Domestic")


gamma_matrix_melted.list <- lapply(gamma_matrix_melted.list, function(df) {
      df$Var1 <- var_labels[as.character(df$Var1)]
      df$Var2 <- y_labels[as.character(df$Var2)]
      return(df) })

# Define the desired order for plotting 
var1_order <- c("Grassland", "Shrubs/Heath", "Wetland", "Woodland", "Elevation", "Distance from Road")
var2_order <- c("Wild", "Domestic")

gamma_matrix_melted.list <- lapply(gamma_matrix_melted.list, function(df) {
  df$Var1 <- factor(df$Var1, levels = var1_order)
  df$Var2 <- factor(df$Var2, levels = var2_order)
  return(df)
})

gamma_plots.list.CA <- lapply(seq_along(gamma_matrix_melted.list), function(i) {
  plot <- ggplot(gamma_matrix_melted.list[[i]], aes(Var1, Var2, fill = value)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "gray", limits = c(-1,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust=0.5)) +
    labs(title = paste(Year[[i]]), x = "", y = "") + 
    coord_fixed() 
  return(plot)
})

(gamma_plots.list.CA[[6]]/gamma_plots.list.CA[[7]]/gamma_plots.list.CA[[9]]) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Gamma Plot - Conditional Abundance Models",
                  theme = theme(
                    plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10))))

## Now presab


postGamma.presab.list <- lapply(models_fitted.presab, function(model) {
  getPostEstimate(model, parName = "Gamma")
})


toPlot_Gamma.list <-  lapply(seq_along(postGamma.presab.list), function(i) {
  result <- postGamma.presab.list[[i]]$support > supportLevel | postGamma.presab.list[[i]]$support < (1 - supportLevel)
  result <- result * postGamma.presab.list[[i]]$mean
  return(result)
})


gamma_matrix_melted.list <- lapply(toPlot_Gamma.list, function(matrix) {
  melted_matrix <- melt(matrix)
  return(melted_matrix)
}) 


gamma_matrix_melted.list <- lapply(gamma_matrix_melted.list, function(df) {
  df$Var1 <- var_labels[as.character(df$Var1)]
  df$Var2 <- y_labels[as.character(df$Var2)]
  return(df) })

gamma_matrix_melted.list <- lapply(gamma_matrix_melted.list, function(df) {
  df$Var1 <- factor(df$Var1, levels = var1_order)
  df$Var2 <- factor(df$Var2, levels = var2_order)
  return(df)
})

gamma_plots.list.presab <- lapply(seq_along(gamma_matrix_melted.list), function(i) {
  plot <- ggplot(gamma_matrix_melted.list[[i]], aes(Var1, Var2, fill = value)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "gray", limits = c(-1,1)) +
    theme_minimal() +
    labs(title = paste(Year[[i]]), x = "", y = "") + 
    coord_fixed() 
  return(plot)
})

(gamma_plots.list.presab[[1]] + gamma_plots.list.presab[[2]]+gamma_plots.list.presab[[3]])/
  (gamma_plots.list.presab[[4]] + gamma_plots.list.presab[[5]] + gamma_plots.list.presab[[6]])/
  (gamma_plots.list.presab[[7]] + gamma_plots.list.presab[[9]]+gamma_plots.list.presab[[10]]) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Gamma Plot - Conditional Abundance Models",
                  theme = theme(
                    plot.title = element_text(hjust = 0.5, size = 11, face = "bold", margin = margin(b = 10))))


## Eta Parameters ####

postEta.presab <- lapply(models_fitted.presab, function(model) {
  getPostEstimate(model, parName ="Eta")
})

xy.list<- lapply(models_fitted.presab, function(model) {
  xy <- model$rL$ID$s
  return(xy)
})


eta <- lapply(postEta.presab, function(model) {
  eta <- as.data.frame(model[[1]])
  return(eta)
})

eta <- lapply(seq_along(eta), function(i) { 
  eta <- cbind(xy.list[[i]],eta[[i]])
  return(eta)
})


Bale_shapefile <- sf::st_read("Shapefiles/BaleMountainsSHP_WDPA_WDOECM_Feb2023_Public_2281_shp", "WDPA_WDOECM_Feb2023_Public_2281_shp-polygons")
Bale_shapefile %>% st_set_crs("+proj=utm +zone=37")
Bale_shapefile<-st_transform(Bale_shapefile, crs = "+proj=utm +zone=37")

min_x = 570000 #minimun x coordinate
min_y = 780000 #minimun y coordinate
x_length = 605694-min_x+2000 #easting amplitude
y_length = 790661.6-min_y+2000

Bale_shapefile<-st_crop(Bale_shapefile,y=c(xmin=(min_x-1000), ymin=(min_y-1000), ymax=(min_y+y_length+1000),xmax=(min_x+x_length+1000))) # crop to size of kriging output
Bale_shapefile %>% st_set_crs("+proj=longlat +datum=WGS84")
Bale_shapefile<-st_transform(Bale_shapefile, crs = "+proj=longlat +datum=WGS84")
Bale_shapefile<- as_Spatial(Bale_shapefile)


Roads <- sf::st_read("Shapefiles/Roads", "Ethiopia_roads_network_clipped")
Roads %>% st_set_crs(st_crs("+proj=utm +zone=37"))
Roads<-st_transform(Roads, crs = "+proj=utm +zone=37")
Roads<-st_crop(Roads,y=c(xmin=(min_x-1000), ymin=(min_y-1000), ymax=(min_y+y_length+1000),xmax=(min_x+x_length+1000))) # crop to size of kriging output
Roads<-st_transform(Roads, crs = "+proj=longlat +datum=WGS84")
Roads<- as_Spatial(Roads)




ggplot(data = eta[[9]], aes(x=Lon, y=Lat, color=V1)) + 
  geom_polygon(data = Bale_shapefile, aes(x = long, y = lat, group = group), colour = "black", fill = NA, size=1)+
  scale_colour_gradient2(low  = "blue",mid="white",high="red", name = "", midpoint=0) +theme(legend.position="bottom")+
  labs(y="",x="")+
  geom_point(size=3)+
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), plot.title = element_text(hjust = 0.5, size=30,face="bold"),
        legend.text = element_text(size=15,face="bold"),
        legend.position = "bottom") +
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                                    style = ggspatial::north_arrow_nautical(
                                      fill = c("black", "white"),
                                      line_col = "black")) +
  geom_path(data=Roads,aes(x = long, y = lat, group = group), colour="darkgrey",linetype="dotted", size=1)

V2<-ggplot(data = eta.normal, aes(x=Longitude, y=Latitude, color=V2)) + 
  geom_polygon(data = Bale_shapefile, aes(x = long, y = lat, group = group), colour = "black", fill = NA, size=1)+
  scale_colour_gradient2(low  = "blue",mid="white",high="red", name = "", midpoint=0) +theme(legend.position="bottom")+
  labs(y="",x="")+
  geom_point(size=3)+
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), plot.title = element_text(hjust = 0.5, size=30,face="bold"),
        legend.text = element_text(size=15,face="bold"),
        legend.position = "bottom")+
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                                    style = ggspatial::north_arrow_nautical(
                                      fill = c("black", "white"),
                                      line_col = "black")) +
  geom_path(data=Roads,aes(x = long, y = lat, group = group), colour="darkgrey",linetype="dotted", size=1)


