library(tidyverse)
library(abind)
library(mapview)
library(sf)
library(ids)
library(foreach) # parallel computing
library(doParallel) # parallel computing
library(kknn) # categorical knn
library(raster)
library(patchwork)


set.seed(1)

da <- read.csv("Raw Data/qrySpacial.csv")


# Data Cleaning -----------------------------------------------------------

focal_species <- c("Cattle",
                  "Menelik's Bushbuck",
                  "Warthog",
                  "Mountain Nyala",
                  "Sheep & goats",
                  "Donkey",
                  "Grey Duiker",
                  "Bohor Reedbuck",
                  "Horse",
                  "Klipspringer",
                  "Mixed Lstock(cows and shoats)",
                  "Domestic dog",
                  "Domestic cat")

Bale_mountains_crs <- "EPSG:20137 - Adindan / UTM zone 37N"

ungulate.da <- da %>%
  filter(SpeciesName %in% focal_species, # Restrict to species of interest
         Year != 1900) %>% # Restrict to correct data
  mutate(SpeciesName =  case_when(SpeciesName == "Mixed Lstock(cows and shoats)" ~ "Mixed Livestock",
                                  TRUE  ~  SpeciesName)) %>% # Simplify name of one entry
  na.omit() # Remove any rows with NAs

# Check for duplicate date rows
ungulate.da %>% janitor::get_dupes() # We have some duplicated data (n=13), let's delete these

ungulate.da <- ungulate.da[!duplicated(ungulate.da),] # dplyr::select all rows that aren't duplicates

ggplot() + geom_bar(data = ungulate.da, aes(as.factor(x=Year)))+
  xlab("Year") +ylab("Count")

ungulate.da$Month <- as.factor(ungulate.da$Month)

ggplot() + geom_bar(data = ungulate.da, aes(as.factor(x=Month)))+
  scale_x_discrete(labels=c("January","February","March","April","May","June", "July", "August","September","October","November","December"),drop=FALSE)  +
  xlab("Year") +ylab("Count")

year.plots<-c()
for (i in 1:14){
ungulate.da.year<-ungulate.da %>% filter(Year==(i+2008))
year.plots[[i]]<-ggplot() + geom_bar(data = ungulate.da.year, aes((x=Month))) +
 xlab("Month") +ylab("Count") + ggtitle(as.character(head(ungulate.da.year$Year))) +
  scale_x_discrete(labels=c("January","February","March","April","May","June", "July", "August","September","October","November","December"),drop=FALSE)  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

(year.plots[[1]]+year.plots[[2]])/(year.plots[[3]]+year.plots[[4]])
(year.plots[[5]]+year.plots[[6]])/(year.plots[[7]]+year.plots[[8]])
(year.plots[[9]]+year.plots[[11]])/(year.plots[[12]]+year.plots[[13]])
year.plots[[14]]

# There's data from end of wet season surveys in 2009, 2010, 2011, 2012, 2013, 2015, 2016, 2017, 2019, 2020, 2021
# There's data from end of dry season surveys from 2014, 2021 and 2022
# 14 survey datasets in total


ungulate.da <-ungulate.da %>% mutate(Survey = case_when(
  Month %in% c(3,4) ~ "Dry",
  Month %in% c(9,10,11) ~ "Wet",
  Month %in% c(1,2,5,6,7,8,12) ~ "Incidental"))


ggplot(data=ungulate.da, aes(as.factor(x=Year), fill=Survey))+
  geom_bar(stat="count", position = position_dodge2(width = 0.9, preserve = "single")) # This just lets us see all data from all years on one graph 

ungulate.da<- ungulate.da %>% dplyr::select(-Survey)
  
ungulate.da <- ungulate.da %>% filter(Month == 3 | Month == 4 | Month ==9| Month ==10 | Month == 11) # Select only data from months where samples took place

ungulate.sf <- sf::st_as_sf(ungulate.da,
                     coords = c("UTMX", "UTMY"),
                     crs = Bale_mountains_crs)

mapview::mapview(ungulate.sf) # At this point, there are clearly some erroneous data points

Bale <- read_sf("Shapefiles/BaleMountainsSHP_WDPA_WDOECM_Feb2023_Public_2281_shp/WDPA_WDOECM_Feb2023_Public_2281_shp-polygons.shp")
Bale <- st_transform(Bale, crs = Bale_mountains_crs)

ungulate.sf <- st_filter(ungulate.sf, Bale) # this step loses 222 data points

ungulate.sf<-ungulate.sf %>% 
  rename(
    Transect = ï..TransectName,
    SpCod = SpeciesCode,
    SpName = SpeciesName
  )


st_write(ungulate.sf, "Processed Data/Bale_surveyobservationsCLEAN_2009_2022.shp", append = FALSE)

mapview::mapview(Bale) +
  mapview::mapview(ungulate.sf)

View(ungulate.sf %>% group_by(SpName) %>% summarise(n = n())) # Look at the numbers of samples by species and sex


ungulate.sf_UTM <- unlist(st_geometry(ungulate.sf)) %>%
  matrix(ncol = 2, byrow = TRUE) %>%
  as_tibble() %>%
  setNames(c("UTMX", "UTMY")) ##### Convert to normal dataframe

ungulate.sf_data <- st_drop_geometry(ungulate.sf)
ungulate.sf_final <- cbind(ungulate.sf_data, ungulate.sf_UTM) #7322 of original 9378 data points preserved


# Create a grid for habitat interpolation

#Specify grid cell size
# dx is the width of a grid cell in meters

dx<-100
dy <- dx #  # dy is the height of a grid cell in meters

# calculate grid width in pixels
width_in_pixels <- ceiling((st_bbox(ungulate.sf)["xmax"] - st_bbox(ungulate.sf)["xmin"]) / dx)


# calculate the height in pixels of the resulting grid
height_in_pixels <- floor( (st_bbox(ungulate.sf)["ymax"] - 
                              st_bbox(ungulate.sf)["ymin"]) / dy)

habitat.grid <- st_make_grid(ungulate.sf, 
                     cellsize = dx,
                     n = c(width_in_pixels, height_in_pixels),
                     what = "centers"
)


plot(habitat.grid)

## Make Grid for binning count data

x <- seq(from = 579700, to = 585000, by = 100)
y <- seq(from = 780000, to = 790300, by = 100)
xy <- expand.grid(x, y)

colnames(xy) <- c("UTMX", "UTMY")

grid <- sf::st_as_sf(xy, coords = c("UTMX", "UTMY"), crs = "EPSG:20137 - Adindan / UTM zone 37N")

mapview::mapview(Bale) + mapview::mapview(grid)

bounds <- data.frame(
  xl = (xy$UTMX - 50),
  xu = (xy$UTMX + 50),
  yl = (xy$UTMY - 50),
  yu = (xy$UTMY + 50),
  xy
)

bounds$id <- random_id(5616, 4)

# Habitat

ungulate.sf_final %>%
  group_by(Habitat) %>%
  summarise(no_rows = length(Habitat))

unique(ungulate.sf_final$Habitat)

ungulate.sf_final$Habitat <- as.factor(ungulate.sf_final$Habitat)

levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "2"] <- "Agriculture" # 2 is originally Agriculture
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "5"] <- "Woodland" # 5 is originally Coffee
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "6"] <- "Woodland" # 6 is originally Closed Woodland
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "7"] <- "Woodland" # 7 is originally Erica Forest
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "8"] <- "Shrubs/Heaths" # 8 is originally Erica Shrub
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "9"] <- "Woodland" # 9 is originally Glade
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "10"] <- "Grassland" # 10 is originally Gaysay grassland
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "11"] <- "Shrubs/Heaths" # 11 is originally Grassland-Shrub mixed
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "12"] <- "Grassland" #  12 is originally Grazing land
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "13"] <- "Shrubs/Heaths" # 13 is originally Helicrisum
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "15"] <- "Woodland" # 15 is originally Open Woodland
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "17"] <- "Shrubs/Heaths" # 17 is originally Shrub
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "19"] <- "Wetland" # 19 is originally Wetland 
levels(ungulate.sf_final$Habitat)[levels(ungulate.sf_final$Habitat) == "20"] <- "Shrubs/Heaths" # 20 is originally Mixed Shrub
levels(ungulate.sf_final$Habitat) # We go from 14 different habitats with unclear definitions to 5 clear categories

ungulate.sf_final <- ungulate.sf_final %>% filter(Habitat != "Agriculture") # Get rid of two points from agriculture
ungulate.sf_final <- droplevels(ungulate.sf_final)



## We will interpolate habitat over the whole grid but do it separately for dry season data and wet season data



# Wet season first 

ALL_WetSeason <- ungulate.sf_final %>% filter(Month == 9 | Month == 10 |Month == 11) # 5008 data points


training_ALL_WetSeason <- data.frame(habitat = ALL_WetSeason$Habitat, 
                                   UTMX = ALL_WetSeason$UTMX, 
                                   UTMY = ALL_WetSeason$UTMY) ## This is the data that the Kriging Interpolation will train on

habitats_result_ALL_WetSeason <- data.frame(habitat = as.factor(NA), 
                                             UTMX = st_coordinates(habitat.grid)[, 1], 
                                             UTMY = st_coordinates(habitat.grid)[, 2]) # This is the empty grid that will be populated with interpolated data

habitats_kknn <- kknn::kknn(habitat ~ ., 
                            train = training_ALL_WetSeason, 
                            test = habitats_result_ALL_WetSeason, 
                            kernel = "epanechnikov", 
                            k = 3) 

habitats_result_ALL_WetSeason<-habitats_result_ALL_WetSeason %>%
  # extract the interpolated dialect at each grid cell with the 
  # kknn::fitted function
  mutate(habitat = fitted(habitats_kknn),
         # only retain the probability of the interpolated habitat,
         # discard the others
         prob = apply(habitats_kknn$prob, 
                      1, 
                      function(x) max(x)))


ggplot(data = habitats_result_ALL_WetSeason,aes(x=UTMX,y=UTMY, col=habitat))+geom_point() 

hist(habitats_result_ALL_WetSeason$prob) # Most points have a high probability of being correct
mean(habitats_result_ALL_WetSeason$prob)

library(raster)
habitats_result_ALL_WetSeason$habitat.code = as.numeric(habitats_result_ALL_WetSeason$habitat) # Add a column changing habitat into a number



habitats_result_ALL_WetSeason.raster <- rasterFromXYZ(habitats_result_ALL_WetSeason[,c("UTMX","UTMY","habitat.code")]) # convert into a raster

habitats_result_ALL_WetSeason.raster[] = factor(levels(habitats_result_ALL_WetSeason$habitat)[habitats_result_ALL_WetSeason.raster[]]) # Converts the numbers to match the original dataset

# 1 is Grassland
# 2 is Shrubs/Heath
# 3 is Wetland
# 4 is Woodland

## Then Dry season

ALL_DrySeason <- ungulate.sf_final %>% filter(Month == 3| Month == 4) # 2312 data points


train_ALL_DrySeason <- data.frame(habitat = ALL_DrySeason$Habitat, 
                                  UTMX = ALL_DrySeason$UTMX, 
                                  UTMY = ALL_DrySeason$UTMY)

habitats_result_ALL_DrySeason <- data.frame(habitat = as.factor(NA), 
                                            UTMX = st_coordinates(habitat.grid)[, 1], 
                                            UTMY = st_coordinates(habitat.grid)[, 2]) 

habitats_kknn <- kknn::kknn(habitat ~ ., 
                            train = train_ALL_DrySeason, 
                            test = habitats_result_ALL_DrySeason, 
                            kernel = "epanechnikov", 
                            k = 3)

habitats_result_ALL_DrySeason<-habitats_result_ALL_DrySeason %>%
  # extract the interpolated dialect at each grid cell with the 
  # kknn::fitted function
  mutate(habitat = fitted(habitats_kknn),
         # only retain the probability of the interpolated habitat,
         # discard the others
         prob = apply(habitats_kknn$prob, 
                      1, 
                      function(x) max(x)))


ggplot(data = habitats_result_ALL_DrySeason,aes(x=UTMX,y=UTMY, col=habitat))+geom_point() 

hist(habitats_result_ALL_DrySeason$prob) # Most points have a high probability of being correct
mean(habitats_result_ALL_DrySeason$prob)


habitats_result_ALL_DrySeason$habitat.code = as.numeric(habitats_result_ALL_DrySeason$habitat)

habitats_result_ALL_DrySeason.raster <- rasterFromXYZ(habitats_result_ALL_DrySeason[,c("UTMX","UTMY","habitat.code")])

habitats_result_ALL_DrySeason.raster[] = factor(levels(habitats_result_ALL_DrySeason$habitat)[habitats_result_ALL_DrySeason.raster[]])

plot(habitats_result_ALL_DrySeason.raster)


## We'll also use elevation as a predictor


elevation <- raster("Rasters/NASA Shuttle Radar Topography_reprojected100m.tif")


## Finally, we will include distance from the B90 road, which runs through the Gaysay area

road <- read_sf("Shapefiles/Roads/Ethiopia_roads_network_clipped.shp")
st_crs(road) <- Bale_mountains_crs
plot(road)


#Let's export the final "clean" data 

write.csv(ungulate.sf_final, "Processed Data/BaleMountainsNP_all survey observations_2009_2022.csv", row.names = FALSE)


# Now we're going to filter points into different seasons and export them as standalone data frames

###  2009 Wet Season

good_points_2009_WetSeason <- ungulate.sf_final %>% filter(Year == 2009 & Month == 9 | Year == 2009 & Month == 10 | Year == 2009 & Month == 11) # 96 data points

for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2009_WetSeason)) {
    x <- good_points_2009_WetSeason[j, 11]
    y <- good_points_2009_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2009_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}

good_points_2009_WetSeason <- good_points_2009_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2009_WetSeason <- unique(good_points_2009_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2009_WetSeason <- pivot_wider(aggregated.all.species_2009_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2009_WetSeason$id)) # Confirm that each grid space only has one row

lon <- final.all.species_2009_WetSeason$UTMX.1
lat <- final.all.species_2009_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2009_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2009_WetSeason$elevation <- raster::extract(elevation,samples)

# Distance to Road for each point
points <- st_as_sf(final.all.species_2009_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2009_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)

write.csv(final.all.species_2009_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2009_WetSeason_res100m.csv", row.names = FALSE)



###  2010_WetSeason


good_points_2010_WetSeason <- ungulate.sf_final %>% filter(Year == 2010 & Month == 9 | Year == 2010 & Month == 10 | Year == 2010 & Month == 11)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2010_WetSeason)) {
    x <- good_points_2010_WetSeason[j, 11]
    y <- good_points_2010_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2010_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2010_WetSeason <- good_points_2010_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2010_WetSeason <- unique(good_points_2010_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2010_WetSeason <- pivot_wider(aggregated.all.species_2010_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2010_WetSeason$id)) # Confirm that each grid space only has one row



lon <- final.all.species_2010_WetSeason$UTMX.1
lat <- final.all.species_2010_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2010_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2010_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2010_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2010_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)



write.csv(final.all.species_2010_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2010_WetSeason_res100m.csv", row.names = FALSE)


### 2011_WetSeason


good_points_2011_WetSeason <- ungulate.sf_final %>% filter(Year == 2011 & Month == 9 | Year == 2011 & Month == 10 | Year == 2011 & Month == 11)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2011_WetSeason)) {
    x <- good_points_2011_WetSeason[j, 11]
    y <- good_points_2011_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2011_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2011_WetSeason <- good_points_2011_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2011_WetSeason <- unique(good_points_2011_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2011_WetSeason <- pivot_wider(aggregated.all.species_2011_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2011_WetSeason$id)) # Confirm that each grid space only has one row

lon <- final.all.species_2011_WetSeason$UTMX.1
lat <- final.all.species_2011_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2011_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2011_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2011_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2011_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)




write.csv(final.all.species_2011_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2011_WetSeason_res100m.csv", row.names = FALSE)


### 2012_WetSeason

good_points_2012_WetSeason <- ungulate.sf_final %>% filter(Year == 2012 & Month == 9 | Year == 2012 & Month == 10 | Year == 2012 & Month == 11)



for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2012_WetSeason)) {
    x <- good_points_2012_WetSeason[j, 11]
    y <- good_points_2012_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2012_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2012_WetSeason <- good_points_2012_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2012_WetSeason <- unique(good_points_2012_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2012_WetSeason <- pivot_wider(aggregated.all.species_2012_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2012_WetSeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2012_WetSeason$UTMX.1
lat <- final.all.species_2012_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2012_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2012_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2012_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2012_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)



write.csv(final.all.species_2012_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2012_WetSeason_res100m.csv", row.names = FALSE)


### 2013_WetSeason

good_points_2013_WetSeason <- ungulate.sf_final %>% filter(Year == 2013 & Month == 9 | Year == 2013 & Month == 10 | Year == 2013 & Month == 11)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2013_WetSeason)) {
    x <- good_points_2013_WetSeason[j, 11]
    y <- good_points_2013_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2013_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2013_WetSeason <- good_points_2013_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2013_WetSeason <- unique(good_points_2013_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2013_WetSeason <- pivot_wider(aggregated.all.species_2013_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2013_WetSeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2013_WetSeason$UTMX.1
lat <- final.all.species_2013_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2013_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2013_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2013_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2013_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)


write.csv(final.all.species_2013_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2013_WetSeason_res100m.csv", row.names = FALSE)


### 2014_DrySeason

good_points_2014_DrySeason <- ungulate.sf_final %>% filter(Year == 2014 & Month == 3 | Year == 2014 & Month == 4)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2014_DrySeason)) {
    x <- good_points_2014_DrySeason[j, 11]
    y <- good_points_2014_DrySeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2014_DrySeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2014_DrySeason <- good_points_2014_DrySeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2014_DrySeason <- unique(good_points_2014_DrySeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2014_DrySeason <- pivot_wider(aggregated.all.species_2014_DrySeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2014_DrySeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2014_DrySeason$UTMX.1
lat <- final.all.species_2014_DrySeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2014_DrySeason$habitat <- raster::extract(habitats_result_ALL_DrySeason.raster,samples)
final.all.species_2014_DrySeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2014_DrySeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2014_DrySeason$distance_from_road<-apply(distance, 1, FUN = min)


write.csv(final.all.species_2014_DrySeason, "Processed Data/BaleMountainsNP_MeanAbundance_2014_DrySeason_res100m.csv", row.names = FALSE)


## 2015_WetSeason

good_points_2015_WetSeason <- ungulate.sf_final %>% filter(Year == 2015 & Month == 9 | Year == 2015 & Month == 10 | Year == 2015 & Month == 11)

for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2015_WetSeason)) {
    x <- good_points_2015_WetSeason[j, 11]
    y <- good_points_2015_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2015_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2015_WetSeason <- good_points_2015_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2015_WetSeason <- unique(good_points_2015_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2015_WetSeason <- pivot_wider(aggregated.all.species_2015_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2015_WetSeason$id)) # Confirm that each grid space only has one row



lon <- final.all.species_2015_WetSeason$UTMX.1
lat <- final.all.species_2015_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2015_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2015_WetSeason$elevation <- raster::extract(elevation,samples)



# Distance to Road for each point
points <- st_as_sf(final.all.species_2015_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2015_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)



write.csv(final.all.species_2015_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2015_WetSeason_res100m.csv", row.names = FALSE)

## 2016_WetSeason

good_points_2016_WetSeason <- ungulate.sf_final %>% filter(Year == 2016 & Month == 9 | Year == 2016 & Month == 10 | Year == 2016 & Month == 11)

for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2016_WetSeason)) {
    x <- good_points_2016_WetSeason[j, 11]
    y <- good_points_2016_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2016_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2016_WetSeason <- good_points_2016_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2016_WetSeason <- unique(good_points_2016_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2016_WetSeason <- pivot_wider(aggregated.all.species_2016_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2016_WetSeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2016_WetSeason$UTMX.1
lat <- final.all.species_2016_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2016_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2016_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2016_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2016_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)


write.csv(final.all.species_2016_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2016_WetSeason_res100m.csv", row.names = FALSE)

### good_points_2017_WetSeason

good_points_2017_WetSeason <- ungulate.sf_final %>% filter(Year == 2017 & Month == 9 | Year == 2017 & Month == 10 | Year == 2017 & Month == 11)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2017_WetSeason)) {
    x <- good_points_2017_WetSeason[j, 11]
    y <- good_points_2017_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2017_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2017_WetSeason <- good_points_2017_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2017_WetSeason <- unique(good_points_2017_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2017_WetSeason <- pivot_wider(aggregated.all.species_2017_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2017_WetSeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2017_WetSeason$UTMX.1
lat <- final.all.species_2017_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2017_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2017_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2017_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2017_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)


write.csv(final.all.species_2017_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2017_WetSeason_res100m.csv", row.names = FALSE)

## 2019_WetSeason

good_points_2019_WetSeason <- ungulate.sf_final %>% filter(Year == 2019 & Month == 9 | Year == 2019 & Month == 10 | Year == 2019 & Month == 11)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2019_WetSeason)) {
    x <- good_points_2019_WetSeason[j, 11]
    y <- good_points_2019_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2019_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2019_WetSeason <- good_points_2019_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2019_WetSeason <- unique(good_points_2019_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2019_WetSeason <- pivot_wider(aggregated.all.species_2019_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2019_WetSeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2019_WetSeason$UTMX.1
lat <- final.all.species_2019_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2019_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2019_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2019_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2019_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)



write.csv(final.all.species_2019_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2019_WetSeason_res100m.csv", row.names = FALSE)

##  2020_WetSeason

good_points_2020_WetSeason <- ungulate.sf_final %>% filter(Year == 2020 & Month == 9 | Year == 2020 & Month == 10 | Year == 2020 & Month == 11)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2020_WetSeason)) {
    x <- good_points_2020_WetSeason[j, 11]
    y <- good_points_2020_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2020_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2020_WetSeason <- good_points_2020_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2020_WetSeason <- unique(good_points_2020_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2020_WetSeason <- pivot_wider(aggregated.all.species_2020_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2020_WetSeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2020_WetSeason$UTMX.1
lat <- final.all.species_2020_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2020_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2020_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2020_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2020_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)


write.csv(final.all.species_2020_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2020_WetSeason_res100m.csv", row.names = FALSE)

##  2021_DrySeason

good_points_2021_DrySeason <- ungulate.sf_final %>% filter(Year == 2021 & Month == 3 | Year == 2021 & Month == 4)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2021_DrySeason)) {
    x <- good_points_2021_DrySeason[j, 11]
    y <- good_points_2021_DrySeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2021_DrySeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2021_DrySeason <- good_points_2021_DrySeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2021_DrySeason <- unique(good_points_2021_DrySeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2021_DrySeason <- pivot_wider(aggregated.all.species_2021_DrySeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2021_DrySeason$id)) # Confirm that each grid space only has one row



lon <- final.all.species_2021_DrySeason$UTMX.1
lat <- final.all.species_2021_DrySeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2021_DrySeason$habitat <- raster::extract(habitats_result_ALL_DrySeason.raster,samples)
final.all.species_2021_DrySeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2021_DrySeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2021_DrySeason$distance_from_road<-apply(distance, 1, FUN = min)



write.csv(final.all.species_2021_DrySeason, "Processed Data/BaleMountainsNP_MeanAbundance_2021_DrySeason_res100m.csv", row.names = FALSE)
### 2021_WetSeason

good_points_2021_WetSeason <- ungulate.sf_final %>% filter(Year == 2021 & Month == 9 | Year == 2021 & Month == 10 | Year == 2021 & Month == 11)

for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2021_WetSeason)) {
    x <- good_points_2021_WetSeason[j, 11]
    y <- good_points_2021_WetSeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2021_WetSeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2021_WetSeason <- good_points_2021_WetSeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2021_WetSeason <- unique(good_points_2021_WetSeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2021_WetSeason <- pivot_wider(aggregated.all.species_2021_WetSeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2021_WetSeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2021_WetSeason$UTMX.1
lat <- final.all.species_2021_WetSeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2021_WetSeason$habitat <- raster::extract(habitats_result_ALL_WetSeason.raster,samples)
final.all.species_2021_WetSeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2021_WetSeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2021_WetSeason$distance_from_road<-apply(distance, 1, FUN = min)



write.csv(final.all.species_2021_WetSeason, "Processed Data/BaleMountainsNP_MeanAbundance_2021_WetSeason_res100m.csv", row.names = FALSE)

### 2022_DrySeason

good_points_2022_DrySeason <- ungulate.sf_final %>% filter(Year == 2022 & Month == 3 | Year == 2022 & Month == 4)


for (i in 1:5616) {
  grid <- bounds[i, ]
  for (j in 1:nrow(good_points_2022_DrySeason)) {
    x <- good_points_2022_DrySeason[j, 11]
    y <- good_points_2022_DrySeason[j, 12]
    if (grid$xl <= x & grid$xu >= x & grid$yl <= y & grid$yu >= y) {
      good_points_2022_DrySeason[j, 13:15] <- grid[, 5:7]
    }
  }
}
good_points_2022_DrySeason <- good_points_2022_DrySeason %>% dplyr::select(-c(11, 12)) # Remove old coordinates and just leave grid coordinates

aggregated.all.species_2022_DrySeason <- unique(good_points_2022_DrySeason[, c(3, 10:13)]) # Take unique records of species with unique totals from each grid space, i.e if one nyala was seen on two different days, only one record will be returned but if two was seen one day and one on a different day, both records will be returned

final.all.species_2022_DrySeason <- pivot_wider(aggregated.all.species_2022_DrySeason, names_from = SpName, values_from = Total, values_fill = 0, values_fn = mean)
summary(unique(final.all.species_2022_DrySeason$id)) # Confirm that each grid space only has one row


lon <- final.all.species_2022_DrySeason$UTMX.1
lat <- final.all.species_2022_DrySeason$UTMY.1
samples<- data.frame(lon,lat)
final.all.species_2022_DrySeason$habitat <- raster::extract(habitats_result_ALL_DrySeason.raster,samples)
final.all.species_2022_DrySeason$elevation <- raster::extract(elevation,samples)


# Distance to Road for each point
points <- st_as_sf(final.all.species_2022_DrySeason, coords = c("UTMX.1", "UTMY.1")) # Convert to sf object
st_crs(points)<-Bale_mountains_crs # Add crs
distance<-st_distance(points, road) # Calculate distance. Roads has three objects so it gives three distances. We want the shortest one.

final.all.species_2022_DrySeason$distance_from_road<-apply(distance, 1, FUN = min)



write.csv(final.all.species_2022_DrySeason, "Processed Data/BaleMountainsNP_MeanAbundance_2022_DrySeason_res100m.csv", row.names = FALSE)

