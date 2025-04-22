# ----Load in the libraries

library(tidyverse)
library(Hmsc)
library(sf)
library(tidysdm)

# ---- Set Seed

set.seed(1)

# ---- Load in the data


years <- c(2010, 2011, 2012, 2013, 2015, 2016, 2017, 2019, 2020, 2021) # List of years where we have wet season surveys

da.wet.list <- lapply(years, function(y) {
  read.csv(paste0("Processed Data/BaleMountainsNP_MeanAbundance_", y, "_WetSeason_res100m.csv"))
}) # Create list of all dataframes for each years

names(da.wet.list) <- paste0("da.wet.", years) # add names to the list of dataframes so we can identify them by year rather than by number

## Drop rows where we have missing predictor variables 

da.wet.list <- lapply(da.wet.list, function(df) {
  df <- df %>% drop_na(habitat) %>% 
    drop_na(elevation) %>% 
    drop_na(distance_from_road)
  return(df)
}
)
  
# Not all dataframes have records for mixed livestock, add an empty column if its absent so we can sum all domestic livestock together

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

## Thin the data to avoid the issues associated with clustering and reduce computational complexity


da.wet.list <-lapply(da.wet.list, function(df) {
  df<-st_as_sf(df, coords = c("UTMX.1", "UTMY.1"), crs = "EPSG:20137 - Adindan / UTM zone 37N")
  return(df)
}) # First convert to sf object


da.wet.list <-lapply(da.wet.list, function(sf) {
  sf <- thin_by_dist(sf, dist_min = 150)
  return(sf)
}) # Thin by distance, no points nearer than 150m


da.wet.list <- lapply(da.wet.list, function(sf) {
  sf <- st_transform(sf, crs= "EPSG:4326") 
  return(sf)
}) # Transform to the appropriate projection for long/lat coordinates

# Extract coordinates from data 
xy.wet.list <- lapply(da.wet.list, function(sf) {
  sf <- unlist(st_geometry(sf))%>% 
    matrix(ncol=2,byrow=TRUE) %>% 
    as_tibble() %>% 
    setNames(c("Lon","Lat"))
}
)

xy.wet.list<- lapply(xy.wet.list, function(df) {
  df <- as.data.frame(df)
  return(df)
}
)

# Now we have the coordinates as a separate object, we can drop the geometry from the data.

da.wet.list <- lapply(da.wet.list, function(sf) {
  sf <- st_drop_geometry(sf) 
  return(sf)
}) # Drop geometry

# Dependent Variables

Y.wet.list <- lapply(da.wet.list, function(df) {
  data.frame(df %>% dplyr::select(Mountain.Nyala,Bohor.Reedbuck, Menelik.s.Bushbuck, Warthog, Livestock))
  }
  ) # Select dependent variables, i.e. abundance data for 4 most important wild ungulate species plus livestock

# We create one dataframe that is purely presence absence and one where all absences are set to NA and we considering abundances only, i.e. conditional abundances (CAs) 

Y.wet.presab.list <- Y.wet.list 
Y.wet.presab.list <- lapply(Y.wet.presab.list, function(df) { 
  df[] <- lapply(df, function(x) ifelse(x > 0, 1, x))  # Replace values greater than 0 with 1
  return(df)
})

Y.wet.CA.list <- Y.wet.list
Y.wet.CA.list <- lapply(Y.wet.CA.list, function(df) { 
  df[] <- lapply(df, function(x) ifelse(x == 0, NA, x))  # Replace values equal to 0 with NA
  return(df)
})

# Independent Variables 

XData.wet.list <- lapply(da.wet.list, function(df) {
  data.frame(df %>% dplyr::select(habitat,elevation, distance_from_road))
}
)

XData.wet.list <- lapply(XData.wet.list, function(df) {
  df$habitat <- as.factor(df$habitat)
  return(df)
})

XData.wet.list <- lapply(XData.wet.list, function(df) {
  df <- df %>%  mutate(habitat = fct_recode(habitat,
                                         Grassland = "1",
                                         Shrubs_Heath = "2",
                                         Wetland = "3",
                                         Woodland = "4")) #Add the actual names of the habitats back in over the original code
})


## Trait Data

TrData <- read.csv("Trait Data/Trait Data.csv")
TrData <- TrData %>% filter(ï..Species %in% colnames(Y.wet.list[[1]])) # Drop unused species from the Trait data
rownames(TrData)<- TrData$ï..Species # Add the species names to the rownames
TrData<- TrData %>% dplyr::select(-1) # Drop column for Species names, HMSC expects these to be in the rownames 

TrFormula <- ~Domestic  # Formula for Traits

## Add Study design

studyDesign.wet.list <- lapply(XData.wet.list, function(df) {
  df <- data.frame(ID = as.factor(1:nrow(df)))
  return(df)
}
)

xy.wet.list <- lapply(seq_along(xy.wet.list), function(i) {
  rownames(xy.wet.list[[i]]) <- studyDesign.wet.list[[i]][, 1]
  return(xy.wet.list[[i]])
})

rL.wet.list <- lapply(xy.wet.list, function(df) {
  df <- HmscRandomLevel(sData = df)
  return(df)
}
) # Define random levels

## Define Models 

XFormula = ~ habitat + elevation + distance_from_road # Define formula for model 

# First conditional abundance models

models.CA.list <- lapply(seq_along(Y.wet.CA.list), function(i) { 
  model <- list() 
  model[[i]] <- Hmsc(Y = Y.wet.CA.list[[i]], 
                     XData = XData.wet.list[[i]], 
                     XFormula = XFormula, 
                     TrData = TrData, 
                     TrFormula = TrFormula, 
                     studyDesign = studyDesign.wet.list[[i]], 
                     ranLevels = list(ID = rL.wet.list[[i]]), 
                     distr = "lognormal poisson")
  return(model[[i]])
})

names(models.CA.list) <- paste0("models.CA.", years) # add names to the list of models so we can identify them by year rather than by number

# Next presence/absence models

models.presab.list <- lapply(seq_along(Y.wet.presab.list), function(i) { 
  model <- list() 
  model[[i]] <- Hmsc(Y = Y.wet.presab.list[[i]], 
                     XData = XData.wet.list[[i]], 
                     XFormula = XFormula, 
                     TrData = TrData, 
                     TrFormula = TrFormula, 
                     studyDesign = studyDesign.wet.list[[i]], 
                     ranLevels = list(ID = rL.wet.list[[i]]), 
                     distr = "probit")
  return(model[[i]])
})

names(models.presab.list) <- paste0("models.presab.", years) # add names to the list of dataframes so we can identify them by year rather than by number

 
# Save the defined models
 
localDir = "."
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)  

unfitted.models.list = list(CA.models = models.CA.list, presab.models = models.presab.list)
save(unfitted.models.list, file = file.path(modelDir, "unfitted_models.RData"))

## Check these models actually run without errors before moving on

for(i in 1:length(models.CA.list)){
  print(i)
  sampleMcmc(models.CA.list[[i]],samples=2)
} # First CA models

for(i in 1:length(models.presab.list)){
  print(i)
  sampleMcmc(models.presab.list[[i]],samples=2)
} # Then presab models
