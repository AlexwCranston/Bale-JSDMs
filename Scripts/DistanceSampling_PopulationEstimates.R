# ---- Load Libraries

library(Distance)
library(dplyr)

# ---- Set Seed

set.seed(1)

# ---- Load Data

filenames.da <- c("Distance Sampling/Data/2010_09_16_to_2010_11_14_LateOther2010.csv",
                  "Distance Sampling/Data/2011_02_16_to_2011_11_14_LateOther2011.csv",
                  "Distance Sampling/Data/2012_06_15_to_2012_11_14_Wet2012_LateOther2012_Combined.csv",
                  "Distance Sampling/Data/2013_06_15_to_2013_11_14_Wet2013_LateOther2013_Combined.csv",
                  "Distance Sampling/Data/2015_09_16_to_2015_11_14_LateOther2015.csv",
                  "Distance Sampling/Data/2016_06_15_to_2016_11_14_Wet2016_LateOther2016_Combined.csv",
                  "Distance Sampling/Data/2017_09_16_to_2017_11_14_LateOther2017.csv",
                  "Distance Sampling/Data/2019_09_16_to_2019_11_14_LateOther2019.csv",
                  "Distance Sampling/Data/2020_09_16_to_2020_11_14_LateOther2020.csv",
                  "Distance Sampling/Data/2021_09_16_to_2021_11_14_LateOther2021.csv") # These are all the dataframes for each survey we're interested in

da.wet.list <- lapply(filenames.da, read.csv)  # Create list of all dataframes for all years
View(da.wet.list[[5]])

effort.da.wet <- c(110.938,
               132.814,
               132.264,
               131.252,
               113.941,
               91.797,
               69.521,
               67.584,
               68.72,
               97.026) # Effort expended in each survey in km

da.wet.list <- lapply(da.wet.list, function(df) {
  df <- df %>%
    filter(!is.na(PerpDistance)) %>%
    filter(!is.na(Total)) %>%
    filter(!is.na(Species))
  return(df)
})   # Drop any records where we don't have Perpendicular Distance, total or species

## Area of study area is 3442 hectares 


conversion.factor <- convert_units("meter", # units for perpendicular distance
                                   "kilometer", # units for effort - i.e. length of transects
                                   "hectare") # units for prediction area, i.e. Gaysay study area

cutpoints<-c(0,5,10,15,20,30,40,50,65,80,100,150,200,250)


# Add object column (required for dht2)

da.wet.list <- lapply(da.wet.list, function(df) {
  df$object <- NA
  df$object[!is.na(df$PerpDistance)] <- 1:sum(!is.na(df$PerpDistance))
  return(df)
})

hist(da.wet.list[[1]] %>% filter(Species=="MN") %>% pull("Total") %>% as.numeric()) # There is one Mountain Nyala record that is an order of magnitude larger than all other observations (500, presumably a typo, should be 50). We drop this.

da.wet.list[[1]] <- da.wet.list[[1]] %>% 
  filter(!(Species == "MN" & Total == 500))

hist(da.wet.list[[1]] %>% filter(Species=="MN") %>% pull("Total") %>% as.numeric())


# Rename columns to meet Distance package expectations

da.wet.list <- lapply(da.wet.list, function(df) {
  df <- df %>%
    rename(
      distance = PerpDistance,
      Sample.Label = Transect,
      Region.Label = ï..Stratum,
      size = Total
    )
  return(df)
}
)


da.wet.list <- lapply(da.wet.list, function(df) {
  df$Species <- as.factor(df$Species)
  return(df)
})

## Detection functions for wild ungulates ####

da.wet.wildONLY.list <- lapply(da.wet.list, function (df) {
  df <- df %>% filter(Species  %in% c("MN","WH","RB","BB")) # Mountain Nyala, Warthog, Reedbuck, Bushbuck
  df <- droplevels(df) 
  return(df)                    
}) 


# Modelling detection probability

# Define the function to run both models and calculate AIC
run_models_and_select_best <- function(df) {
  # Run the Hazard-rate model
  hr_model <- ds(data = df, transect = "line", key = "hr", truncation = 250, adjustment = "cos", convert_units = conversion.factor, formula = ~ Species)

  # Run the Half-normal model
  hn_model <- ds(data = df, transect = "line", key = "hn", truncation = 250, adjustment = "cos", convert_units = conversion.factor, formula = ~Species)
  
  # Calculate AIC for both models
  hr_aic <- hr_model$ddf$criterion
  hn_aic <- hn_model$ddf$criterion
  
  # Compare AIC and select the best model
  if (hr_aic < hn_aic) {
    best_model <- hr_model
  } else {
    best_model <- hn_model
  }
  
  # Return a list with both models and the best model based on AIC
  return(list(
    hr_model = hr_model,
    hn_model = hn_model,
    hr_aic = hr_aic,
    hn_aic = hn_aic,
    best_model = best_model
  ))
}

wild.results_list <- lapply(da.wet.wildONLY.list, run_models_and_select_best)

# Let's look at the results

wild.results_list[[9]] # Hazard-rate key function is the better fit every year except 2020, model 9

# 

plot(wild.results_list[[8]][["best_model"]], main="Detection function", breaks=cutpoints)

for (i in 1:10) {
  result.gof  <- gof_ds(wild.results_list[[i]][["best_model"]]) ## All years pass goodness of fit test except 2020, model 8
  p_value <- result.gof$dsgof$CvM$p
  print(paste0("Model ", i, " p-value: ", p_value)) }


plot(wild.results_list[[1]][[5]], showpoints=FALSE, main="Gaysay line transects\nspecies as covariate", breaks=cutpoints)
add.df.covar.line(wild.results_list[[1]][["best_model"]], data=da.wet.wildONLY.list[[1]] %>% filter(Species=="MN"), 
                 lwd=3, lty=1, col="blue")
add.df.covar.line(wild.results_list[[1]][["best_model"]], data=da.wet.wildONLY.list[[1]] %>% filter(Species=="BB"), 
                  lwd=3, lty=1, col="darkgreen")
add.df.covar.line(wild.results_list[[1]][["best_model"]], data=da.wet.wildONLY.list[[1]] %>% filter(Species=="RB"), 
                  lwd=3, lty=1, col="brown")
add.df.covar.line(wild.results_list[[1]][["best_model"]], data=da.wet.wildONLY.list[[1]] %>% filter(Species=="WH"), 
                  lwd=3, lty=1, col="salmon")


estimates.wild.list <- lapply(seq_along(wild.results_list), function(i) {
  dht2(ddf=wild.results_list[[i]][["best_model"]], flatfile=da.wet.wildONLY.list[[i]],
       strat_formula = ~Species, convert_units = conversion.factor,
       stratification = "object")
})

abundance_wild_estimates <- data.frame(
  Year = integer(),
  Species = character(),
  Estimate = numeric(),
  SE = numeric(),
  CV = numeric(),
  LCI = numeric(),
  UCI = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(estimates.wild.list)) {
  Year <- c(2010,2011,2012,2013,2015,2016,2017,2019,2020,2021)
  data <- as.data.frame(estimates.wild.list[[i]])
  for (j in 1:nrow(data)) {
    species <- data$Species[j]
    estimate <- data$Abundance[j]
    se <- data$Abundance_se[j]
    cv <- data$Abundance_CV[j]
    lci <- data$LCI[j]
    uci <- data$UCI[j]
    abundance_wild_estimates <- rbind(abundance_wild_estimates, data.frame(
      Year = Year[i],  
      Species = species,
      Estimate = estimate,
      SE = se,
      CV = cv,
      LCI = lci,
      UCI = uci
    ))
  }
}   

abundance_wild_estimates$DomesticY <- "N"
abundance_wild_estimates <- abundance_wild_estimates %>%
  mutate(Species = recode_factor(Species, "Total" = "TotalWild"))

## Detection functions for wild ungulates ####

da.wet.domesticONLY.list <- lapply(da.wet.list, function (df) {
  df <- df %>% filter(Species  %in% c("CT","HO","ST")) # Cattle, Horses, Sheep and Goats
  df <- droplevels(df) 
  return(df)                    
}) 


# Modelling detection probability

domestic.results_list <- lapply(da.wet.domesticONLY.list, run_models_and_select_best)

# We get a warning in certain models (models 3 and 5) due to spike in data at zero distance, i.e. on the transect - for these we force select of half normal as the best model to avoid overestimating abundance
# Warning in ddf.ds(dsmodel = dsmodel, data = data, meta.data = meta.data,  : 
# Estimated hazard-rate scale parameter close to 0 (on log scale). Possible problem in data (e.g., spike near zero distance).
#  

domestic.results_list[[3]][["best_model"]] <- domestic.results_list[[3]][["hn_model"]] 
domestic.results_list[[5]][["best_model"]] <- domestic.results_list[[5]][["hn_model"]] 


# Let's look at the results

domestic.results_list[[4]] # Hazard-rate key function is the better fit every year except 2011, 2012 & 2015

# 

plot(domestic.results_list[[5]][["best_model"]], main="Detection function", breaks=cutpoints)

for (i in 1:10) {
  result.gof  <- gof_ds(domestic.results_list[[i]][["best_model"]]) ## All years pass goodness of fit test except 2012
  p_value <- result.gof$dsgof$CvM$p
  print(paste0("Model ", i, " p-value: ", p_value)) }


plot(domestic.results_list[[5]][[5]], showpoints=FALSE, main="Gaysay line transects\nspecies as covariate", breaks=cutpoints)
add.df.covar.line(domestic.results_list[[5]][["best_model"]], data=da.wet.domesticONLY.list[[5]] %>% filter(Species=="CT"), 
                  lwd=3, lty=1, col="blue")
add.df.covar.line(domestic.results_list[[5]][["best_model"]], data=da.wet.domesticONLY.list[[5]] %>% filter(Species=="HO"), 
                  lwd=3, lty=1, col="darkgreen")
add.df.covar.line(domestic.results_list[[5]][["best_model"]], data=da.wet.domesticONLY.list[[5]] %>% filter(Species=="ST"), 
                  lwd=3, lty=1, col="salmon")


estimates.domestic.list <- lapply(seq_along(domestic.results_list), function(i) {
  dht2(ddf=domestic.results_list[[i]][["best_model"]], flatfile=da.wet.domesticONLY.list[[i]],
       strat_formula = ~Species, convert_units = conversion.factor,
       stratification = "object")
})

abundance_domestic_estimates <- data.frame(
  Year = integer(),
  Species = character(),
  Estimate = numeric(),
  SE = numeric(),
  CV = numeric(),
  LCI = numeric(),
  UCI = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(estimates.domestic.list)) {
  Year <- c(2010,2011,2012,2013,2015,2016,2017,2019,2020,2021)
  data <- as.data.frame(estimates.domestic.list[[i]])
  for (j in 1:nrow(data)) {
    species <- data$Species[j]
    estimate <- data$Abundance[j]
    se <- data$Abundance_se[j]
    cv <- data$Abundance_CV[j]
    lci <- data$LCI[j]
    uci <- data$UCI[j]
    abundance_domestic_estimates <- rbind(abundance_domestic_estimates, data.frame(
      Year = Year[i],  
      Species = species,
      Estimate = estimate,
      SE = se,
      CV = cv,
      LCI = lci,
      UCI = uci
    ))
  }
}   

abundance_domestic_estimates$DomesticY <- "Y"

abundance_domestic_estimates <- abundance_domestic_estimates %>%
  mutate(Species = recode_factor(Species, "Total" = "TotalDomestic"))

# Export the final estimates

abundance_estimates <- rbind(abundance_domestic_estimates, abundance_wild_estimates)

write.csv(abundance_estimates, "Processed Data/AbundanceEstimates.csv", row.names = FALSE)

# Let's have a quick look at the data

ggplot(abundance_estimates %>% filter(Species == "CT"|Species == "ST"|Species == "HO"), aes(x = Year, y = Estimate, color = Species)) +
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
  labs(
    title = "Abundance of Livestock vs Wild Ungulates across Years",
    x = "Year",
    y = "Estimated Abundance",
    caption = "Error bars represent 95% confidence interval (LCI, UCI)"
  ) +
  theme_minimal() +  # Use a minimal theme for clean visualization
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )

