library(tidyverse)
library(Hmsc)
library(corrplot)
library(reshape2)


set.seed(1)

# ---- Load in Data

data <- read.csv("Processed Data/AbundanceEstimates.csv")

data <- data %>% filter(Species != "TotalDomestic" & Species != "TotalWild")

# Two models had poor fit, we set the Estimates of these to NA 
data <- data %>%
  mutate(Estimate = ifelse(DomesticY == "Y" & Year == 2012, NA, Estimate),
         SE = ifelse(DomesticY == "Y" & Year == 2012, NA, SE),
         CV = ifelse(DomesticY == "Y" & Year == 2012, NA, CV),
         LCI = ifelse(DomesticY == "Y" & Year == 2012, NA, LCI),
         UCI = ifelse(DomesticY == "Y" & Year == 2012, NA, UCI))
data <- data %>%
  mutate(Estimate = ifelse(DomesticY == "N" & Year == 2019, NA, Estimate),
         SE = ifelse(DomesticY == "N" & Year == 2019, NA, SE),
         CV = ifelse(DomesticY == "N" & Year == 2019, NA, CV),
         LCI = ifelse(DomesticY == "N" & Year == 2019, NA, LCI),
         UCI = ifelse(DomesticY == "N" & Year == 2019, NA, UCI))




ggplot(data %>% filter(Species %in% c("TotalDomestic", "TotalWild")), aes(x = Year, y = Estimate, color = Species)) +
  geom_point(size = 3) + 
  geom_line(data = data %>% 
              filter(Species %in% c("TotalDomestic", "TotalWild")) %>%
              filter(!is.na(Estimate)), size = 1) +
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




data_model <- data %>%
  select(Year, Species, Estimate) %>%  # Select only necessary columns
  pivot_wider(names_from = Species, values_from = Estimate)  # Pivot


## Modelling 

Y <- data_model %>%
  select(CT, ST, HO, MN, BB, RB, WH)

Y  <- Y  %>%
  rename(
    Cattle = CT,
    Sheep...goats = ST,
    Horse = HO,
    Mountain.Nyala = MN,
    Menelik.s.Bushbuck = BB,
    Bohor.Reedbuck = RB,
    Warthog = WH
  )

xy = as.matrix(data_model %>%
                 select(Year))

TrData<- read.csv("Trait Data/Trait Data.csv")

rownames(TrData) <- TrData$ï..Species

TrData <- TrData[colnames(Y), ]

TrData <- TrData %>%
  dplyr::select("Domestic") 
  
TrFormula=~Domestic

studyDesign = data.frame(ID = as.factor(1:10))
rownames(xy)=studyDesign[,1]
rL = HmscRandomLevel(sData = xy)

XFormula = ~1

XData<- as.data.frame(c(1,1,1,1,1,1,1,1,1,1))

m = Hmsc(Y = Y, XData = XData,
         XFormula = XFormula,
         TrData = TrData,
         TrFormula = TrFormula,
         studyDesign = studyDesign,
         ranLevels = list(ID = rL), distr = "normal")


# Save the defined models

localDir = "."
modelDir = file.path(localDir, "models")
if(!dir.exists(modelDir)) dir.create(modelDir)  


save(m, file = file.path(modelDir, "unfitted_Time_Series_model.RData"))


## Check these models actually run without errors before moving on

sampleMcmc(m,samples=2)

nChains = 4
samples = 200 # We have a modest number of samples because this would make the model objects too large and lead to computationally intensive post-processing of results which would be slow and unnecessary
thin = 100 # We have a high thinning interval because we need the chains to run for a long time in order to achieve convergence
transient = round(0.5*samples*thin)
m.time <- sampleMcmc(m, thin =   thin, 
                                 samples = samples,transient = transient,
                                 nChains = nChains, nParallel = 4)
# Save Models

save(m.time, file = file.path(modelDir,paste("TimeSeriesmodel_thin_", as.character(thin),
                                             "_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             ".Rdata",sep = "")))

filename = file.path(modelDir,paste("TimeSeriesmodel_thin_", as.character(thin),
                              "_samples_", as.character(samples),
                              "_chains_",as.character(nChains),
                              ".Rdata",sep = ""))


load(filename)

# Evaluate Convergence

mpost = convertToCodaObject(m.time, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
psrf_beta = gelman.diag(mpost$Beta, multivariate=FALSE)$psrf


tmp = mpost$Omega[[1]]
psrf_Omega = gelman.diag(tmp, multivariate = FALSE)$psrf


## Check Associations

OmegaCor.time = computeAssociations(m.time)

supportLevel = 0.95

toPlot.time = ((OmegaCor.time[[1]]$support > supportLevel) + (OmegaCor.time[[1]]$support < (1-supportLevel)) > 0) * OmegaCor.time[[1]]$mean

melted_matrix.time <- melt(toPlot.time)

corrplot.time <- ggplot(melted_matrix.time, aes(Var1, Var2, fill = value)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "gray", limits = c(-1,1)) + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Correlation Plot - Time Series", x = "", y = "")) + 
    coord_fixed() +
    guides(fill = "none")
  return(plot)



