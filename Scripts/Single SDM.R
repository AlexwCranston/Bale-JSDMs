# ----Load in the libraries

library(tidyverse)
library(Hmsc)
library(abind)
library(sf)
library(corrplot)
library(ggplot2)
library(gplots)

# ---- Set Seed

set.seed(1)

# ---- Load in the data

da<-read.csv("Processed Data/BaleMountainsNP_MeanAbundance_2020_WetSeason_res100m.csv")
da<-na.omit(da)


# Dependent Variables

Y=data.frame(da%>% select(c(5))) # Select dependent variables, i.e. abundance data for mountain nyala
Y.presab <- Y
Y.presab[Y.presab > 0] <- 1



# Independent Variables # habitat.CCI

XData <- da %>% select(7,14:17)

XData$habitat<-as.factor(XData$habitat)
XData$habitat.CCI <-as.factor(XData$habitat.CCI)



# Coordinates

xy = as.data.frame(cbind(da$UTMX.1,da$UTMY.1))

p1 <- st_as_sf(xy, coords = c("V1", "V2"), crs = "EPSG:20137 - Adindan / UTM zone 37N")
p2 <- st_transform(p1, crs= "EPSG:4326")
p2 <- unlist(st_geometry(p2))%>% 
  matrix(ncol=2,byrow=TRUE) %>% 
  as_tibble() %>% 
  setNames(c("Lon","Lat"))
xy<-as.data.frame(p2)

## Defining the models

studyDesign = data.frame(ID = as.factor(1:254))
rownames(xy)=studyDesign[,1]

XFormula = ~ habitat + elevation + distance_from_road + Cattle
XFormula.no.Cattle = ~ habitat + elevation + distance_from_road 

m.normal = Hmsc(Y = Y, XData = XData, XFormula = XFormula, studyDesign = studyDesign, distr = "lognormal poisson") 
m.normal.presab = Hmsc(Y = Y.presab, XData = XData, XFormula = XFormula, studyDesign = studyDesign, distr = "probit") 

m.no.Cattle = Hmsc(Y = Y, XData = XData, XFormula = XFormula.no.Cattle, studyDesign = studyDesign, distr = "lognormal poisson") 
m.no.Cattle.presab = Hmsc(Y = Y.presab, XData = XData, XFormula = XFormula.no.Cattle, studyDesign = studyDesign, distr = "probit") 

# Running the models
Start_time <- Sys.time()
nChains = 4
samples = 200 # We have a modest number of samples because this would make the model objects too large and lead to computationally intensive post-processing of results which would be slow and unnecessary
thin = 10 # We have a high thinning interval because we need the chains to run for a long time in order to achieve convergence
transient = round(0.5*samples*thin)
models <- list(m.normal,m.normal.presab,m.no.Cattle,m.no.Cattle.presab)
for (i in 1:4){
  models[[i]] = sampleMcmc(models[[i]], thin =   thin, 
                           samples = samples,transient = transient,
                           nChains = nChains, nParallel = 4) }
end_time <- Sys.time()
end_time-start_time 


# Results 

preds = computePredictedValues(models[[1]])
results<-evaluateModelFit(hM=models[[1]], predY=preds)
results<-as.data.frame(results)
rownames(results)<-colnames(Y)

preds.presab = computePredictedValues(models[[2]])
results.presab<-evaluateModelFit(hM=models[[2]], predY=preds.presab)
results.presab<-as.data.frame(results.presab)
rownames(results.presab)<-colnames(Y.presab)


preds.no.Cattle = computePredictedValues(models[[3]])
results.no.Cattle<-evaluateModelFit(hM=models[[3]], predY=preds.no.Cattle)
results.no.Cattle<-as.data.frame(results.no.Cattle)
rownames(results.no.Cattle)<-colnames(Y)

preds.presab.no.Cattle = computePredictedValues(models[[4]])
results.presab.no.Cattle<-evaluateModelFit(hM=models[[4]], predY=preds.presab.no.Cattle)
results.presab.no.Cattle<-as.data.frame(results.presab.no.Cattle)
rownames(results.presab.no.Cattle)<-colnames(Y.presab)

computeWAIC(models[[2]]) ] # 0.5964629
computeWAIC(models[[4]]) ] # 0.596367


#Variance Partitioning

VP.normal = computeVariancePartitioning(models[[1]])
vals.normal = VP.normal$vals
mycols = rainbow(nrow(VP.normal$vals))
plotVariancePartitioning(hM=models[[1]], VP=VP.normal, legend.text = TRUE, main = "Proportion of explained variance, spatial", las = 2, cex.names = 0.75)

VP.presab = computeVariancePartitioning(models[[2]])
vals.presab = VP.presab$vals
mycols = rainbow(nrow(VP.presab$vals))
plotVariancePartitioning(hM=models[[2]], VP=VP.presab, legend.text = TRUE, main = "Proportion of explained variance, spatial", las = 2, cex.names = 0.75)


VP.no.Cattle  = computeVariancePartitioning(models[[3]])
vals.no.Cattle  = VP.no.Cattle $vals
mycols = rainbow(nrow(VP.no.Cattle $vals))
plotVariancePartitioning(hM=models[[2]], VP=VP.no.Cattle, legend.text = TRUE, main = "Proportion of explained variance, spatial", las = 2, cex.names = 0.75)
