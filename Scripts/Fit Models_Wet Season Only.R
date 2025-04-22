# ----Load in the libraries

library(tidyverse)
library(Hmsc)

# ---- Set Seed

set.seed(1)

# ---- Set Directory

localDir = "."
modelDir = file.path(localDir, "models")


# ---- Load in the unfitted models

load(file=file.path(modelDir,"unfitted_models.RData"))

# ---- Load in the data

models.CA.list<-unfitted.models.list[["CA.models"]]
models.presab.list<-unfitted.models.list[["presab.models"]]


nm = length(models.CA.list)

samples_list = c(5,250,250,250)
thin_list = c(1,1,10,100) 
nChains = 4
nParallel = nChains

## presab first


for (Lst in seq_along(samples_list)) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  models_fitted <- list()
  filename = file.path(modelDir,paste("models.presab_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),  
                                      ".Rdata",sep = ""))
  print(filename)
  start_time <- Sys.time()
  for (mi in 1:nm) {
    m = models.presab.list[[mi]]
    print(paste0("model = ", names(models.presab.list)[[mi]]))
    m = sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains,
                   nParallel = nParallel) 
    models_fitted[[mi]] = m
  }
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  print(paste0("Completed in ", round(elapsed, 2), " ", attr(elapsed, "units")))
  save(models_fitted,file=filename)
}

## CA next

for (Lst in seq_along(samples_list)) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  models_fitted <- list()
  filename = file.path(modelDir,paste("models.CA_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),  
                                      ".Rdata",sep = ""))
  print(filename)
  start_time <- Sys.time()
  for (mi in 1:nm) {
      m = models.CA.list[[mi]]
      print(paste0("model = ", names(models.CA.list)[[mi]]))
      m = sampleMcmc(m, samples = samples, thin=thin,
                     adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains,
                     nParallel = nParallel) 
      models_fitted[[mi]] = m
  }
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  print(paste0("Completed in ", round(elapsed, 2), " ", attr(elapsed, "units")))
  save(models_fitted,file=filename)
  }

