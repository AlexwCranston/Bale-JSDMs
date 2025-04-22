# ---- Load Libraries

library(Hmsc)
library(colorspace)
library(vioplot)

# ---- Set Directory

localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

# ---- Set Parameters

showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100
showRho = FALSE
showAlpha = TRUE


# ---- Evaluate Convergence


samples_list = c(5,250,250, 250)
thin_list = c(1,1,10, 100)
nst = length(thin_list)
nChains = 4

text.file.CA = file.path(resultDir,"/MCMC_convergence_CA.txt")
cat("MCMC Convergence statistics\n\n",file=text.file.CA,sep="")

text.file.presab = file.path(resultDir,"/MCMC_convergence_presab.txt")
cat("MCMC Convergence statistics\n\n",file=text.file.presab,sep="")

ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL  
ma.rho = NULL
na.rho = NULL

# ---- First CA models

Lst = 1
while(Lst <= nst){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename = file.path(modelDir,paste("models.CA_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),  
                                      ".Rdata",sep = ""))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  if(file.exists(filename)){
    load(filename)
    load(filename_unfitted)
    models.CA.list <- unfitted.models.list[[1]]
    cat(c("\n",filename,"\n\n"),file=text.file.CA,sep="",append=TRUE)
    nm = length(models_fitted)
    for(j in 1:nm){
      names(models_fitted)[j] <- names(models.CA.list)[j]
      mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      nr = models_fitted[[j]]$nr
      cat(c("\n",names(models_fitted)[j],"\n\n"),file=text.file.CA,sep="",append=TRUE)
      if(showBeta){
        psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\nbeta\n\n",file=text.file.CA,sep="",append=TRUE)
        cat(tmp[,1],file=text.file.CA,sep="\n",append=TRUE)
      }
      if(showGamma){
        psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\ngamma\n\n",file=text.file.CA,sep="",append=TRUE)
        cat(tmp[,1],file=text.file.CA,sep="\n",append=TRUE)
        
      }
      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat("\nrho\n\n",file=text.file.CA,sep="",append=TRUE)
        cat(psrf[1],file=text.file.CA,sep="\n",append=TRUE)
      }
      if(showOmega & nr>0){
        cat("\nomega\n\n",file=text.file.CA,sep="",append=TRUE)
        for(k in 1:nr){
          cat(c("\n",names(models_fitted[[j]]$ranLevels)[k],"\n\n"),file=text.file.CA,sep="",append=TRUE)
          tmp = mpost$Omega[[k]]
          z = dim(tmp[[1]])[2]
          if(z > maxOmega){
            sel = sample(1:z, size = maxOmega)
            for(i in 1:length(tmp)){
              tmp[[i]] = tmp[[i]][,sel]
            }
          }
          psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
          tmp = summary(psrf)
          cat(tmp[,1],file=text.file.CA,sep="\n",append=TRUE)
        }
      }
      if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models_fitted[[j]]$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n",file=text.file.CA,sep="\n",append=TRUE)
            cat(c("\n",names(models_fitted[[j]]$ranLevels)[k],"\n\n"),file=text.file.CA,sep="",append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=text.file.CA,sep="\n",append=TRUE)            
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

# ---- Next Presab Models

Lst = 1
while(Lst <= nst){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename = file.path(modelDir,paste("models.presab_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),  
                                      ".Rdata",sep = ""))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  if(file.exists(filename)){
    load(filename)
    load(filename_unfitted)
    models.CA.list <- unfitted.models.list[[1]]
    cat(c("\n",filename,"\n\n"),file=text.file.presab,sep="",append=TRUE)
    nm = length(models_fitted)
    for(j in 1:nm){
      names(models_fitted)[j] <- names(models.CA.list)[j]
      mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      nr = models_fitted[[j]]$nr
      cat(c("\n",names(models_fitted)[j],"\n\n"),file=text.file.presab,sep="",append=TRUE)
      if(showBeta){
        psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\nbeta\n\n",file=text.file.presab,sep="",append=TRUE)
        cat(tmp[,1],file=text.file.presab,sep="\n",append=TRUE)
      }
      if(showGamma){
        psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\ngamma\n\n",file=text.file.presab,sep="",append=TRUE)
        cat(tmp[,1],file=text.file.presab,sep="\n",append=TRUE)
      }
      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat("\nrho\n\n",file=text.file.presab,sep="",append=TRUE)
        cat(psrf[1],file=text.file.presab,sep="\n",append=TRUE)
      }
      if(showOmega & nr>0){
        cat("\nomega\n\n",file=text.file.presab,sep="",append=TRUE)
        for(k in 1:nr){
          cat(c("\n",names(models_fitted[[j]]$ranLevels)[k],"\n\n"),file=text.file.presab,sep="",append=TRUE)
          tmp = mpost$Omega[[k]]
          z = dim(tmp[[1]])[2]
          if(z > maxOmega){
            sel = sample(1:z, size = maxOmega)
            for(i in 1:length(tmp)){
              tmp[[i]] = tmp[[i]][,sel]
            }
          }
          psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
          tmp = summary(psrf)
          cat(tmp[,1],file=text.file.presab,sep="\n",append=TRUE)
        }
      }
      if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models_fitted[[j]]$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n",file=text.file.presab,sep="\n",append=TRUE)
            cat(c("\n",names(models_fitted[[j]]$ranLevels)[k],"\n\n"),file=text.file.presab,sep="",append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=text.file.presab,sep="\n",append=TRUE)            
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

# --- Make plots for CA


pdf(file= file.path(resultDir,"/MCMC_convergence_CA.pdf"))

if(showBeta){
  filename = file.path(modelDir,paste("models.CA_thin_100_samples_250_chains_4.Rdata"))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  load(filename)
  load(filename_unfitted)
  models.CA.list <- unfitted.models.list[[1]]
  names(models_fitted) <- names(models.CA.list)
  par(mfrow=c(2,1))
  nm = length(models_fitted)
  ma.beta = NULL
  na.beta = names(models_fitted)
  for(j in 1:nm){
    mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    if(is.null(ma.beta)){
      ma.beta = psrf[,1]
    } else {
      ma.beta = cbind(ma.beta,psrf[,1])
    }
  }
  par(las = 2, mar = c(5, 4, 4, 2))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}


if(showGamma){
  filename = file.path(modelDir,paste("models.CA_thin_100_samples_250_chains_4.Rdata"))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  load(filename)
  load(filename_unfitted)
  models.CA.list <- unfitted.models.list[[1]]
  names(models_fitted) <- names(models.CA.list)
  par(mfrow=c(2,1))
  nm = length(models_fitted)
  ma.gamma = NULL
  na.gamma = names(models_fitted)
  for(j in 1:nm){
    mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
    if(is.null(ma.gamma)){
      ma.gamma = psrf[,1]
    } else {
      ma.gamma = cbind(ma.gamma,psrf[,1])
    }
  }
  par(las = 2, mar = c(5, 4, 4, 2))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}


if(showOmega) {
  filename = file.path(modelDir,paste("models.CA_thin_100_samples_250_chains_4.Rdata"))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  load(filename)
  load(filename_unfitted)
  models.CA.list <- unfitted.models.list[[1]]
  names(models_fitted) <- names(models.CA.list)
  par(mfrow=c(2,1))
  nm = length(models_fitted)
  ma.omega = NULL
  na.omega = names(models_fitted)
  for(j in 1:nm){
    mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    tmp = mpost$Omega[[1]]
    psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
    if(is.null(ma.omega)){
      ma.omega = psrf[,1]
    } else {
      ma.omega = cbind(ma.omega,psrf[,1])
    }
  }
  par(las = 2, mar = c(5, 4, 4, 2))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}

dev.off()
message("PDF device closed, file ready to open.")


# --- Make plots for presab


pdf(file= file.path(resultDir,"/MCMC_convergence_presab.pdf"))

if(showBeta){
  filename = file.path(modelDir,paste("models.presab_thin_100_samples_250_chains_4.Rdata"))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  load(filename)
  load(filename_unfitted)
  models.presab.list <- unfitted.models.list[[2]]
  names(models_fitted) <- names(models.presab.list)
  par(mfrow=c(2,1))
  nm = length(models_fitted)
  ma.beta = NULL
  na.beta = names(models_fitted)
  for(j in 1:nm){
    mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    if(is.null(ma.beta)){
      ma.beta = psrf[,1]
    } else {
      ma.beta = cbind(ma.beta,psrf[,1])
    }
  }
  par(las = 2, mar = c(5, 4, 4, 2))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}


if(showGamma){
  filename = file.path(modelDir,paste("models.presab_thin_100_samples_250_chains_4.Rdata"))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  load(filename)
  load(filename_unfitted)
  models.presab.list <- unfitted.models.list[[2]]
  names(models_fitted) <- names(models.presab.list)
  par(mfrow=c(2,1))
  nm = length(models_fitted)
  ma.gamma = NULL
  na.gamma = names(models_fitted)
  for(j in 1:nm){
    mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
    if(is.null(ma.gamma)){
      ma.gamma = psrf[,1]
    } else {
      ma.gamma = cbind(ma.gamma,psrf[,1])
    }
  }
  par(las = 2, mar = c(5, 4, 4, 2))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}


if(showOmega) {
  filename = file.path(modelDir,paste("models.presab_thin_100_samples_250_chains_4.Rdata"))
  filename_unfitted = file.path(modelDir,paste("unfitted_models.Rdata"))
  load(filename)
  load(filename_unfitted)
  models.presab.list <- unfitted.models.list[[2]]
  names(models_fitted) <- names(models.presab.list)
  par(mfrow=c(2,1))
  nm = length(models_fitted)
  ma.omega = NULL
  na.omega = names(models_fitted)
  for(j in 1:nm){
    mpost = convertToCodaObject(models_fitted[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    tmp = mpost$Omega[[1]]
    psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
    if(is.null(ma.omega)){
      ma.omega = psrf[,1]
    } else {
      ma.omega = cbind(ma.omega,psrf[,1])
    }
  }
  par(las = 2, mar = c(5, 4, 4, 2))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}

dev.off()
message("PDF device closed, file ready to open.")
















