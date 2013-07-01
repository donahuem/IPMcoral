x<-matrix(ncol=1)
index=c(0.001,0.01,0.1,1)
for (i in 1:length(index)){
  dummy.fit$coefficients[1]=index[[i]]
  fo@fitFec[[1]]<-dummy.fit
  delta=1e-04  
  nfec <- 0
  fec.coeff.names <- c()
  for (i in 1:length(fo@fitFec)) {
    nfec <- nfec + length(fo@fitFec[[i]]$coefficients)
    fec.coeff.names <- c(fec.coeff.names, paste("reprod", i, names(fo@fitFec[[i]]$coefficients)))
  }
  nclonal <- 0
  clonal.coeff.names <- c()
  for (i in 1:length(co@fitFec)) {
    nclonal <- nclonal + length(co@fitFec[[i]]$coefficients)
    clonal.coeff.names <- c(clonal.coeff.names, paste("clonal",i, names(co@fitFec[[i]]$coefficients)))
  }

  npar <- length(gr1@fit$coefficients) + 1 + length(sv1@fit$coefficients) + 
    length(fo@offspringRel$coefficients) + 1 + (sum(!is.na(fo@fecConstants))) + 
    nfec + length(co@offspringRel$coefficients) + 1 + (sum(!is.na(co@fecConstants)))+nclonal
  elam <- rep(0, npar)
  nmes <- c(paste("grow", names(gr1@fit$coefficients)), "sd growth", 
            paste("surv", names(sv1@fit$coefficients)), paste("offspring rel", 
            names(fo@offspringRel$coefficients)), "sd offspring", 
            fec.coeff.names, paste("offspring rel",names(co@offspringRel$coefficients)), "sd offspring", 
            clonal.coeff.names)
  names(elam) <- nmes[1:npar]
  slam <- elam
  ####################No perturbation######################
  Pmatrix <- createIPMPmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, growObj = gr1, 
                              survObj = sv1)
  Cmatrix <- createIPMCmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, clonalObj = co)
  Fmatrix <- createIPMFmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, fecObj = fo)
  IPM <- Pmatrix + Fmatrix + Cmatrix
  lambda1 <- Re(eigen(IPM)$value[1])
  #########################Growthperturbation#######################
  for (param.test in 1:length(gr1@fit$coefficients)) {
    gr1@fit$coefficients[param.test] <- gr1@fit$coefficients[param.test] * 
      (1 + delta)
    Pmatrix <- createIPMPmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, growObj = gr1, 
                                survObj = sv1)
    Cmatrix <- createIPMCmatrix(nBigMatrix = 100, minSize = minsize, 
                                maxSize = maxsize, clonalObj = co)
    
    Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, fecObj = fo)
    IPM <- Pmatrix + Fmatrix + Cmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    gr1@fit$coefficients[param.test] <- gr1@fit$coefficients[param.test]/(1 + 
                                                                            delta)
    slam[param.test] <- (lambda2 - lambda1)/(gr1@fit$coefficients[param.test] * 
                                               delta)
    elam[param.test] <- (lambda2 - lambda1)/(lambda1 * delta)
  }
  slamparam.test <- param.test + 1
  #############################Growth SD perturbation##########################
  sd.store <- gr1@sd
  gr1@sd <- gr1@sd * (1 + delta)
  Pmatrix <- createIPMPmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, growObj = gr1, survObj = sv1)
  Cmatrix <- createIPMCmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, clonalObj = co)
  Fmatrix <- createIPMFmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, fecObj = fo)
  IPM <- Pmatrix + Fmatrix + Cmatrix
  lambda2 <- Re(eigen(IPM)$value[1])
  gr1@sd <- sd.store
  slam[param.test] <- (lambda2 - lambda1)/(gr1@sd * delta)
  elam[param.test] <- (lambda2 - lambda1)/(lambda1 * delta)
  count <- param.test
  #####################Survival Perturbation########################
  for (param.test in 1:length(sv1@fit$coefficients)) {
    sv1@fit$coefficients[param.test] <- sv1@fit$coefficients[param.test] * 
      (1 + delta)
    Pmatrix <- createIPMPmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, growObj = gr1, 
                                survObj = sv1)
    Cmatrix <- createIPMCmatrix(nBigMatrix = 100, minSize = minsize, 
                                maxSize = maxsize, clonalObj = co)
    
    Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, fecObj = fo)
    IPM <- Pmatrix + Fmatrix + Cmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    sv1@fit$coefficients[param.test] <- sv1@fit$coefficients[param.test]/(1 + 
                                                                            delta)
    slam[param.test + count] <- (lambda2 - lambda1)/(sv1@fit$coefficients[param.test] * 
                                                       delta)
    elam[param.test + count] <- (lambda2 - lambda1)/(lambda1 * 
                                                       delta)
  }
  count <- count + param.test
  
  ################################Fecundity Rel#############################
  for (param.test in 1:length(fo@offspringRel$coefficients)) {
    fo@offspringRel$coefficients[param.test] <- fo@offspringRel$coefficients[param.test] * 
      (1 + delta)
    Pmatrix <- createIPMPmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, growObj = gr1, 
                                survObj = sv1)
    Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, fecObj = fo)
    
    Cmatrix <- createIPMCmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, clonalObj = co)
    IPM <- Pmatrix + Fmatrix+ Cmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    fo@offspringRel$coefficients[param.test] <- fo@offspringRel$coefficients[param.test]/(1 + 
                                                                                            delta)
    slam[param.test + count] <- (lambda2 - lambda1)/(fo@offspringRel$coefficients[param.test] * 
                                                       delta)
    elam[param.test + count] <- (lambda2 - lambda1)/(lambda1 * 
                                                       delta)
  }
  count <- count + param.test
  #################################FecunditySDOffspringSize######################
  fo@sdOffspringSize <- fo@sdOffspringSize * (1 + delta)
  Pmatrix <- createIPMPmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, growObj = gr1,
                              survObj = sv1)
  Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                              minSize = minsize, maxSize = maxsize, fecObj = fo)
  Cmatrix <- createIPMCmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, clonalObj = co)
  IPM <- Pmatrix + Fmatrix + Cmatrix
  lambda2 <- Re(eigen(IPM)$value[1])
  fo@sdOffspringSize <- fo@sdOffspringSize/(1 + delta)
  slam[count + 1] <- (lambda2 - lambda1)/(fo@sdOffspringSize * 
                                            delta)
  elam[count + 1] <- (lambda2 - lambda1)/(lambda1 * delta)
  chs <- which(!is.na(fo@fecConstants), arr.ind = TRUE)[, 
                                                        2]
  count <- count + 1
  ###############################FecundityFit############################
  for (i in 1:length(fo@fitFec)) {
    for (param.test in 1:length(fo@fitFec[[i]]$coefficients)) {
      fo@fitFec[[i]]$coefficients[param.test] <- fo@fitFec[[i]]$coefficients[param.test] * 
        (1 + delta)
      Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                                  minSize = minsize, maxSize = maxsize, fecObj = fo)         
      Cmatrix1 <- createIPMCmatrix(nBigMatrix = 100, 
                                   minSize = minsize, maxSize = maxsize, clonalObj = co)
      IPM <- Pmatrix + Fmatrix + Cmatrix1
      lambda2 <- Re(eigen(IPM)$value[1])
      fo@fitFec[[i]]$coefficients[param.test] <- fo@fitFec[[i]]$coefficients[param.test]/(1 + 
                                                                                            delta)
      slam[param.test + count] <- (lambda2 - lambda1)/(fo@fitFec[[i]]$coefficients[param.test] * 
                                                         delta)
      elam[param.test + count] <- (lambda2 - lambda1)/(lambda1 * 
                                                         delta)
    }
    count <- count + param.test
  }
  #############################Clonal perturbation Rel######################
  for (param.test in 1:length(co@offspringRel$coefficients)) {
    co@offspringRel$coefficients[param.test] <- co@offspringRel$coefficients[param.test] * 
      (1 + delta)
    Pmatrix <- createIPMPmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, growObj = gr1, 
                                survObj = sv1)
    Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, fecObj = fo)
    
    Cmatrix <- createIPMCmatrix(nBigMatrix = 100, 
                                minSize = minsize, maxSize = maxsize, clonalObj = co)
    IPM <- Pmatrix + Fmatrix+ Cmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    co@offspringRel$coefficients[param.test] <- co@offspringRel$coefficients[param.test]/(1 + 
                                                                                            delta)
    slam[param.test + count] <- (lambda2 - lambda1)/(co@offspringRel$coefficients[param.test] * 
                                                       delta)
    elam[param.test + count] <- (lambda2 - lambda1)/(lambda1 * 
                                                       delta)
  }
  count <- count + param.test
  ######################Clonal PerturbationSDOffspringSize############################
  co@sdOffspringSize <- co@sdOffspringSize * (1 + delta)
  Pmatrix <- createIPMPmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, growObj = gr1,
                              survObj = sv1)
  Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                              minSize = minsize, maxSize = maxsize, fecObj = fo)
  Cmatrix <- createIPMCmatrix(nBigMatrix = 100, minSize = minsize, 
                              maxSize = maxsize, clonalObj = co)
  IPM <- Pmatrix + Fmatrix + Cmatrix
  lambda2 <- Re(eigen(IPM)$value[1])
  co@sdOffspringSize <- co@sdOffspringSize/(1 + delta)
  slam[count + 1] <- (lambda2 - lambda1)/(co@sdOffspringSize * 
                                            delta)
  elam[count + 1] <- (lambda2 - lambda1)/(lambda1 * delta)
  chs <- which(!is.na(co@fecConstants), arr.ind = TRUE)[, 
                                                        2]
  count <- count + 1
  ################################ClonalFit#######################################
  for (i in 1:length(co@fitFec)) {
    for (param.test in 1:length(co@fitFec[[i]]$coefficients)) {
      co@fitFec[[i]]$coefficients[param.test] <- co@fitFec[[i]]$coefficients[param.test] * 
        (1 + delta)
      Fmatrix <- createIPMFmatrix(nBigMatrix = 100, 
                                  minSize = minsize, maxSize = maxsize, fecObj = fo)         
      Cmatrix1 <- createIPMCmatrix(nBigMatrix = 100, 
                                   minSize = minsize, maxSize = maxsize, clonalObj = co)
      IPM <- Pmatrix + Fmatrix + Cmatrix1
      lambda2 <- Re(eigen(IPM)$value[1])
      co@fitFec[[i]]$coefficients[param.test] <- co@fitFec[[i]]$coefficients[param.test]/(1 + 
                                                                                            delta)
      slam[param.test + count] <- (lambda2 - lambda1)/(co@fitFec[[i]]$coefficients[param.test] * 
                                                         delta)
      elam[param.test + count] <- (lambda2 - lambda1)/(lambda1 * 
                                                         delta)
    }
    count <- count + param.test
    
  }
  x<-rbind(x,as.matrix(elam))
}