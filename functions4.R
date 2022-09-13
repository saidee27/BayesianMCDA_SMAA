
MCDA <- function(datamatrix, Nvector, Ntrt, Nendpts, weights, most, least, piece1, piece2, piece3, topvar, variance, a_in, b_in, contvars = contvars, binvars = binvars) {
  
  xi=datamatrix/Nvector
  
  pvf <- function(x, most, least, piece1, piece2 , piece3) {
    partialval = as.numeric(most > least) * (as.numeric(x <= piece1)*(((x-least)/(piece1-least)) * 0.25) + as.numeric(x > piece1)*as.numeric(x <= piece2)*(0.25+(((x-piece1)/(piece2-piece1)) * 0.25)) + as.numeric(x > piece2)*as.numeric(x <= piece3)*(0.5+(((x-piece2)/(piece3-piece2)) * 0.25)) + as.numeric(x > piece3)*(0.75 + ((x-piece3)/(most-piece3)) * 0.25)) + as.numeric(most < least) * (as.numeric(x >= piece3)*(((x-least)/(piece3-least)) * 0.25) + as.numeric(x >= piece2)*as.numeric(x < piece3)*(0.25+(((x-piece3)/(piece2-piece3)) * 0.25)) + as.numeric(x >= piece1)*as.numeric(x < piece2)*(0.5+(((x-piece2)/(piece1-piece2)) * 0.25)) + as.numeric(x < piece1)*(0.75 + ((x-piece1)/(most-piece1)) * 0.25))
    
    return(partialval)}
  values=pvf(t(xi), most, least, piece1, piece2, piece3)
  # Utility score
  us <- function (v, w) { return (v * w)}
  
  data1= us(values,weights)
  
  
  
  aflat=datamatrix+1
  bflat=Nvector-datamatrix+1
  
  alphaflat = datamatrix+0.001
  betaflat = Nvector + 0.001
  
  sigmaflat = Nvector/variance + 1/1000 
  
  muflat = ((datamatrix / variance) + (0 / 1000))/sigmaflat
  
  xiX = array(0, c(11, Ntrt, Nendpts))
  
  meanpost = seq(0,1,by=0.1)
  
  for (i in 1:Ntrt) {
    for(j in 1:Nendpts) {
      if (binvars[j] == 1) {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 1
          b_in[i,j] = 1
        }
        a=datamatrix[i,j]+a_in[i,j]
        b=Nvector[i,j]-datamatrix[i,j]+b_in[i,j]-a_in[i,j]
        xiX[,i,j]=meanpost * (a/(a+b)) + (1-meanpost)* (aflat[i,j]/(aflat[i,j]+bflat[i,j]))
      } else if (contvars[j] == 1) {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 0
          b_in[i,j] = 1000
        }
        sigma = Nvector[i,j]/variance[i,j] + 1/b_in[i,j] 
        mu = ((datamatrix[i,j] / variance[i,j]) + (a_in[i,j] / b_in[i,j]))/sigma
        xiX[,i,j]=meanpost * mu + (1-meanpost)* muflat[i,j]
      } else {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 0.001
          b_in[i,j] = 0.001
        }
        alpha = datamatrix[i,j]+a_in[i,j]
        beta = Nvector[i,j] + b_in[i,j]
        xiX[,i,j]=meanpost * (alpha/(beta)) + (1-meanpost)* (alphaflat[i,j]/(betaflat[i,j]))
      }
      
    }
  }
  
  valuesX = array(0, c(11, Nendpts, Ntrt))
  for (i in 1:11) {
    valuesX[i,,]=pvf(t(xiX[i,,]), most, least, piece1, piece2, piece3)
  }
  
  # Utility score
  usX <- function (v, w) { return (t(v) %*% w)}
  
  data2_t=matrix(NA,Ntrt,11)
  
  for (i in 1:11) {
    data2_t[,i] = usX(valuesX[i,,], weights)
  }
  
  data2 = t(data2_t)
  
  return(list(data1=data1,data2=data2))

}

pMCDA <- function(datamatrix, Nvector, Ntrt, Nendpts, a_in, b_in, variance, weights, most, least, piece1, piece2, piece3, contvars, binvars, wmix) {
  
  nsim=100000 # nb of simulations to obtain the posterior distributions
  
  aflat=datamatrix+1
  bflat=Nvector-datamatrix+1
  
  alphaflat = datamatrix+0.001
  betaflat = Nvector + 0.001

  sigmaflat = Nvector/variance + 1/1000 
  
  muflat = ((datamatrix / variance) + (0 / 1000))/sigmaflat
  
  samplepost = rbinom(nsim,1,wmix)
  
  # Distribution of the parameters
  xi = array(0, c(nsim, Ntrt, Nendpts))
  for (i in 1:Ntrt) {
    for(j in 1:Nendpts) {
      if (binvars[j] == 1) {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 1
          b_in[i,j] = 1
        }
        a=datamatrix[i,j]+a_in[i,j]
        b=Nvector[i,j]-datamatrix[i,j]+b_in[i,j]-a_in[i,j]
        xi[,i,j]=samplepost * rbeta(nsim, a, b) + (1-samplepost)* rbeta(nsim, aflat[i,j], bflat[i,j])
      } else if (contvars[j] == 1) {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 0
          b_in[i,j] = 1000
        }
        sigma = Nvector[i,j]/variance[i,j] + 1/b_in[i,j] 
        mu = ((datamatrix[i,j] / variance[i,j]) + (a_in[i,j] / b_in[i,j]))/sigma
        xi[,i,j]=samplepost * rnorm(nsim, mean = mu, sd = sqrt(sigma)) + (1-samplepost) * rnorm(nsim, mean = muflat[i,j], sd = sqrt(sigmaflat[i,j]))
      } else {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 0.001
          b_in[i,j] = 0.001
        }
        alpha = datamatrix[i,j]+a_in[i,j]
        beta = Nvector[i,j] + b_in[i,j]
        xi[,i,j]=samplepost * rgamma(nsim, shape = alpha, rate = beta) + (1-samplepost) * rgamma(nsim, shape = alphaflat[i,j], rate = betaflat[i,j])
      }
      
    }
  }
  
  
  pvf <- function(x, most, least, piece1, piece2 , piece3) {
    partialval = as.numeric(most > least) * (as.numeric(x <= piece1)*(((x-least)/(piece1-least)) * 0.25) + as.numeric(x > piece1)*as.numeric(x <= piece2)*(0.25+(((x-piece1)/(piece2-piece1)) * 0.25)) + as.numeric(x > piece2)*as.numeric(x <= piece3)*(0.5+(((x-piece2)/(piece3-piece2)) * 0.25)) + as.numeric(x > piece3)*(0.75 + ((x-piece3)/(most-piece3)) * 0.25)) + as.numeric(most < least) * (as.numeric(x >= piece3)*(((x-least)/(piece3-least)) * 0.25) + as.numeric(x >= piece2)*as.numeric(x < piece3)*(0.25+(((x-piece3)/(piece2-piece3)) * 0.25)) + as.numeric(x >= piece1)*as.numeric(x < piece2)*(0.5+(((x-piece2)/(piece1-piece2)) * 0.25)) + as.numeric(x < piece1)*(0.75 + ((x-piece1)/(most-piece1)) * 0.25))
    
    return(partialval)}
  values = array(0, c(nsim, Nendpts, Ntrt))
  for (i in 1:nsim) {
    values[i,,]=pvf(t(xi[i,,]), most, least, piece1, piece2, piece3)
  }
  
  rm(xi)
  
  # Utility score
  us <- function (v, w) { return (t(v) %*% w)}
  
  uss=matrix(NA,Ntrt,nsim)
  
  difference=vector(length=nsim)
  
  for (i in 1:nsim) {
    uss[,i] = us(values[i,,], weights)
  }
  
  return(t(uss))
  
}



SMAA <- function(datamatrix, Nvector, Ntrt, Nendpts, a_in, b_in, variance, weights, most, least, piece1, piece2, piece3, contvars, binvars, wmix) {
  
  nsim=100000 # nb of simulations to obtain the posterior distributions
  
  
  aflat=datamatrix+1
  bflat=Nvector-datamatrix+1
  
  alphaflat = datamatrix+0.001
  betaflat = Nvector + 0.001
  
  sigmaflat = Nvector/variance + 1/1000 
  
  muflat = ((datamatrix / variance) + (0 / 1000))/sigmaflat
  
  samplepost = rbinom(nsim,1,wmix)
  
  # Distribution of the parameters
  xi = array(0, c(nsim, Ntrt, Nendpts))
  for (i in 1:Ntrt) {
    for(j in 1:Nendpts) {
      if (binvars[j] == 1) {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 1
          b_in[i,j] = 1
        }
        a=datamatrix[i,j]+a_in[i,j]
        b=Nvector[i,j]-datamatrix[i,j]+b_in[i,j]-a_in[i,j]
        xi[,i,j]=samplepost * rbeta(nsim, a, b) + (1-samplepost)* rbeta(nsim, aflat[i,j], bflat[i,j])
      } else if (contvars[j] == 1) {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 0
          b_in[i,j] = 1000
        }
        sigma = Nvector[i,j]/variance[i,j] + 1/b_in[i,j] 
        mu = ((datamatrix[i,j] / variance[i,j]) + (a_in[i,j] / b_in[i,j]))/sigma
        xi[,i,j]=samplepost * rnorm(nsim, mean = mu, sd = sqrt(sigma)) + (1-samplepost) * rnorm(nsim, mean = muflat, sd = sqrt(sigmaflat))
      } else {
        if (is.na(a_in[i,j])) {
          a_in[i,j] = 0.001
          b_in[i,j] = 0.001
        }
        alpha = datamatrix[i,j]+a_in[i,j]
        beta = Nvector[i,j] + b_in[i,j]
        xi[,i,j]=samplepost * rgamma(nsim, shape = alpha, rate = beta) + (1-samplepost) * rgamma(nsim, shape = alphaflat[i,j], rate = betaflat[i,j])
      }
      
    }
  }
  
  
  pvf <- function(x, most, least, piece1, piece2 , piece3) {
    partialval = as.numeric(most > least) * (as.numeric(x <= piece1)*(((x-least)/(piece1-least)) * 0.25) + as.numeric(x > piece1)*as.numeric(x <= piece2)*(0.25+(((x-piece1)/(piece2-piece1)) * 0.25)) + as.numeric(x > piece2)*as.numeric(x <= piece3)*(0.5+(((x-piece2)/(piece3-piece2)) * 0.25)) + as.numeric(x > piece3)*(0.75 + ((x-piece3)/(most-piece3)) * 0.25)) + as.numeric(most < least) * (as.numeric(x >= piece3)*(((x-least)/(piece3-least)) * 0.25) + as.numeric(x >= piece2)*as.numeric(x < piece3)*(0.25+(((x-piece3)/(piece2-piece3)) * 0.25)) + as.numeric(x >= piece1)*as.numeric(x < piece2)*(0.5+(((x-piece2)/(piece1-piece2)) * 0.25)) + as.numeric(x < piece1)*(0.75 + ((x-piece1)/(most-piece1)) * 0.25))
    
    return(partialval)}
  values = array(0, c(nsim, Nendpts, Ntrt))
  for (i in 1:nsim) {
    values[i,,]=pvf(t(xi[i,,]), most, least, piece1, piece2, piece3)
  }
  
  rm(xi)
  
  c = mergeConstraints(lowerRatioConstraint(length(weights),1,2,100/(weights[[2]][2])),upperRatioConstraint(length(weights),1,2,100/(weights[[2]][1])))
  for(j in 3:length(weights)) { c = mergeConstraints(c,lowerRatioConstraint(length(weights),1,j,100/(weights[[j]][2])),upperRatioConstraint(length(weights),1,j,100/(weights[[j]][1])) )}
  c = mergeConstraints(c,simplexConstraints(length(weights)))
  w = hitandrun(c, n.samples = nsim)
  
  # Utility score
  us <- function (v, w) { return (t(v) %*% w)}
  
  uss=matrix(NA,Ntrt,nsim)
  
  difference=vector(length=nsim)
  
  for (i in 1:nsim) {
    uss[,i] = us(values[i,,], w[i,])
  }
  
  return(t(uss))
  
}
