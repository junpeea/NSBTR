########################################################################################
########################################################################################
# Thesis : Bayesian time series regression with
# nonparametric modeling of autocorrelation
# Project name : NBTR_Func (Function & Libraries)

# Date : 201712xx
# Author : YB JUN
########################################################################################
########################################################################################
require(coda);require(MASS);require(mvtnorm);require(KernSmooth)
require(geoR);require(fields);require(spam);require(fftwtools)
require(matrixcalc);require(lattice);require(LatticeKrig);require(NORMT3)
require(RandomFields);require(clusterGeneration);require(splines);require(dplyr)
require(bayesGARCH);require(rugarch)

SampleSetSeed = function(n.train,n.test,seed){
  # seed   : seed number to control SampleSetSeed
  # n.train : number of train datasets
  # n.test  : number of test datasets
  set.seed(seed)
  seednum.vec1 = sample(n.train:(100*n.train),n.train)
  set.seed(seed)
  seednum.vec2 = sample((100*n.train+1):(100*n.train+10*n.test),n.test)
  
  seednum = c(seednum.vec1,seednum.vec2)
  return(seednum)
}

CompAR2SpectDen2 =function(a1,a2,sigma2,n){
  #compute AR2 spectral density  with w in (0,1/2)
  #for AR1, let a2=0
  #  eps_t -a1 eps_{t-1} -a2 eps_{t-2} = sqrt(sigma2)* z_t
  #plot(CompAR2SpectDen(0.5,-0.3,0.5^2,n))
  # here sigma2 is variance for the inovation term
  k <- seq(0,floor(n/2),1)
  w.j <- k/n
  
  res= (1-a1*cos(2*pi*w.j) -a2*cos(2*pi*2*w.j))^2
  res=res+(a1*sin(2*pi*w.j) +a2*sin(2*pi*2*w.j))^2
  return(sigma2/res)
  
}

CompARMASpectDen =function(a1,b1,sigma2,n){
  #compute AR2 spectral density  with w in (0,1/2)
  #eps_t -a1 eps_{t-1} = b1 sqrt(sigma2)* z_(t-1)+sqrt(sigma2)* z_t
  # plot(CompAR2SpectDen(0.5,-0.3,0.5^2,n))
  
  k <- seq(0,floor(n/2),1)
  w.j <- k/n
  
  res1= (1-a1*cos(2*pi*w.j))^2+(a1*sin(2*pi*w.j))^2
  res2= (1+b1*cos(2*pi*w.j))^2+(b1*sin(2*pi*w.j))^2
  return(sigma2*res2/res1)
  
}

CompLogPerioFunc = function(epsi){
  #sub-routine to get log of peridogram
  
  n <- length(epsi)
  tvec <- seq(1,n,1)
  k <- seq(0,floor(n/2),1)
  w.j <- k/n
  
  peri.temp = function(w.j)
  { 
    ((sum(epsi * cos(2* pi * tvec * w.j)))^2 + (sum(epsi * sin(2* pi * tvec * w.j)))^2)/n
  }
  
  peri <- log(sapply(w.j,peri.temp))
  return(peri)
}

  
CompQo=function(n){
  Q = matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      s = (i-1)*(j-1)*2*pi/n
      Q[i,j] = complex(real=cos(s),imaginary=sin(s))/sqrt(n)
    }
  }
  return(Q)
}

mydist = function(m,n){
  w.j <- seq(0,m,1)/n
  # if(n%%2==0){w.j <- seq(0,m-1,1)/n
  # }else{w.j <- seq(0,m-1,1)/n}
  dist0=as.matrix(dist(cbind(w.j,w.j)))
  return(dist0)
}


AR.mysd = function(truepar){
  # current version : AR(2) only!!!!!
  sig0 = NULL
  ar = truepar$ar
  sig0 = sqrt( (1+ar[2]) * ((1-ar[2])^2 - ar[1]^2) / (1-ar[2])  )
  # sig0 = sqrt( (1-ar[2]) / ((1+ar[2]) * ((1-ar[2])^2 - ar[1]^2)))
  return(sig0)
}

mysd = function(modelkey,truepar){
  # current version : AR(2) only!!!!!
  sig0 = NULL
  if(modelkey[2]==2){
    ar = truepar$ar
    sig0 = sqrt( (1+ar[2]) * ((1-ar[2])^2 - ar[1]^2) / (1-ar[2]) )
  }
  return(sig0)
}

Null_Initialize_Parameters = function(data00,seed){
  
  # known paramters (Reference: CARTER and KOHN (1997))
  kl.w0 = c(-4.63,-2.87,-1.44,-0.33,0.76)
  kl.wj = c(-2.20,-0.80,-0.55,-0.035,0.48)
  vl.w0 = c(8.75,1.95,0.88,0.45,0.41)
  vl.wj = c(1.93,1.01,0.69,0.60,0.29)
  pl.w0 = c(0.13,0.16,0.23,0.22,0.25)
  pl.wj = c(0.19,0.11,0.27,0.25,0.18)
  
  X <- cbind(data00$X1,data00$X2,data00$X3)
  Y <- data00$Y
  n <- nrow(data00); m <- floor(n/2)
  
  # initialize beta
  XX = as.matrix(t(X)%*%X)
  betahat <- chol2inv(chol(XX)) %*% (t(X) %*% Y)
  # perturbation
  seednum = seed + 10; set.seed(seednum)
  betahat <- betahat + rnorm(length(betahat),0,5e-01)
  
  # initial residual
  resi <- Y - X %*% betahat
  
  # initial sig2
  sig2 <- 1
  scaled.resi <- resi / sqrt(sig2)
  
  # initial phi
  phi <- CompLogPerioFunc(scaled.resi)
  # initial lambda
  w.j <- seq(0,m,1)/n
  lambda <- ksmooth(w.j,exp(phi),x.points=w.j,bandwidth=0.1)
  
  # initial theta
  theta  <- log(lambda$y)
  
  # initial label
  xi <- phi - theta
  label = vector(length=(m+1),mode="numeric")
  for(i in 1:(m+1)){
    if(i > 1 & i < (m+1)){kl = kl.wj; vl = vl.wj
    }else{                kl = kl.w0; vl = vl.w0}
    logdens = log(pl.wj) + dnorm(xi[i],mean=kl,sd=sqrt(vl),log=TRUE)
    logdens.scale = logdens - max(logdens)
    plabel = exp(logdens.scale)/sum(exp(logdens.scale))
    label[i] = sample(1:5,prob=plabel,size=1)
  }
  klabel = c(kl.w0[label[1]],kl.wj[label[2:m]],kl.w0[label[(m+1)]])
  vlabel = c(vl.w0[label[1]],vl.wj[label[2:m]],vl.w0[label[(m+1)]])
  
  # initialize ParaList
  ParaList = list()
  ParaList$beta       <- betahat
  ParaList$resi       <- resi
  ParaList$phi        <- phi
  ParaList$lambda     <- lambda$y
  ParaList$theta      <- theta
  ParaList$nu.theta   <- rep(0,(m+1))
  ParaList$rho.theta  <- 1
  ParaList$sig2       <- sig2
  ParaList$label      <- label
  ParaList$klabel     <- klabel
  ParaList$vlabel     <- vlabel
  ParaList$seed       <- seednum
  
  return(ParaList)
}

Initialize_Parameters = function(data00,seed){
  
  # known paramters (Reference: CARTER and KOHN (1997))
  kl.w0 = c(-4.63,-2.87,-1.44,-0.33,0.76)
  kl.wj = c(-2.20,-0.80,-0.55,-0.035,0.48)
  vl.w0 = c(8.75,1.95,0.88,0.45,0.41)
  vl.wj = c(1.93,1.01,0.69,0.60,0.29)
  pl.w0 = c(0.13,0.16,0.23,0.22,0.25)
  pl.wj = c(0.19,0.11,0.27,0.25,0.18)
  
  X <- cbind(data00$X1,data00$X2,data00$X3)
  Y <- data00$Y
  n <- nrow(data00); m <- floor(n/2)
  
  # initialize beta
  XX = as.matrix(t(X)%*%X)
  betahat <- chol2inv(chol(XX)) %*% (t(X) %*% Y)
  # perturbation
  seednum = seed + 10; set.seed(seednum)
  betahat <- betahat + rnorm(length(betahat),0,5e-01)
  
  # initial residual
  resi <- Y - X %*% betahat
  
  # initial eta
  nu.eta  = log( smooth.spline(c(1:n),resi^2,df=5)$y )
  nu.eta[is.na(nu.eta)] <- -10
  logsig2 = nu.eta
  # plot(nu.eta)
  # plot(log(truepar$sd[1:n.train]^2),ylim=c(-1.5,1.0))
  # logsig2 = log(truepar$sd[1:n]^2) # need to get rid of... later
    # plot(resi,ylim=c(-3,3),ylab="")
    # par(new=T);plot(log(truepar$sd[1:n]^2),type='l',col=2,ylab="",ylim=c(-3,3),lty=2)
    # par(new=T);plot(logsig2,type='l',col=4,ylab="",ylim=c(-3,3),lty=2)
  scaled.resi <- resi / sqrt(exp(logsig2))
  # initial delta
  d = 3; k = 10
  Bmatrix   = bs(seq(1,(n+k),by=1),knots=c(0.25,0.5,0.75)*(n+k),Boundary.knots=c(1,(n+k)),degree=d,intercept=TRUE)
  Bmatrix00 = Bmatrix[1:n,]
  # plot(seq(1,n,by=1),Bmatrix00[,1],xlim=c(1,n),type='n',xlab='x',ylab='')
  # for (i in 1:dim(Bmatrix00)[2]) lines(seq(1,n,by=1),Bmatrix00[,i])
  
  delta <- rep(0,dim(Bmatrix00)[2])
  
  # initial phi
  phi <- CompLogPerioFunc(scaled.resi)
  # initial lambda
  w.j <- seq(0,m,1)/n
  lambda <- ksmooth(w.j,exp(phi),x.points=w.j,bandwidth=0.1)
  
  # initial theta
  theta  <- log(lambda$y)
  
  # initial label
  xi <- phi - theta
  label = vector(length=(m+1),mode="numeric")
  for(i in 1:(m+1)){
    if(i > 1 & i < (m+1)){kl = kl.wj; vl = vl.wj
    }else{                kl = kl.w0; vl = vl.w0}
    logdens = log(pl.wj) + dnorm(xi[i],mean=kl,sd=sqrt(vl),log=TRUE)
    logdens.scale = logdens - max(logdens)
    plabel = exp(logdens.scale)/sum(exp(logdens.scale))
    label[i] = sample(1:5,prob=plabel,size=1)
  }
  klabel = c(kl.w0[label[1]],kl.wj[label[2:m]],kl.w0[label[(m+1)]])
  vlabel = c(vl.w0[label[1]],vl.wj[label[2:m]],vl.w0[label[(m+1)]])
  
  # initialize ParaList
  ParaList = list()
  ParaList$beta       <- betahat
  ParaList$resi       <- resi
  ParaList$phi        <- phi
  ParaList$lambda     <- lambda$y
  ParaList$theta      <- theta
  ParaList$nu.theta   <- rep(0,(m+1))
  ParaList$rho.theta  <- 1
  ParaList$Bmatrix00  <- Bmatrix00
  ParaList$delta      <- delta
  ParaList$eta        <- logsig2
  ParaList$nu.eta     <- nu.eta
  ParaList$rho.eta    <- 1
  ParaList$tau.eta    <- 1
  ParaList$label      <- label
  ParaList$klabel     <- klabel
  ParaList$vlabel     <- vlabel
  ParaList$seed       <- seednum
  
  return(ParaList)
}

Null_Update1_beta_phi = function(Parset_,data00,beta.mu0=0,beta.sig20=100,a=1,b=1){
  
  X <- cbind(data00$X1,data00$X2,data00$X3)
  Y <- data00$Y
  n <- nrow(data00); m <- ceiling(n/2)
  seednum <- Parset_$seed
  
  lambda_m <- Parset_$lambda
  if(n%%2==0){lambda_n <- c(lambda_m[seq(1,(m+1),1)],lambda_m[seq(m,2,by=-1)])
  }else{lambda_n <- c(lambda_m[seq(1,m,1)],lambda_m[seq(m,2,by=-1)])}
  
  sig2   <- Parset_$sig2 
  inv_covar_n <- diag(1/sqrt(sig2),n)%*%Qo%*%diag(1/lambda_n)%*%Qc%*%diag(1/sqrt(sig2),n)
  
  beta.sig2.star <- solve((t(X)%*%Re(inv_covar_n)%*%X + 1/beta.sig20))
  beta.mu.star   <- beta.sig2.star%*%(t(X)%*%Re(inv_covar_n)%*%Y + 1/beta.sig20*beta.mu0)
  seednum <- seednum + 10; set.seed(seednum)
  library(MASS);  newbeta <- as.matrix(mvrnorm(1,beta.mu.star,beta.sig2.star),length(beta.mu.star),1)
  
  newresi        <- Y - X%*%newbeta
  scaled.newresi <- newresi / sqrt(sig2)
  newphi  <- Re(CompLogPerioFunc(scaled.newresi))
  
  ## Update sig2 ##
  new.a = a+n/2; new.b = b + t(newresi)%*%Re(inv_covar_n)%*%(newresi)/2
  seednum <- seednum + 10; set.seed(seednum)
  newsig2 = 1/rgamma(1,new.a,new.b)
  
  newParset      <- Parset_
  newParset$beta <- newbeta
  newParset$resi <- newresi
  newParset$phi  <- newphi
  newParset$sig2 <- newsig2
  newParset$seed <- seednum
  return(newParset)
  
}


Update1_beta_phi = function(Parset_,data00,beta.mu0=0,beta.sig20=100,a=1,b=1){
  
  X <- cbind(data00$X1,data00$X2,data00$X3)
  Y <- data00$Y
  n <- nrow(data00); m <- ceiling(n/2)
  seednum <- Parset_$seed
  
  lambda_m <- Parset_$lambda
  if(n%%2==0){lambda_n <- c(lambda_m[seq(1,(m+1),1)],lambda_m[seq(m,2,by=-1)])
  }else{lambda_n <- c(lambda_m[seq(1,m,1)],lambda_m[seq(m,2,by=-1)])}
  
  
  eta  <- Parset_$eta
  sig2   <- exp(eta) 
  inv_covar_n <- diag(1/sqrt(sig2))%*%Qo%*%diag(1/lambda_n)%*%Qc%*%diag(1/sqrt(sig2))
  
  # beta.sig2.star <- solve(t(X)%*%inv_covar_n%*%X + 1/beta.sig20)
  # beta.sig2.star <- chol2inv(chol((t(X)%*%Re(inv_covar_n)%*%X + 1/beta.sig20)))
  beta.sig2.star <- solve((t(X)%*%Re(inv_covar_n)%*%X + 1/beta.sig20))
  beta.mu.star   <- beta.sig2.star%*%(t(X)%*%Re(inv_covar_n)%*%Y + 1/beta.sig20*beta.mu0)
  seednum <- seednum + 10; set.seed(seednum)
  library(MASS);  newbeta <- as.matrix(mvrnorm(1,beta.mu.star,beta.sig2.star),length(beta.mu.star),1)
  
  newresi        <- Y - X%*%newbeta
  scaled.newresi <- newresi / sqrt(sig2)
  newphi  <- Re(CompLogPerioFunc(scaled.newresi))
  
  ## Update kappa ##
  # new.a = a+n/2; new.b = b + t(newresi)%*%Re(inv_covar_n)%*%(newresi)/2
  # seednum <- seednum + 10; set.seed(seednum)
  # newkappa = sum(1/rgamma(n,new.a,new.b))
  newkappa = kappa
  
  newParset      <- Parset_
  newParset$beta <- newbeta
  newParset$resi <- newresi
  newParset$phi  <- newphi
  newParset$seed <- seednum
  return(newParset)
  
}
# Parset1 <- Update1_beta_phi(Parset0,data00)
# par(mfrow=c(1,2))
# plot(exp(Parset0$phi)/sum(exp(Parset0$phi)),ylim=c(0,0.1),type='o',xlab="freq",ylab="",main="spd");par(new=T);plot(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],mysd(modelkey,truepar)^2,n)/sum(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],mysd(modelkey,truepar)^2,n)),type='l',col=2,ylim=c(0,0.1),xlab="freq",ylab="",main="spd")
# plot(exp(Parset1$phi)/sum(exp(Parset0$phi)),ylim=c(0,0.1),type='o',xlab="freq",ylab="",main="spd");par(new=T);plot(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],mysd(modelkey,truepar)^2,n)/sum(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],mysd(modelkey,truepar)^2,n)),type='l',col=2,ylim=c(0,0.1),xlab="freq",ylab="",main="spd")

Update2_xi_mixture = function(Parset_,data00){
  
  # known paramters (Reference: CARTER and KOHN (1997))
  kl.w0 = c(-4.63,-2.87,-1.44,-0.33,0.76)
  kl.wj = c(-2.20,-0.80,-0.55,-0.035,0.48)
  vl.w0 = c(8.75,1.95,0.88,0.45,0.41)
  vl.wj = c(1.93,1.01,0.69,0.60,0.29)
  pl.w0 = c(0.13,0.16,0.23,0.22,0.25)
  pl.wj = c(0.19,0.11,0.27,0.25,0.18)
  
  n <- nrow(data00); m <- floor(n/2)
  phi     <- Parset_$phi
  theta   <- Parset_$theta
  seednum <- Parset_$seed
  
  xi = Re(phi - theta)
  phi.mix  = rep(0,(m+1))
  phi.mean = rep(0,(m+1))
  phi.var  = rep(0,(m+1))
  
  for(i in 1:(m+1)){
    if(i==1 | i==(m+1)){
      log.alpha0 = array(0,5)
      for(j in 1:5){
        seednum = seednum + 10; set.seed(seednum)
        log.alpha0[j] <- log(pl.w0[j]) + dnorm(xi[i],kl.w0[j],sqrt(vl.w0[j]),log=TRUE)
      }
      log.alpha = log.alpha0 - max(log.alpha0)
      pi.hat = exp(log.alpha)/sum(exp(log.alpha))
      phi.mix[i]  <- sample(c(1:5),size=1,prob=pi.hat)
      phi.mean[i] <- kl.w0[phi.mix[i]]
      phi.var[i]  <- vl.w0[phi.mix[i]]
    }else{
      log.alpha0 = array(0,5)
      for(j in 1:5){
        seednum = seednum + 10; set.seed(seednum)
        log.alpha0[j] <- log(pl.wj[j]) + dnorm(xi[i],kl.wj[j],sqrt(vl.wj[j]),log=TRUE)
      }
      log.alpha = log.alpha0 - max(log.alpha0)
      pi.hat = exp(log.alpha)/sum(exp(log.alpha))
      phi.mix[i]  <- sample(c(1:5),size=1,prob=pi.hat)
      phi.mean[i] <- kl.wj[phi.mix[i]]
      phi.var[i]  <- vl.wj[phi.mix[i]]
    }
  }
  newParset <- Parset_
  newParset$label  <- phi.mix
  newParset$klabel <- phi.mean
  newParset$vlabel <- phi.var
  newParset$seed   <- seednum
  return(newParset)
  
}
# Parset2 <- Update2_xi_mixture(Parset1,data00)


Update3_theta = function(Parset_,data00){
  
  n <- nrow(data00); m <- floor(n/2)
  rho     <- Parset_$rho.theta
  phi     <- Parset_$phi
  label   <- Parset_$label
  klabel  <- Parset_$klabel
  vlabel  <- Parset_$vlabel
  seednum <- Parset_$seed
  
  nu.mm    = Parset_$nu.theta
  tau.mm   = exp(-rho*mydist(m,n))
  tau.star = chol2inv(chol((solve(tau.mm)+diag(1/vlabel))))
  nu.star  = tau.star%*%diag(1/vlabel)%*%(phi-klabel-nu.mm)+nu.mm
  seednum <- seednum + 10; set.seed(seednum)
  library(MASS); newtheta = mvrnorm(1,mu=nu.star,Sigma=tau.star)
  
  # rescaling process
  newlambda = ConvThetaToLambdaFunc(n,newtheta)
  newtheta  = ConvLambdaToThetaFunc(m,newlambda)

  newParset        <- Parset_
  newParset$lambda <- newlambda[1:(m+1)]
  newParset$theta  <- newtheta
  newParset$seed   <- seednum
  return(newParset)
  
}

# Parset3 <- Update3_theta(Parset2,data00)
# par(mfrow=c(1,2))
# plot(Parset2$lambda,ylim=c(0,2.5),type='o',xlab="freq",ylab="",main="spd");par(new=T);plot(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],mysd(modelkey,truepar)^2,n),ylim=c(0,2.5),type='l',col=2,xlab="freq",ylab="",main="spd")
# plot(Parset3$lambda,ylim=c(0,2.5),type='o',xlab="freq",ylab="",main="spd");par(new=T);plot(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],mysd(modelkey,truepar)^2,n),ylim=c(0,2.5),type='l',col=2,xlab="freq",ylab="",main="spd")

Update4_theta_prior  = function(Parset_,data00,rho.len = 100){
  
  n <- nrow(data00); m <- floor(n/2)
  theta   <- Parset_$theta
  nu      <- Parset_$nu.theta
  rho     <- Parset_$rho.theta
  seednum <- Parset_$seed
  
  # Update nu.theta -> pass!
  
  # Update rho.theta
  rho.vec = exp(seq(log(1),log(10),length.out=rho.len))
  logprob0 = vector(len=rho.len,mode="numeric")
  for(k in 1:rho.len){
    tau.mm2 = exp(-rho.vec[k]*mydist(m,n))
    inv.tau.mm2 = chol2inv(chol(tau.mm2))
    # inv.tau.nn2 = solve(tau.nn2)
    logprob0[k] = -0.5*t(Re(theta))%*%inv.tau.mm2%*%Re(theta)
    logprob0[k] = logprob0[k] + 0.5*determinant(inv.tau.mm2)$mod
  }
  logprob = logprob0-max(logprob0)
  prob.v   = exp(logprob)/sum(exp(logprob))
  cumprob  = cumsum(prob.v[1:rho.len])
  seednum <- seednum + 10; set.seed(seednum)
  uni <- runif(1)
  new.rho  = rho.vec[which(cumprob>uni)][1]
  if(is.na(new.rho)) new.rho <- rho
  
  newParset <- Parset_
  newParset$rho.theta  <- new.rho
  newParset$seed       <- seednum
  return(newParset)
  
}
# Parset4 <- Update4_theta_prior(Parset3,data00)

Update5_eta  = function(Parset_,data00){
  
  n <- nrow(data00); m <- ceiling(n/2)
  betahat  <- Parset_$beta
  resi     <- Parset_$resi
  Bmatrix00<- Parset_$Bmatrix00
  delta    <- Parset_$delta
  eta      <- Parset_$eta
  nu.eta   <- Parset_$nu.eta
  seednum  <- Parset_$seed
  
  ######################################################################### loglikelihood
  
  loglikelihood = function(model, data00){
    X = cbind(data00$X1,data00$X2,data00$X3)
    Y = data00$Y
    eta = Bmatrix00 %*% model %>% as.vector
    d = dmvnorm(t(Y-X%*%betahat),sigma=diag(exp(eta)),log=TRUE)
    d
  }

  ############################################################################# generate
  generate = function(x, sigma){
    w = ceiling(runif(1) * length(x))
    x[w] = x[w] + rnorm(1, 0, sigma[w])
    return(x)
  }
  
  ########################################################################### logproposal
  # Proposal distribution is defined as multivariate normal, with mean
  # zero and standard deviations sigma:
  logproposal = function(x1, x2, sigma){
    -0.5 * sum(((x1) - (x2))^2/(sigma+1e-12)^2)
  }
  
  ############################################################################# logprior
  require(mvtnorm)
  logprior = function(m){
    dmvnorm(m,mean=rep(0,length(m)),sigma=1*t(Bmatrix00)%*%Bmatrix00,log=TRUE)
  }
  
  ############################################################################# MCMC (M-H)
  require(TDD)
  m1 = delta
  sigma = rep(1,length(m1))
  x = Metropolis(loglikelihood, sigma, m1, niter = 1000, gen = generate,
                 logproposal = logproposal, logprior = logprior, burn = 0, save_int = 1,
                 data = data00)

  newdelta = colMeans(x$m)
  neweta  = Bmatrix00 %*% newdelta
  
  newParset      <- Parset_
  newParset$delta<- newdelta
  newParset$eta  <- as.vector(neweta)
  newParset$seed <- seednum
  return(newParset)
}
# Parset5 <- Update5_eta(Parset0,data00)
# par(mfrow=c(1,2))
# plot(log(exp(Parset0$eta)/sum(exp(Parset0$eta))),type='o',ylim=c(-7,-4),xlab="freq",ylab="",main="");par(new=T);plot(log(truepar$sd[1:n.train]^2/sum(truepar$sd[1:n.train]^2)),type='l',col=2,ylim=c(-7,-4),lty=2,xlab="freq",ylab="",main="")
# plot(log(exp(Parset5$eta)/sum(exp(Parset5$eta))),type='o',ylim=c(-7,-4),xlab="freq",ylab="",main="");par(new=T);plot(log(truepar$sd[1:n.train]^2/sum(truepar$sd[1:n.train]^2)),type='l',col=2,ylim=c(-7,-4),lty=2,xlab="freq",ylab="",main="")

Update6_eta_prior  = function(Parset_,data00,a=0,b=100,c=1,d=1e-02){

  n <- nrow(data00); m <- ceiling(n/2)
  resi    <- Parset_$resi
  kappa   <- Parset_$kappa
  eta     <- Parset_$eta
  nu      <- Parset_$nu.eta
  rho     <- Parset_$rho.eta
  tau     <- Parset_$tau.eta
  seednum <- Parset_$seed
  
  m.nu <- rep(a,n)           # mean hyperparameter for nu
  v.nu <- b*diag(rep(1,n))   # variance hyperparameter for nu
  dist0 = mydist(n-1,n)
  Sigma  <- (1/tau)*exp(-rho*dist0)
  
  ## Update nu.eta ##
  new.v <- chol2inv(chol(solve(v.nu)+solve(Sigma)))
  new.m <- new.v %*% (solve(v.nu)%*%m.nu + solve(Sigma)%*%eta)
  seednum <- seednum + 10; set.seed(seednum)
  new.nu.eta <- as.vector(rmvnorm(1,mean=new.m,sigma=new.v))
  
  ## Update rho.lsig2 ##
  rho.len = 100
  rho.vec = exp(seq(log(1),log(100),length.out=rho.len))
  logprob0 = vector(len=rho.len,mode="numeric")
  for(k in 1:rho.len){
    tau.nn2 = (1/tau)*exp(-rho.vec[k]*dist0)
    inv.tau.nn2 = chol2inv(chol(tau.nn2))
    # inv.tau.nn2 = solve(tau.nn2)
    logprob0[k] = -0.5*t(Re(eta-new.nu.eta))%*%inv.tau.nn2%*%Re(eta-new.nu.eta)
    logprob0[k] = logprob0[k] + 0.5*determinant(inv.tau.nn2)$mod
  }

  logprob = logprob0-max(logprob0)
  prob.v   = exp(logprob)/sum(exp(logprob))
  cumprob  = cumsum(prob.v[1:rho.len])
  seednum <- seednum + 10; set.seed(seednum)
  uni <- runif(1)
  new.rho  = rho.vec[which(cumprob>uni)][1]
  if(is.na(new.rho)) new.rho <- rho

  ## Update tau.lsig2 ##
  Sigma   <- (1/tau)*exp(-new.rho*dist0); inv_Sigma = chol2inv(chol(Sigma))
  new.c = c+n/2; new.d = d + t(eta-new.nu.eta)%*%inv_Sigma%*%(eta-new.nu.eta)/2
  seednum <- seednum + 10; set.seed(seednum)
  new.tau = rgamma(1,new.c,new.d)
  
  newParset <- Parset_
  newParset$nu.eta  <- new.nu.eta
  newParset$rho.eta <- new.rho
  newParset$tau.eta <- new.tau
  newParset$seed       <- seednum
  return(newParset)
  
}


ConvLambdaToThetaFunc=function(m,lambda){
  #function to convert lambda to theta
  #lambda: spectral density 1:n
  #theta=log(lambda) 1:m=floor(n/2)
  
  theta=as.matrix(log(lambda[1:(m+1)]))
  
  return(theta)
}



ConvThetaToLambdaFunc =function(n,theta){
  #function to convert theta to lambda
  # lambda: spectral density
  # theta=log(lambda)
  #library(Matrix)
  
  
  #compute exp(theta) for 1:n
  e.theta = matrix(0,n,1)
  
  if(n%%2 ==0){
    for(i in 1: (n/2+1))
    {
      e.theta[i] = exp(theta[i])
    }
    for(i in 1:(n/2-1))
      
    {
      e.theta[(n/2+1)+i] = exp(rev(theta)[i+1])
    }
    
    
  }
  else{
    
    for(i in 1: (floor(n/2)+1))
    {
      e.theta[i] = exp(theta[i])
    }
    
    for(i in 1:(n/2))
      
    {
      e.theta[(n/2+1)+i] = exp(rev(theta)[i+1])
    }
    
    
    
  }
  #try to rescale so that sum to n(need to check this)
  #\int_0^1 lambda =1 for the spectral density 
  #when it is for correlation function
  #discrete version is 1/n \sum lambda(w_j) =1
  #\sum lambda(w_j)=n
  
  #rescaled version
  return(n*e.theta/sum(e.theta))
  
  #unscaled version #do not converge
  #return(e.theta)
  
}

regression_matrix  <- function(data,p,constant){
  nrow <- as.numeric(dim(data)[1])
  nvar <- as.numeric(dim(data)[2])
  
  Y1 <- as.matrix(data, ncol = nvar)
  X <- embed(Y1, p+1)
  X <- X[,(nvar+1):ncol(X)]
  if(constant == TRUE){
    X <-cbind(rep(1,(nrow-p)),X)
  }
  Y = matrix(Y1[(p+1):nrow(Y1),])
  nvar2 = ncol(X)
  return = list(Y=Y,X=X,nvar2=nvar2,nrow=nrow) 
}

regression_matrix2  <- function(data,covariates,p,constant){
  nrow <- as.numeric(dim(data)[1])
  nvar <- as.numeric(dim(data)[2])
  
  Y1 <- as.matrix(data, ncol = nvar)
  X <- embed(Y1, p+1)
  X <- X[,(nvar+1):ncol(X)]
  if(constant == TRUE){
    X <-cbind(rep(1,(nrow-p)),X)
  }
  Y = matrix(Y1[(p+1):nrow(Y1),])
  X <- cbind(X,covariates[(p+1):nrow(Y1),])
  nvar2 = ncol(X)
  return = list(Y=Y,X=X,nvar2=nvar2,nrow=nrow) 
}

ar_companion_matrix <- function(beta){
  #check if beta is a matrix
  if (is.matrix(beta) == FALSE){
    stop('error: beta needs to be a matrix')
  }
  # dont include constant
  k = nrow(beta) - 1
  FF <- matrix(0, nrow = k, ncol = k)
  
  #insert identity matrix
  FF[2:k, 1:(k-1)] <- diag(1, nrow = k-1, ncol = k-1)
  
  temp <- t(beta[2:(k+1), 1:1])
  #state space companion form
  #Insert coeffcients along top row
  FF[1:1,1:k] <- temp
  return(FF)
}

gibbs_sampler <- function(X,Y,B0,sigma0,sigma2,theta0,D0,reps,out,out1){
  
  for(i in 1:reps){
    if (i %% 1000 == 0){
      print(sprintf("Interation: %d", i))
    }
    M = solve(solve(sigma0) + as.numeric(1/sigma2) * t(X) %*% X) %*%
      (solve(sigma0) %*% B0 + as.numeric(1/sigma2) * t(X) %*% Y)
    
    V = solve(solve(sigma0) + as.numeric(1/sigma2) * t(X) %*% X)
    
    chck = -1; chck.iter = 0
    while(chck < 0 & chck.iter < 10){   # check for stability
      
      # B <- M + t(rnorm(p+1+1) %*% chol(V))
      B <- M + t(rnorm(p+1+2) %*% chol(V))
      
      # Check : not stationary for 3 lags
      b = ar_companion_matrix(B)
      ee <- max(sapply(eigen(b)$values,abs))
      if( ee<=1){
        chck=1
      }
      chck.iter <- chck.iter + 1
    }
    # compute residuals
    resids <- Y- X%*%B
    # D1 = D0 + t(resids) %*% resids
    D1 = D0 + t(resids) %*% resids / length(resids)
    
    # keeps samples after burn period
    out[i,] <- t(matrix(c(t(B),sigma2)))
    
    
    #draw from Inverse Gamma
    z0 = rnorm(T0,1); #z0 = rnorm(T1,1)
    z0z0 = t(z0) %*% z0
    sigma2 = D1/z0z0
    
    # keeps samples after burn period
    out[i,] <- t(matrix(c(t(B),sigma2)))
    
    # compute 2 year forecasts
    yhat = rep(0,horizon)
    end = as.numeric(length(Y))
    yhat[1:2] = Y[(end-1):end,]
    cfactor = sqrt(sigma2)
    X_mat = c(1,rep(0,p))
    for(m in (p+1):horizon){
      for (lag in 1:p){
        #create X matrix with p lags
        X_mat[(lag+1)] = yhat[m-lag]
      }
      # Use X matrix to forecast yhat
      # yhat[m] = X_mat %*% B + rnorm(1) * cfactor
      # yhat[m] = c(X_mat,as.numeric(data0[(n-2+m),2:2])) %*% B + rnorm(1) * cfactor
      yhat[m] = c(X_mat,as.numeric(data0[(n-2+m),2:3])) %*% B + rnorm(1) * cfactor
    }
    
    out1[i,] <- yhat
  }
  return = list(out,out1)
}

MyBayesARCH = function(X,Y,n.pred,reps){
  n.dat = length(Y)
  results1 = list()
  results1[[1]] = matrix(NA,reps,(3+n.dat))
  results1[[2]] = matrix(NA,reps,n.pred)
  y = data00$Y
  e <- as.vector(y - X%*%c(1,1,1))
  set.seed(1)
  for(iter in 1:reps){
    if(iter %% 100 == 0){
      print(sprintf("Interation: %d", iter))
      cat("Beta", colMeans(results1[[1]][,1:3],na.rm=T),"\n")
    }
    if(iter%%100==1){
      temp = tryCatch(capture.output(MCMC <- bayesGARCH(e, control = list(n.chain = 3, l.chain = 100))), error = function(e){1})
      adj.iter <- 0
      while(class(temp) == "numeric"){
        adj.iter = adj.iter + 1
        temp = tryCatch(capture.output(MCMC <- bayesGARCH(e, control = list(n.chain = 3, l.chain = max(0,(100-10*adj.iter))))), error = function(e){1})
        if(adj.iter > 10) break
      }
    }
    
    invisible(capture.output(smpl <- formSmpl(MCMC, l.bi = 0)));  temp = summary(smpl)
    spec = ugarchspec(fixed.pars = list("omega"=temp$statistics[1,1],"alpha1"=temp$statistics[2,1],"beta1"=temp$statistics[3,1]))
    temp = try(fit <- ugarchfit(data = e, spec = spec, solver="hybrid"))
    if(is.null(fit@fit$fitted.values)==TRUE){
      # temp = try(fit <- ugarchfit(data = e, spec = spec, solver="lbfgs"))
      temp = try(fit <- ugarchfit(data = e, spec = ugarchspec(), solver="lbfgs"))
    }
    
    ytilde <- y - fit@fit$fitted.values
    fit0 = lm(ytilde~-1+X[,1]+X[,2]+X[,3])
    e <- as.numeric(fit0$residuals)
    
    par(mfrow=c(1,1))
    plot(fit@fit$fitted.values,type='o',main=iter,ylim=c(-3,3))
    
    results1[[1]][iter,1:3] <- as.numeric(fit0$coefficients)
    results1[[1]][iter,(3+1):(3+n)] <- sigma(fit)^2
    pred = ugarchforecast(fit,n.ahead=n.pred)
    results1[[2]][iter,] <- Xtilde[(n+1):(n+n.pred),]%*%fit0$coefficients + pred@forecast$seriesFor
    
  }
  return(results1)
}


