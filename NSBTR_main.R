rm(list = ls())
########################################################################################################
########################################################################################################
# Thesis : Bayesian time series regression with nonparametric modeling of autocorrelation
# Project name : NSBTR_Main

                               RF.structure   = c("Gaussian")
                               True.Corr.type = c("AR","ARMA","FARIMA")
                               
                               model.Mean  = c("Linear","Nonlinear")
                               model.Sigma = c("fixed","time-varying")
                               model.Corr  = c("UK","NSBTR")
                               
                               modelkey   = c(1,1,2,2)
                                  
fname0 = paste0(model.Mean[modelkey[2]],'_',
                'var(',substr(model.Sigma[modelkey[3]],1,1),")",'_',
                model.Corr[modelkey[4]])
print(paste("                                     ",fname0))

#"                                       Linear_var(t)_NSBTR"

                                                                                      # Date : 201712xx
                                                                                      # Author : YB JUN
setwd("C:/Users/user/Desktop/WORK2019/190107_NSBTR")
source(file="NSBTR_func.R")
########################################################################################################
########################################################################################################

iniseed = 123
n.data  = 100; n.chain=3           ; n.iter = 10000; n.burn = 9000#;  n.iter = 1000; n.burn = 900;
n.train = 200; m = floor(n.train/2); n.test = 5

n = n.train + n.test
truepar = list(beta=c(1,2,3), ar = c(0.5,-0.3), ma = c(-0.6,0.2488), d = 0.25, sd = sin(pi*c(1:n)/n)+0.5)

n.beta    = length(truepar$beta)
n.theta   = (m+1) # symmetricity of periodogram
n.delta   = 7
n.eta     = n.train


########################################################################### Generate Simulation Datasets
########################################################################################################
start = Sys.time()
setwd("C:/Users/user/Desktop/Studydata/LinearAR(t)/Dataset")
seednum = SampleSetSeed(n.train,n.test,iniseed)
n = n.train + n.test

# Gererate_X
set.seed(iniseed)
X0 = matrix(1,n,1);                         set.seed(iniseed);
X1 = matrix(rnorm(n,0,1)+rnorm(n,0,2),n,1); set.seed(iniseed);
X2 = matrix(rexp(n,1),n,1)   
X  = cbind(X0,X1,X2)


for(nd in 1:n.data){
  seed = seednum[nd]; set.seed(seed)
  # Generate_Eps
  if      (modelkey[1]==1){ e =  arima.sim(model=list(ar = truepar$ar), n = n, sd = AR.mysd(truepar))
  }else if(modelkey[1]==2){ e =  arima.sim(model=list(ar = truepar$ar, ma = truepar$ma), n = n, sd = AR.mysd(truepar))
  }else if(modelkey[1]==3){
    # library(fArma)
    # e =  armaSim(model=list(d = 0.25), n = n, sd = truepar$sd)
  }
  # Generate_Y
  if      (modelkey[2]==1){ Y = X%*%truepar$beta + truepar$sd * e  
  }
  data0 = data.frame(X,e,Y)
  fname = paste(fname0,'_SimData',"(",nd,")",sep='')
  save(data0, file=fname)
}

cat("Directory:",getwd())
print(Sys.time() - start) # Time difference of 0.2170122 secs
########################################################################################


####################################################################### Initial Settings
########################################################################################
start = Sys.time()
n = n.train

nd = 1
setwd("C:/Users/user/Desktop/Studydata/LinearAR(t)/Dataset")
fname = paste(fname0,'_SimData',"(",nd,")",sep='')
load(file=fname)
data00 = data0[1:n.train,]

iniseed.chain = SampleSetSeed(n.chain,1,iniseed)
ParaList = vector(len=n.chain,mode="list")
ParaList[[1]] = Initialize_Parameters(data00,iniseed.chain[1])
ParaList[[2]] = Initialize_Parameters(data00,iniseed.chain[2])
ParaList[[3]] = Initialize_Parameters(data00,iniseed.chain[3])

fname2 = paste0(fname,"_ParaList")
setwd("C:/Users/user/Desktop/Studydata/LinearAR(t)/InitialList")
save(ParaList, file=fname2)

cat("Directory:",getwd())
print(Sys.time() - start) # Time difference of 0.431025 secs
########################################################################################


##################################################################### MCMC Gibbs Sampler
########################################################################################
start = Sys.time()
# source(file="NBTR_func.R")
setwd("C:/Users/user/Desktop/Studydata/LinearAR(t)/Dataset")
nd <- 1

# Set Estimation-result objects: EstiMat
EstiMat = list()
for(chain in 1:n.chain){
  EstiMat[[chain]] = as.mcmc(matrix(0,n.iter,(n.beta+n.theta+n.delta+n.eta)))
  colnames(EstiMat[[chain]]) <- c(paste0('beta' ,c(0:(n.beta-1))),paste0('theta',c(0:m)),
                                  paste0('delta',c(0:(n.delta-1))),paste0('logsig2.',c(0:(n.train-1))))
}
EstiMat = mcmc.list(EstiMat)

fname1 = paste(fname0,'_SimData',"(",nd,")",sep='')
fname2 = paste(fname0,'_SimData',"(",nd,")","_ParaList",sep='')
setwd("C:/Users/user/Desktop/Studydata/LinearAR(t)/Dataset")    ;load(file=fname1); data00 = data0[1:n.train,]
setwd("C:/Users/user/Desktop/Studydata/LinearAR(t)/InitialList");load(file=fname2);

Qo = CompQo(n.train)
Qc = solve(Qo)

# Gibbs Iteration
chain <- 1
for(chain in 1:n.chain){
  
  Parset0 <- ParaList[[chain]]
  
  setwd("C:/Users/user/Desktop/Studydata/LinearAR(t)")
  png(file="Initial_Results.png",width=480*2)
  par(mfrow=c(1,2))
  plot(Parset0$theta,ylim=c(-1.5,1.5),type='o',xlab="freq",ylab="",main="spd");par(new=T);plot(log(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],AR.mysd(truepar)^2,n.train)),type='l',col=2,xlab="freq",ylab="",main="spd",ylim=c(-1.5,1.5))
  plot(log(exp(Parset0$nu.eta)/sum(exp(Parset0$nu.eta))),type='o',xlab="freq",ylab="",main="eta",ylim=c(-7,-4));par(new=T);plot(log(truepar$sd[1:n.train]^2/sum(truepar$sd[1:n.train]^2)),type='l',col=2,xlab="freq",ylab="",main="eta",ylim=c(-7,-4))
  dev.off()

  iter <- 1
  for(iter in 1:n.iter){
    
    Parset1 <- Update1_beta_phi(Parset0,data00)
    Parset2 <- Update2_xi_mixture(Parset1,data00)
    Parset3 <- Update3_theta(Parset2,data00)
    Parset4 <- Update4_theta_prior(Parset3,data00)
    Parset5 <- Update5_eta(Parset4,data00)
    # Parset5 <- Update5_eta(Parset4,data00,resol=1e-02)
    # Parset6 <- Update6_eta_prior(Parset5,data00)
    newParset <- Parset5
    
    cat("iter:",iter,"  ","chain:",chain,"  ","\n")
    cat("beta=",newParset$beta,"\n")
    # cat("delta=",newParset$delta,"\n")
    cat("rho.theta=",newParset$rho.theta,"\n")
    
    EstiMat[[chain]][iter,] = c(newParset$beta,newParset$theta,newParset$delta,newParset$eta)
    
    if(iter%%10==0) cat(iter,date(),'\n')
    
    setwd("C:/Users/user/Desktop/Studydata")
    png(file=paste0("Results_chain(",chain,"_iter(",iter,").png"),width=480*2)
    par(mfrow=c(1,2))
    mse.lambda = round(mean((newParset$lambda - CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],AR.mysd(truepar)^2,n.train))^2),4)
    mse.eta  = round(mean((newParset$eta - log(truepar$sd[1:n.train]^2))^2),4)
    plot(newParset$lambda,ylim=c(0,2.5),type='o',xlab="freq",ylab="",main="");par(new=T);plot(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],AR.mysd(truepar)^2,n),ylim=c(0,2.5),type='l',col=2,xlab="freq",ylab="",main=paste0("spd(",mse.lambda,")"))
    plot(log(exp(newParset$eta)/sum(exp(newParset$eta))),type='o',ylim=c(-7,-4),xlab="freq",ylab="",main="");par(new=T);plot(log(truepar$sd[1:n.train]^2/sum(truepar$sd[1:n.train]^2)),type='l',col=2,ylim=c(-7,-4),lty=2,xlab="freq",ylab="",main=paste0("eta(",mse.eta,")"))
    dev.off()
    
    Parset0 <- newParset
    
  }
  
}

cat("Directory:",getwd())
print(Sys.time() - start)
save(EstiMat,file="190131_LCY_EstiMat.Rdata")
# save.image("C:/Users/user/Desktop/WORK2019/190207_NSBTR/190219_200sin.RData")
########################################################################################
# n.iter=100  Time difference of 3.294225 hours
# n.iter=1000 Time difference of 


# ##################################################################### Estimation Results
# ################################################  ########################################
# 
# 
# meanbeta  = matrix(0,1,n.beta)
# meantheta = matrix(0,1,(m+1))
# lower.ci  = matrix(0,1,(m+1))
# upper.ci  = matrix(0,1,(m+1))
# meaneta   = matrix(0,1,n.eta)
# eta.lower.ci  = matrix(0,1,n.eta)
# eta.upper.ci  = matrix(0,1,n.eta)
# # predoutput= matrix(0,indlen,n.test)  
# # obsoutput = matrix(0,indlen,n.test)
# 
# n.burn = 500
# nd <- 1
# # load(file=paste0("DATA",ind[i],"_paraMcMC3.Rdata"))
# 
# 
# # burn-in process
# require(coda)
# EstiMat.S = mcmc.list()
# for(chain in 1:n.chain){
#   EstiMat.S[[chain]] <- mcmc(EstiMat[[chain]][(n.burn+1):n.iter,])
#   
# }
# EstiMat.S <- mcmc.list(EstiMat.S)
# 
# # par(mfrow=c(1,2))
# # plot(temp$statistics[(3+1):(3+101),1],ylim=c(-1.5,1.5),type='o');par(new=T);plot(log(CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],mysd(modelkey,truepar)^2,n)),type='l',col=2,xlab="freq",ylab="",main="spd",ylim=c(-1.5,1.5))
# # plot(temp$statistics[(3+101+1):nrow(temp$statistics),1],ylim=c(-1.5,1.5),type='o');par(new=T);plot(log(truepar$sd^2),type='l',col=2,xlab="freq",ylab="",main="sig2",ylim=c(-1.5,1.5))
# 
# # setwd("C:/Users/user/Desktop/")
# # pdf("iteration_plot.pdf")
# # plot(EstiMat.S)
# # dev.off()
# 
# # save summaries (theta)
# sum.post = summary(EstiMat.S)
# meanbeta[nd,]  <- Re(sum.post$statistics[c(1:n.beta),1])
# meantheta[nd,] <- Re(sum.post$statistics[c((n.beta+1):(n.beta+m+1)),1])
# lower.ci[nd,]  <- sum.post$quantiles[c((n.beta+1):(n.beta+m+1)),1]
# upper.ci[nd,]  <- sum.post$quantiles[c((n.beta+1):(n.beta+m+1)),5]
# 
# 
# # plot theta
# est.mean <- meantheta[nd,]
# est.lci  <- lower.ci[nd,]
# est.uci  <- upper.ci[nd,]
# w.j      <- seq(0,m,1)/n.train
# 
# plot(w.j,exp(est.mean),type='l',ylim=c(0,4),ylab="",main="ESTIMATED SPD")
# par(new=T);plot(w.j,exp(est.lci),type='l',ylim=c(0,4),col=2,lty=2,ylab="")
# par(new=T);plot(w.j,exp(est.uci),type='l',ylim=c(0,4),col=2,lty=2,ylab="")
# par(new=T);plot(w.j,CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],AR.mysd(truepar)^2,n.train),type='l',col=4,ylim=c(0,4),lty=2,xlab="",ylab="")
# 
# 
# plot(w.j,est.mean,type='l',ylim=c(-3,1.5),ylab="",main="ESTIMATED THETA")
# par(new=T);plot(w.j,est.lci,type='l',ylim=c(-3,1.5),col=2,lty=2,ylab="")
# par(new=T);plot(w.j,est.uci,type='l',ylim=c(-3,1.5),col=2,lty=2,ylab="")
# par(new=T);plot(w.j,log( CompAR2SpectDen2(truepar$ar[1],truepar$ar[2],AR.mysd(truepar)^2,n.train) ),type='l',col=4,ylim=c(-3,1.5),lty=2,xlab="",ylab="")
# 
# 
# meaneta[nd,]   <- sum.post$statistics[c((n.beta+m+2):nrow(sum.post$statistics)),1]
# eta.lower.ci[nd,]  <- sum.post$quantiles[c((n.beta+m+2):nrow(sum.post$statistics)),1]
# eta.upper.ci[nd,]  <- sum.post$quantiles[c((n.beta+m+2):nrow(sum.post$statistics)),5]
# est.mean <- meaneta[nd,]
# est.lci  <- eta.lower.ci[nd,]
# est.uci  <- eta.upper.ci[nd,]
# 
# plot(c(1:n),exp(est.mean),type='l',ylim=c(0,10),ylab="",main="ESTIMATED MAR.VAR")
# par(new=T);plot(c(1:n),exp(est.lci),type='l',ylim=c(0,10),col=2,lty=2,ylab="")
# par(new=T);plot(c(1:n),exp(est.uci),type='l',ylim=c(0,10),col=2,lty=2,ylab="")
# par(new=T);plot(c(1:n),truepar$sd[1:n.train]^2,type='l',col=4,ylim=c(0,10),lty=2,xlab="",ylab="")
# 
# 
# plot(c(1:n),est.mean,type='l',ylim=c(-2,2),ylab="",main="ESTIMATED ETA")
# par(new=T);plot(c(1:n),est.lci,type='l',ylim=c(-2,2),col=2,lty=2,ylab="")
# par(new=T);plot(c(1:n),est.uci,type='l',ylim=c(-2,2),col=2,lty=2,ylab="")
# par(new=T);plot(c(1:n),log( truepar$sd[1:n.train]^2 ),type='l',col=4,ylim=c(-2,2),lty=2,xlab="",ylab="")
# 
# mean( ( log( truepar$sd[1:n.train]^2 ) - est.mean )^2 )












