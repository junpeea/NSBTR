#source('functions_cor.R')

### this is for y=\eta(x) + \sigma \eps case  #############################
### so that spectral density is for correlation############################
### that is \int lambda = 1
#length(global.var$beta.mu0)

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

CompAR2SpectDen =function(a1,a2,sigma2,n){
#compute AR2 spectral density  with w in (0,1/2)
#for AR1, let a2=0
#  eps_t -a1 eps_{t-1} -a2 eps_{t-2} = sqrt(sigma2)* z_t
#plot(CompAR2SpectDen(0.5,-0.3,0.5^2,n))

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


CompFARIMASpectDen=function(d,sigma2,n){
#compute spectral density of FARIMA(0,d,0) with var=sigma2
# with w in (0,1/2)
#plot(CompFARIMASpectDen(0.25,0.5^2,200))


		k <- seq(0,floor(n/2),1)
		w.j <- k/n

		res= ((1-cos(2*pi*w.j))^2+(sin(2*pi*w.j))^2)^d
		return(sigma2/res)


}



CompLogPerioFunc = function(epsi){
#sub-routine to get log of peridogram

	n <- length(epsi)
	t.peri <- seq(1,n,1)
	k <- seq(0,floor(n/2),1)
	w.j <- k/n

	peri.temp = function(w.j)
	{ 
		((sum(epsi * cos(2* pi * t.peri * w.j)))^2 + (sum(epsi * sin(2* pi * t.peri * w.j)))^2)/n
	}

	peri <- log(sapply(w.j,peri.temp))
	return(peri)
}




CompQFunc_old=function(n){

		c.q  = s.q = matrix(0, ncol = floor(n/2), nrow = n)
		QQ <- matrix(rep(0,n*(n-1)),nrow=n)

		for(i in 1:n){
			for(j in 1:dim(c.q)[2])
				{
				c.q[i,j] <-	(sqrt(2))* cos((2 * pi * j * i)/n)
				s.q[i,j] <-	(sqrt(2))* sin((2 * pi * j * i)/n)
				}
			     }

#print(c.q)
#print(s.q)


		if(n%%2 ==0){


			for(i in 1:n-1){
	
				if (i %%2 !=0) QQ[,i] = c.q[,(i+1)/2]
				else QQ[,i] = s.q[,(i+1)/2]
				}

		QQ[,n-1] <- 2^(-1/2) * QQ[,n-1]

		QQ <- n^(-1/2) * cbind(rep(1,length(n)),QQ)

		}



	else{

		for(i in 1:n-1){
			if (i %%2 !=0) QQ[,i] <- c.q[,(i+1)/2]
			else QQ[,i] <- s.q[,(i+1)/2]
		}



			QQ <- n^(-1/2) * cbind(rep(1,length(n)),QQ)
		}


	return(QQ)
}



CompQ0Func=function(n){
#as of 02/26/2013, Q_{u,v}=(1/sqrt(n))*exp(i(u-1)(v-1)2pi/n)
#Note that it is a complex-valued matrix.
#we need to calculate \Gamma_n =Q \Lambda_n Q^*
#Note that the resulting matrix is real-valued matrix
#Let Q_0 = (u-1)*(v-1)*2*pi/n matrix
#Then, \Gamma_n = (1/sqrt(n))*cosQ_0 \Lambda_n (1/sqrt(n))*cos Q_0^T
#               + (1/sqrt(n))*sinQ_0 \Lambda_n (1/sqrt(n))*sin Q_0^T
#so, let's just compute Q_0 and then do \Gamma_n calculation in the main code
#

		
		umat=matrix(0:(n-1),n,1)%*%matrix(1,1,n)
		vmat=t(umat)
		QQ=umat*vmat*2*pi/n		


	return(QQ)
}


CompThetaPriorFunc=function(n,tau2,rho){
#function to set mean and covariance matrix of prior of theta
#input:
#global.var: information that contain all necessary parameter settings
#output:
#nu.n: prior mean 
#tau.n: prior covariance
	
## creating matrix tau.nn
## covariance kernel is of the following form
## sigma(w1,w2) = [exp(-10(w1-w2)^2)]/tau^2
## please have a look at the matrix tau.nn value!
##########################################
	k <- seq(0,floor(n/2),1)
	w.j <- k/n

### this is the most sensitive part
### the value of tau.nn formula 
### const * exp(-con*(w.j[i] - w.j[j])^2)
### intially const =1000, con = 100
      
	#const=1000;con=100
	tau.nn = matrix(0,length(w.j),length(w.j))
	dist0=as.matrix(dist(cbind(w.j,w.j)))
	tau.nn=(1/tau2)*exp(-rho*dist0)
#     distance should be abs() instead of ()^2
#	for(i in 1:length(w.j)){
#		for(j in 1:length(w.j)){
#		tau.nn[i,j] = (1/tau2)*exp(-rho*(w.j[i] - w.j[j])^2)
#		}
#	
	

	#prior mean, nu.n 
	nu.n=matrix(0, dim(tau.nn)[1],1)

	return(list(tau.nn=tau.nn, nu.n=nu.n))
}

	

ConvLambdaToThetaFunc=function(m,lambda){
#function to convert lambda to theta
#lambda: spectral density 1:n
#theta=log(lambda) 1:m=floor(n/2)

theta=as.matrix(log(lambda[1:m]))

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


#x.f=cbind(1,data0$data[(n-kstep+1):n,2:n.beta])
#Q.mat.f=global.var$Q.mat.f
#eta.x.f=as.matrix(x.f,kstep,n.beta) %*% as.matrix(ConvMat[[chain]][4918,1:n.beta])
#n=200;n.theta=101
#theta=ConvMat[[chain]][4918,(n.beta+1):(n.beta+n.theta)]
#lambda=ConvThetaToLambdaFunc(n,theta)
#gamma.inv.nn= global.var$Q.mat %*% diag(1/as.vector(lambda)) %*% t(global.var$Q.mat)



ForecastFunc=function(all,x.f,kstep,Q0.mat.f){
#function to forecast y value in f-step ahead
#x.f=c(1,2,3);kstep=1



n=dim(all$data)[1]
epsi.0=all$epsi.0

cosQ0=cos(global.var$Q0.mat)
sinQ0=sin(global.var$Q0.mat)
gamma.inv.nn.1= cosQ0%*%diag(1/as.vector(all$lambda))%*%t(cosQ0)/n
gamma.inv.nn.2= sinQ0%*%diag(1/as.vector(all$lambda))%*%t(sinQ0)/n
gamma.inv.nn=(all$tau2.e)*(gamma.inv.nn.1+gamma.inv.nn.2)
#gamma.inv.nn= (all$tau2.e)*global.var$Q.mat %*% diag(1/as.vector(all$lambda)) %*% t(global.var$Q.mat)


if(global.var$RegModel=='linear'){

	n.beta=length(all$beta)
	#compute eta
	eta.x.f=as.matrix(x.f,kstep,n.beta) %*% all$beta
}

if(global.var$RegModel=='additive'){
	n.x=dim(all$data)[2]-1
	n.beta=length(all$beta)
	max.mat=t(matrix(all$data.max[-1],n.x,1)%*%matrix(2,1,kstep))
	min.mat=t(matrix(all$data.min[-1],n.x,1)%*%matrix(2,1,kstep))
	x.f.star0=(x.f[,-1]-min.mat)/(max.mat-min.mat)
	
	x.f.star=matrix(1,kstep,1)
	for(j in 1:n.x){
		tmp.bs=bs(all$data.s[,j+1],knots=global.var$knots.p)
		x.f.star=cbind(x.f.star,predict(tmp.bs,x.f.star0[,j]))
	}
	#compute eta
	eta.x.f=as.matrix(x.f.star,kstep,n.beta) %*% all$beta
}

#this only works for two covariates since it is using bivariate spline
#
if(global.var$RegModel=='nonparametric'){
	nx=dim(x.f)[2]
	if(nx<=2){
		stop('less than two covariates!!')
	}
	x.f.mat=x.f[,2:nx]
	x.mean.mat=t(as.matrix(all$data.mean[2:nx]))%x%matrix(1,kstep,1)
	x.sd.mat=t(as.matrix(all$data.sd[2:nx]))%x%matrix(1,kstep,1)
	x.f.mat.s=(x.f.mat-x.mean.mat)/x.sd.mat
		
	require(akima) #package for bivariate spline
	eta.x.f=interpp(all$data.s[,2],all$data.s[,3],all$eta,xo=x.f.mat.s[,1],yo=x.f.mat.s[,2])$z


}


y.f=c()
for(ks in 1:kstep){

f=ks+n

#compute Gamma.f for new fourier frequency range
fourier.f= (0:(f-1))/f
lambda.f=diag(spline((0:(n-1))/n,all$lambda,xout=fourier.f)$y)

cosQ0=cos(Q0.mat.f[[ks]])
sinQ0=sin(Q0.mat.f[[ks]])
Gamma.f.1= cosQ0%*%lambda.f%*%t(cosQ0)/n
Gamma.f.2= sinQ0%*%lambda.f%*%t(sinQ0)/n
Gamma.f=(1/all$tau2.e)*(Gamma.f.1+Gamma.f.2)
#Gamma.f=(1/all$tau2.e)*Q.mat.f[[ks]] %*% lambda.f %*% t(Q.mat.f[[ks]])

#compute h.f
h.f=as.matrix(rev(Gamma.f[(ks+1):f,1]))
#compute eps.f
eps.f=t(h.f)%*%gamma.inv.nn%*%epsi.0
y.f[ks]=eta.x.f[ks]+eps.f

}


return(y.f)

}

#all=all[[chain]]; x.f=cbind(1,data[(n-kstep+1):n,2:n.beta]);Q0.mat.f=global.var$Q0.mat.f


ForecastFunc2=function(all,x.f,kstep,global.var){
  #function to forecast y value in f-step ahead
  #x.f=c(1,2,3);kstep=1
  
  
  #Q0.mat.f=global.var$Q0.mat.f
  
  n=dim(all$data)[1]
  epsi.0=all$epsi.0
  
  
    
    
  cosQ0=cos(global.var$Q0.mat)
  sinQ0=sin(global.var$Q0.mat)
  gamma.inv.nn.1= cosQ0%*%diag(1/as.vector(all$lambda))%*%t(cosQ0)/n
  gamma.inv.nn.2= sinQ0%*%diag(1/as.vector(all$lambda))%*%t(sinQ0)/n
  gamma.inv.nn=(all$tau2.e)*(gamma.inv.nn.1+gamma.inv.nn.2)
  #gamma.inv.nn= (all$tau2.e)*global.var$Q.mat %*% diag(1/as.vector(all$lambda)) %*% t(global.var$Q.mat)
  
  
  if(global.var$RegModel=='linear'){
    
    n.beta=length(all$beta)
    #compute eta
    eta.x.f=as.matrix(x.f,kstep,n.beta) %*% all$beta
  }
  
  if(global.var$RegModel=='additive'){
    n.x=dim(all$data)[2]-1
    n.beta=length(all$beta)
    max.mat=t(matrix(all$data.max[-1],n.x,1)%*%matrix(2,1,kstep))
    min.mat=t(matrix(all$data.min[-1],n.x,1)%*%matrix(2,1,kstep))
    x.f.star0=(x.f[,-1]-min.mat)/(max.mat-min.mat)
    
    x.f.star=matrix(1,kstep,1)
    for(j in 1:n.x){
      tmp.bs=bs(all$data.s[,j+1],knots=global.var$knots.p)
      x.f.star=cbind(x.f.star,predict(tmp.bs,x.f.star0[,j]))
    }
    #compute eta
    eta.x.f=as.matrix(x.f.star,kstep,n.beta) %*% all$beta
  }
  
  #this only works for two covariates since it is using bivariate spline
  #
  if(global.var$RegModel=='nonparametric'){
    nx=dim(x.f)[2]
    if(nx<=2){
      stop('less than two covariates!!')
    }
    x.f.mat=x.f[,2:nx]
    x.mean.mat=t(as.matrix(all$data.mean[2:nx]))%x%matrix(1,kstep,1)
    x.sd.mat=t(as.matrix(all$data.sd[2:nx]))%x%matrix(1,kstep,1)
    x.f.mat.s=(x.f.mat-x.mean.mat)/x.sd.mat
    
    require(akima) #package for bivariate spline
    eta.x.f=interpp(all$data.s[,2],all$data.s[,3],all$eta,xo=x.f.mat.s[,1],yo=x.f.mat.s[,2])$z
    
    
  }
  
  
  step.size=1000
  w.j.n=seq(0,floor(step.size/2),1)/step.size
  lambda.n=spline(seq(0,floor(n/2),1)/n,all$lambda[1:(floor(n/2)+1)],xout=w.j.n)$y
  
  
  y.f=c()
  for(ks in 1:kstep){
    
    f=ks+n
  
    umat=matrix(1,length(w.j.n),1)%*%t(matrix(rev(seq(ks,ks+n-1,by=1))))
    w.j.mat=matrix(w.j.n)%*%matrix(1,1,dim(umat)[2])
    cos.mat=cos(2*pi*umat*w.j.mat)
    h.f=(1/all$tau2.e)*t(t(matrix(lambda.n))%*%cos.mat)*2/step.size
    
    
    #compute Gamma.f for new fourier frequency range
    #fourier.f= (0:(f-1))/f
    #lambda.f=diag(spline((0:(n-1))/n,all$lambda,xout=fourier.f)$y)
    
    #cosQ0=cos(Q0.mat.f[[ks]])
    #sinQ0=sin(Q0.mat.f[[ks]])
    #Gamma.f.1= cosQ0%*%lambda.f%*%t(cosQ0)/n
    #Gamma.f.2= sinQ0%*%lambda.f%*%t(sinQ0)/n
    #Gamma.f=(1/all$tau2.e)*(Gamma.f.1+Gamma.f.2)
    
    
    #compute h.f
    #h.f=as.matrix(rev(Gamma.f[(ks+1):f,1]))
    #compute eps.f
    eps.f=t(h.f)%*%gamma.inv.nn%*%epsi.0
    y.f[ks]=eta.x.f[ks]+eps.f
    
  }
  
  
  return(y.f)
  
}


ForecastOnlyFunc=function(beta,theta,epsi.0,tau2.e,x.f,kstep,Q0.mat.f){
#function to forecast y value in f-step ahead
#x.f=c(1,2,3);kstep=1



n=length(epsi.0)
lambda=ConvThetaToLambdaFunc(n,theta)

cosQ0=cos(global.var$Q0.mat)
sinQ0=sin(global.var$Q0.mat)
gamma.inv.nn.1= cosQ0%*%diag(1/as.vector(lambda))%*%t(cosQ0)/n
gamma.inv.nn.2= sinQ0%*%diag(1/as.vector(lambda))%*%t(sinQ0)/n
gamma.inv.nn=(tau2.e)*(gamma.inv.nn.1+gamma.inv.nn.2)
#gamma.inv.nn= (tau2.e)*global.var$Q.mat %*% diag(1/as.vector(lambda)) %*% t(global.var$Q.mat)

#if(global.var$RegModel=='linear'){

	n.beta=length(beta)
	#compute eta
	eta.x.f=as.matrix(x.f,kstep,n.beta) %*% matrix(beta,n.beta,1)
#}


y.f=c()
for(ks in 1:kstep){

f=ks+n

#compute Gamma.f for new fourier frequency range
fourier.f= (0:(f-1))/f
lambda.f=diag(spline((0:(n-1))/n,lambda,xout=fourier.f)$y)

cosQ0=cos(Q0.mat.f[[ks]])
sinQ0=sin(Q0.mat.f[[ks]])
Gamma.f.1= cosQ0%*%lambda.f%*%t(cosQ0)/n
Gamma.f.2= sinQ0%*%lambda.f%*%t(sinQ0)/n
Gamma.f=(1/tau2.e)*(Gamma.f.1+Gamma.f.2)
#Gamma.f=(1/tau2.e)*Q.mat.f[[ks]] %*% lambda.f %*% t(Q.mat.f[[ks]])

#compute h.f
h.f=as.matrix(rev(Gamma.f[(ks+1):f,1]))
#compute eps.f
eps.f=t(h.f)%*%gamma.inv.nn%*%epsi.0
y.f[ks]=eta.x.f[ks]+eps.f

}

return(y.f)

}


ForecastFunc_Par=function(all,x.f,kstep,global.var){
#function to forecast y value in f-step ahead
#y.f[,1]= prediction with conditional mean
#y.f[,2]= prediction with predictive y

n.beta=length(all$beta)
#compute eta
eta.x.f=as.matrix(x.f,kstep,n.beta) %*% all$beta
n=dim(all$data)[1]
new.epsi=as.matrix(c(all$epsi,array(0,kstep)))
new.epsi2=cbind(new.epsi,new.epsi)
phi=all$ar.phi
sigma2=all$sigma2
ar.order=global.var$ar.order
y.f=matrix(0,kstep,2)
seednum=all$seednum

seednum=seednum+10;set.seed(seednum)
err.Z=rnorm(kstep,0,1)

for(ks in 1:kstep){
new.epsi2[n+ks,1]= sum(phi*new.epsi2[(n+ks-1):(n+ks-ar.order),1])
new.epsi2[n+ks,2]= sum(phi*new.epsi2[(n+ks-1):(n+ks-ar.order),2]) + sqrt(sigma2)*err.Z[ks]
}


y.f[,1]=eta.x.f+new.epsi2[(n+1):(n+kstep),1]
y.f[,2]=eta.x.f+new.epsi2[(n+1):(n+kstep),2]

return(y.f)

}


#ForecastFunc_Par2(all[[chain]],cbind(1,data[(n-kstep+1):n,2:n.beta]),kstep,global.var)


ForecastFunc_Par2=function(all,x.f,kstep,global.var){
#function to forecast y value in f-step ahead using parametric model
# and conditional mean approach same as nonparametric case
# i.e. y_f = \hat\eta(x.f)+E(\epsilon_f | \epsilon_0)

n.beta=length(all$beta)
#compute residual
epsi=all$epsi

eta.x.f=as.matrix(x.f,kstep,n.beta) %*% all$beta
n=dim(all$data)[1]
gamma.v=array(0,n+kstep)
phi=all$ar.phi
sigma2=all$sigma2
ar.order=global.var$ar.order

if(ar.order==2){#AR(2) process
#first calculate gamma(0), gamma(1), gamma(2), where
#gamma(k) is the auto-covariance ft at lag k
tmp.X=matrix(0,3,3)
diag(tmp.X)=1
tmp.X[1,2]=tmp.X[2,1]=tmp.X[3,2]=-phi[1]
tmp.X[1,3]=tmp.X[2,3]=tmp.X[3,1]=-phi[2]
tmp.A=matrix(c(sigma2,0,0),3,1)
tmp.gamma=ginv(tmp.X)%*%tmp.A

gamma.v[1:2]=tmp.gamma[1:2]
for(i in 3:(n+kstep)){
#gamma(0),gamma(1), ... 
gamma.v[i]=phi[1]*gamma.v[i-1]+phi[2]*gamma.v[i-2]
}
} 

if(ar.order==1){#AR(1) process
#first calculate gamma(0), gamma(1), where
#gamma(k) is the auto-covariance ft at lag k
tmp.X=matrix(0,2,2)
diag(tmp.X)=1
tmp.X[1,2]=tmp.X[2,1]=-phi[1]
tmp.A=matrix(c(sigma2,0),2,1)
tmp.gamma=ginv(tmp.X)%*%tmp.A

gamma.v[1]=tmp.gamma[1]
for(i in 2:(n+kstep)){
#gamma(0),gamma(1), ... 
gamma.v[i]=phi[1]*gamma.v[i-1]
}
} 

#create Gamma matrix (cov matrix)
gamma.mat=matrix(0,n,n)
for(i in 1:n){
if(i==1){
gamma.mat[i,]=gamma.v[1:n]
}else{
gamma.mat[i,]=c(rev(gamma.v[2:i]),gamma.v[1:(n-i+1)])
}
}

gamma.inv.nn=ginv(gamma.mat)

y.f=c()
for(ks in 1:kstep){

f=ks+n

#compute h.f
h.f=as.matrix(rev(gamma.v[(ks+1):f]))
#compute eps.f
eps.f=t(h.f)%*%gamma.inv.nn%*%epsi
y.f[ks]=eta.x.f[ks]+eps.f

}

return(y.f)

}




ForecastOnlyFunc_Par=function(beta,epsi,phi,sigma2,x.f,kstep,global.var){
#function to forecast y value in f-step ahead using parametric model
# and conditional mean approach same as nonparametric case
# i.e. y_f = \hat\eta(x.f)+E(\epsilon_f | \epsilon_0)

n.beta=length(global.var$beta.mu0)

eta.x.f=as.matrix(x.f,kstep,n.beta) %*%matrix(beta,n.beta,1)
n=length(epsi)
gamma.v=array(0,n+kstep)
ar.order=global.var$ar.order

if(ar.order==2){#AR(2) process
#first calculate gamma(0), gamma(1), gamma(2), where
#gamma(k) is the auto-covariance ft at lag k
tmp.X=matrix(0,3,3)
diag(tmp.X)=1
tmp.X[1,2]=tmp.X[2,1]=tmp.X[3,2]=-phi[1]
tmp.X[1,3]=tmp.X[2,3]=tmp.X[3,1]=-phi[2]
tmp.A=matrix(c(sigma2,0,0),3,1)
tmp.gamma=solve(tmp.X)%*%tmp.A

gamma.v[1:2]=tmp.gamma[1:2]
for(i in 3:(n+kstep)){
#gamma(0),gamma(1), ... 
gamma.v[i]=phi[1]*gamma.v[i-1]+phi[2]*gamma.v[i-2]
}
} 

if(ar.order==1){#AR(1) process
#first calculate gamma(0), gamma(1), where
#gamma(k) is the auto-covariance ft at lag k
tmp.X=matrix(0,2,2)
diag(tmp.X)=1
tmp.X[1,2]=tmp.X[2,1]=-phi[1]
tmp.A=matrix(c(sigma2,0),2,1)
tmp.gamma=solve(tmp.X)%*%tmp.A

gamma.v[1]=tmp.gamma[1]
for(i in 2:(n+kstep)){
#gamma(0),gamma(1), ... 
gamma.v[i]=phi[1]*gamma.v[i-1]
}
} 

#create Gamma matrix (cov matrix)
gamma.mat=matrix(0,n,n)
for(i in 1:n){
if(i==1){
gamma.mat[i,]=gamma.v[1:n]
}else{
gamma.mat[i,]=c(rev(gamma.v[2:i]),gamma.v[1:(n-i+1)])
}
}

gamma.inv.nn=solve(gamma.mat)

y.f=c()
for(ks in 1:kstep){

f=ks+n

#compute h.f
h.f=as.matrix(rev(gamma.v[(ks+1):f]))
#compute eps.f
eps.f=t(h.f)%*%gamma.inv.nn%*%epsi
y.f[ks]=eta.x.f[ks]+eps.f

}

return(y.f)

}




GenerateDataFunc=function(datasize,RegModel,ErrModel,beta=c(1,2,3),
ar.coef=c(0.5,-0.3), ma.coef=c(-0.6,0.2488),err.sd=0.5,farima.d=0.25,seednum=c()){
#function to generate data for the simulation 
#input: 
#datasize: size of data
#kstep: step for forecasting
#RegModel: model for regression function
#beta: parameter for regression function
#ErrModel: model for error process
#seednum: a seed number to produce same data
#what model is generated?
#RegModel =="linear":    eta= beta0+beta1 x1 + beta2 x2
#	    =="nonlinear": eta= beta0+beta1 x1*x2
#         =="nonlinear2": eta= beta0+beta1 sin(x1) + beta2 log(x2)
#ErrModel =="ar":    epsilon ~ AR(2)
#         =="farima"   epsilon ~ FARIMA(0,d=0.25,0)
#As a default, x1 ~ N(0,1), x2 ~ exp(1)
#                            	
	require(stats)
	if(is.null(seednum)){
		seednum=floor(proc.time()[3]/60)
	}

	#generate x1 and x2
	set.seed(1234)
	x1= matrix(rnorm(datasize,0,1),datasize,1)
	set.seed(2345)
	x2= matrix(rexp(datasize),datasize,1)	

	#compuate regression function
	if(RegModel=='linear'){
		eta=beta[1]+beta[2]*x1+beta[3]*x2
		beta0=beta
	}
	if(RegModel=='nonlinear'){
		eta=beta[1]+beta[2]*x1*x2
		beta0=beta
	}

	if(RegModel=='nonlinear2'){
		eta=beta[1]+beta[2]*sin(x1)+ beta[3]*log(x2)
		beta0=beta
	}


	#generate error process
	if(ErrModel=='ar'){ 
		seednum=seednum+10;set.seed(seednum)
		eps=arima.sim(model=list(ar=ar.coef,ma=c(0,0)),n=datasize,sd=err.sd)
		err.par=list(ar.coef=ar.coef, sd=err.sd)
	}
	if(ErrModel=='farima'){
		seednum=seednum+10
		set.seed(seednum)
		eps=farima.generate(500, datasize, farima.d, err.sd, seednum)

		err.par=list(d=farima.d, sd=err.sd)
	}


	if(ErrModel=='arma'){
		seednum=seednum+10;set.seed(seednum)
		eps=arima.sim(model=list(ar=ar.coef[1],ma=ma.coef[1]),n=datasize,sd=err.sd)
		err.par=list(ar.coef=ar.coef[1], ma.coef=ma.coef[1], sd=err.sd)
	}

	#compute y
      y=eta+eps

	res=data.frame(cbind(y,x1,x2))

	res2=list(data=res, eta=eta, eps=eps
                , beta=beta0, RegModel=RegModel, ErrModel=ErrModel
                , ErrPar=err.par,seednum=seednum)

	return(res2)
}

InitializeFunc=function(data,global.var){
#function to get initial parameter values
#Input:
#data: data matrix
#global.var: global variables

	nn=length(data[,1])-global.var$kstep
	y=matrix(data[1:nn,1],nn,1)

	
	#contain all necessary component
	all=list(seednum=global.var$seednum
                ,data=data[1:nn,])


	if(global.var$RegModel=='linear'){#linear regression function
		x=as.matrix(cbind(1,data[1:nn,-1]))
		XX=as.matrix(t(x)%*%x)
		all$beta = chol2inv(chol(XX)) %*% (t(x) %*% y)
		#compute initial eta
		all$eta= x %*% all$beta #for linear model
	}
	if(global.var$RegModel=='additive'){ #additive regression function
		np=dim(data)[2]-1
		#knots.p= seq(0.2,0.8,by=0.2)
		max.mat=t(as.matrix(apply(data,2,FUN=max))%x%matrix(1,1,nn))
		min.mat=t(as.matrix(apply(data,2,FUN=min))%x%matrix(1,1,nn))
		all$data.s=(all$data-min.mat)/(max.mat-min.mat)
		all$data.max=apply(data,2,FUN=max)
		all$data.min=apply(data,2,FUN=min)
		
		x=matrix(1,nn,1)
		for(j in 1:np){
			x=cbind(x,bs(all$data.s[,j+1],knots=global.var$knots.p))
		}
		XX=as.matrix(t(x)%*%x)
		all$beta = chol2inv(chol(XX)) %*% (t(x) %*% y)
		#compute initial eta
		all$eta= x %*% all$beta #for linear model

		all$x=x
		all$XX=XX
		all$n.spline=(dim(x)[2]-1)/2
		#all$knots.p=knots.p
	}
	if(global.var$RegModel=='nonparametric'){

		#compute initial eta
		all$eta= global.var$eta.mu0

		all$eta.mu0=global.var$eta.mu0
		all$eta.tau2=global.var$eta.tau2
		all$eta.rho=global.var$eta.rho

		########################################################################################
		#all$data.s: scaled data after taking mean and sd
		#########################################################################################
		mean.mat=t(as.matrix(apply(all$data,2,FUN=mean))%x%matrix(1,1,nn))
		sd.mat=t(as.matrix(apply(all$data,2,FUN=sd))%x%matrix(1,1,nn))
		all$data.s=(all$data-mean.mat)/sd.mat
		all$data.mean=apply(all$data,2,FUN=mean)
		all$data.sd=apply(all$data,2,FUN=sd)
		
		#compute all$eta.rho.dist.inv: to update parameters in Eta's Cov matrix 
		dist0=as.matrix(dist(all$data.s[,-1]))
		exp.dist0=array(exp(-t(as.matrix(global.var$eta.rho.vec))%x% dist0),c(dim(dist0)[1],dim(dist0)[2],length(global.var$eta.rho.vec)))
		exp.dist0.inv=exp.dist0
		for(ri in 1:length(global.var$eta.rho.vec)){
			exp.dist0.inv[,,ri]=chol2inv(chol(exp.dist0[,,ri]))
		}
		all$eta.rho.dist.inv=exp.dist0.inv




	}


	j <- seq(0,floor(nn/2),1)
	w.j <- j/nn
	dist0=as.matrix(dist(cbind(w.j,w.j)))
	exp.dist0=array(exp(-t(as.matrix(global.var$rho.vec))%x% dist0),c(dim(dist0)[1],dim(dist0)[2],length(global.var$rho.vec)))
	exp.dist0.inv=exp.dist0
	for(ri in 1:length(global.var$rho.vec)){
		exp.dist0.inv[,,ri]=chol2inv(chol(exp.dist0[,,ri]))

	}
	all$rho.dist.inv=exp.dist0.inv

	
	#all$epsi.0 is original residual
	all$epsi.0 = y - all$eta

	# precision parameter for the error
	all$tau2.e=as.numeric(1/var(all$epsi.0))

	#all$epsi is scaled residual
	all$epsi = (all$epsi.0)*sqrt(all$tau2.e)

	all$log.peri = CompLogPerioFunc(all$epsi)

	#initial for theta: smooth curve using log.peri
	require(KernSmooth)
	tmp.lambda=ksmooth(w.j,exp(all$log.peri),x.points=w.j,bandwidth=0.1)
	all$theta=log(tmp.lambda$y)
#	all$theta= as.matrix(rep(mean(all$log.peri), length(all$log.peri)))
	j <- seq(0,floor(nn/2),1)
	w.j <- j/nn
	m=length(j)

	
	#set up initial guess for prior mean (used mean of log.peri)
	#all$nu.n=all$theta 
	#tmp1=ksmooth(w.j,exp(all$log.peri),x.points=w.j,bandwidth=0.1)
	#par(mfrow=c(1,2))
	#plot(w.j,exp(all$log.peri),col='red')
	#lines(w.j,(tmp1$y),col='blue')
	#plot(w.j,(all$log.peri))
	#lines(w.j,log(tmp1$y))	
	#plot(w.j,all$log.peri-log(tmp1$y),ylim=c(-10,10))


	#contain mixture component information 
	all = UpdateMixtureFunc(all)


	
	#compute prior mean and covariance for theta
	all$tau2=global.var$tau2
	all$rho=global.var$rho
	theta.pr=CompThetaPriorFunc(nn,all$tau2,all$rho)
	all$tau.nn=theta.pr$tau.nn
	

	#convert theta on w_0, ... w_{n/2+1} to lambda on w_0, ... w_{n-1}
	all$lambda = ConvThetaToLambdaFunc(nn,all$theta)
	all$theta=ConvLambdaToThetaFunc(length(all$theta),all$lambda)
	all$nu.n=all$theta

	return(all)


}


InitializeFunc_Par=function(data,global.var,ar.order){
#function to get initial parameter values
#Input:
#data: data matrix
#global.var: global variables
#ar.order: order of AR model
	
	nn=length(data[,1])-global.var$kstep
	y=matrix(data[1:nn,1],nn,1)
	if(global.var$RegModel=='linear'){#linear regression function
		x=as.matrix(cbind(1,data[1:nn,-1]))
		XX=as.matrix(t(x)%*%x)
		ini.beta = chol2inv(chol(XX)) %*% (t(x) %*% y)
		epsi.temp = y - x %*% ini.beta
	}
	if(global.var$RegModel=='nonlinear'){
		#use spline
		#need ini.beta and epsi.temp
	}
	
	#initial ar coef
	ar.phi=matrix(0,ar.order,1)

	#initial sigma2
	seednum=global.var$seednum
	seednum=seednum+10;set.seed(seednum)
	sigma2.inv=1 #rgamma(1,shape=global.var$sigma2.hyp[1],
             #        scale=global.var$sigma2.hyp[2])
	#initial eps -p+1  eps -p+2 ... eps 0
	seednum=seednum+10;set.seed(seednum)
	eps0=as.matrix(rnorm(global.var$ar.order,0,sqrt(global.var$tau2)))
	#contain all necessary component
	all=list(seednum=seednum
                	,beta=ini.beta, epsi=epsi.temp
			,ar.phi=ar.phi
			,sigma2=1/sigma2.inv
			,eps0=eps0
			,data=data[1:nn,])

	return(all)


}




PerturbFunc=function(all0,n.chain,global.var){
#function to perturb initial values
	res=list(c())
	seednum=all0$seednum
	for(chain in 1:n.chain){
		res[[chain]]=all0

		if (global.var$RegModel=='linear'){
			seednum=seednum+10;set.seed(seednum)
			temp.beta=all0$beta+rnorm(length(all0$beta),0,1)	
			res[[chain]]$beta=as.matrix(temp.beta)
			x=as.matrix(cbind(1,all0$data[,-1]))
			XX=as.matrix(t(x)%*%x)
			res[[chain]]$eta= x %*% as.matrix(temp.beta) #for linear model
			res[[chain]]$epsi.0 = as.matrix(all0$data[,1]) -res[[chain]]$eta 
			res[[chain]]$epsi = sqrt(res[[chain]]$tau2.e)*(res[[chain]]$epsi.0)
			res[[chain]]$log.peri = CompLogPerioFunc(res[[chain]]$epsi)
		}
		if(global.var$RegModel=='additive'){
			seednum=seednum+10;set.seed(seednum)
			temp.beta=all0$beta+rnorm(length(all0$beta),0,1)	
			res[[chain]]$beta=as.matrix(temp.beta)
			res[[chain]]$eta= all0$x %*% as.matrix(temp.beta) #for linear model
			res[[chain]]$epsi.0 = as.matrix(all0$data[,1]) -res[[chain]]$eta 
			res[[chain]]$epsi = sqrt(res[[chain]]$tau2.e)*(res[[chain]]$epsi.0)
			res[[chain]]$log.peri = CompLogPerioFunc(res[[chain]]$epsi)
		}
		if (global.var$RegModel=='nonparametric'){
			seednum=seednum+10;set.seed(seednum)
			res[[chain]]$eta=all0$eta+rnorm(length(all0$eta),0,1)
			res[[chain]]$epsi.0 = as.matrix(all0$data[,1]) -res[[chain]]$eta 
			res[[chain]]$epsi = sqrt(res[[chain]]$tau2.e)*(res[[chain]]$epsi.0)
			res[[chain]]$log.peri = CompLogPerioFunc(res[[chain]]$epsi)
		}
	
		seednum=seednum+10;set.seed(seednum)
		temp.theta=all0$theta+0.1*rnorm(length(all0$theta),0,1)
		res[[chain]]$theta=as.matrix(temp.theta)
		res[[chain]]$lambda = ConvThetaToLambdaFunc(length(res[[chain]]$eta),res[[chain]]$theta)

		res[[chain]]$seednum=seednum+123*chain
	}
	return(res)

}

PerturbFunc_Par=function(all0,n.chain,global.var){
#function to perturb initial values
	res=list(c())
	seednum=all0$seednum
	for(chain in 1:n.chain){
		seednum=seednum+10;set.seed(seednum)
		temp.beta=all0$beta+rnorm(length(all0$beta),0,1)	
		seednum=seednum+10;set.seed(seednum)
		temp.sigma2=all0$sigma2+abs(rnorm(1,0,1))
		res[[chain]]=all0
		res[[chain]]$beta=as.matrix(temp.beta)
		#res[[chain]]$sigma2=temp.sigma2
		res[[chain]]$seednum=seednum+123*chain
	}
	return(res)

}



UpdateMixtureFunc = function(all,global.var)
{
# sub-routine for finding mixture component of normal's
# xi = logperi - theta

      xi = all$log.peri - all$theta
	seednum=all$seednum
	#this is for log X_1^2 for w_j =0 or \pi
	phi.prob.temp1 <-c(0.13,0.16,0.23,0.22,0.25)
	phi.mean.temp1 <- c(-4.63, -2.87, -1.44, -0.33, 0.76)
	phi.var.temp1 <- c(8.75, 1.95, 0.88, 0.45, 0.41)
	#this is for log (X_1^2/2) for w_j \not =0 or \pi
	phi.prob.temp2 <-c(0.19,0.11,0.27,0.25,0.18)
	phi.mean.temp2 <- c(-2.20, -0.80, -0.55, -0.035, 0.48)
	phi.var.temp2 <- c(1.93, 1.01, 0.69, 0.60, 0.29)


	phi.mix <- rep(0,length(xi))
	phi.mean <- rep(0,length(xi))
	phi.var <- as.matrix(diag(0,length(xi)))
	for(i in 1:length(xi))
	{
		if(i ==1) { # this is for w_j=0
		#dnorm input is mean and sd NOT var
		log.alpha=array(0,5)
		for(j in 1:5){
			seednum=seednum+10;set.seed(seednum)
			log.alpha[j] <- (log(phi.prob.temp1[j]) 
                             + dnorm(xi[i], phi.mean.temp1[j], sqrt(phi.var.temp1[j]),log=TRUE))
		}
		log.alpha1= log.alpha-max(log.alpha)
		pi.hat = exp(log.alpha1)/sum(exp(log.alpha1))
#		sum.alpha <- alpha.1 + alpha.2 + alpha.3 + alpha.4 + alpha.5
#		pi.hat <- c(alpha.1/sum.alpha, alpha.2/sum.alpha, alpha.3/sum.alpha, alpha.4/sum.alpha, alpha.5/sum.alpha)
		seednum=seednum+10;set.seed(seednum)
		phi.mix[i] <- sample(x=(1:5), size=1, prob = pi.hat)
		phi.mean[i] <- phi.mean.temp1[phi.mix[i]]
		phi.var[i,i] <- phi.var.temp1[phi.mix[i]]
		}else{
		#dnorm input is mean and sd NOT var
		log.alpha=array(0,5)
		for(j in 1:5){
			seednum=seednum+10;set.seed(seednum)
			log.alpha[j] <- (log(phi.prob.temp2[j]) 
                             + dnorm(xi[i], phi.mean.temp2[j], sqrt(phi.var.temp2[j]),log=TRUE))
		}
		log.alpha1= log.alpha-max(log.alpha)
		pi.hat = exp(log.alpha1)/sum(exp(log.alpha1))
#		sum.alpha <- alpha.1 + alpha.2 + alpha.3 + alpha.4 + alpha.5
#		pi.hat <- c(alpha.1/sum.alpha, alpha.2/sum.alpha, alpha.3/sum.alpha, alpha.4/sum.alpha, alpha.5/sum.alpha)
		seednum=seednum+10;set.seed(seednum)
		phi.mix[i] <- sample(x=(1:5), size=1, prob = pi.hat)
		phi.mean[i] <- phi.mean.temp2[phi.mix[i]]
		phi.var[i,i] <- phi.var.temp2[phi.mix[i]]

		}
	}
	all$phi.mix=phi.mix
	all$phi.mean=phi.mean
	all$phi.var=phi.var
	all$seednum=seednum

	return(all)
}

UpdateEtaFunc = function(all,global.var){
# sub-routine for updating eta
# two options 1. update beta, 2. update eta directly

	require(MASS)
	seednum=all$seednum
	y=as.matrix(all$data[,1])
	n=length(y)
#old code for calculate gamma.inv.nn
#global.var$Q.mat=CompQFunc_old(dim(data0$data)[1]-kstep)
#gamma.inv.nn.old= (all$tau2.e)*global.var$Q.mat %*% diag(1/as.vector(all$lambda)) %*% t(global.var$Q.mat)
	cosQ0=cos(global.var$Q0.mat)
	sinQ0=sin(global.var$Q0.mat)
	gamma.inv.nn.1= cosQ0%*%diag(1/as.vector(all$lambda))%*%t(cosQ0)/n
	gamma.inv.nn.2= sinQ0%*%diag(1/as.vector(all$lambda))%*%t(sinQ0)/n
	gamma.inv.nn=(all$tau2.e)*(gamma.inv.nn.1+gamma.inv.nn.2)

	#option 1: update beta
	if(global.var$RegModel=='linear'){
		beta.mu0=global.var$beta.mu0
		beta.sigma20=global.var$beta.sigma20
		x=as.matrix(cbind(1,all$data[,-1]))
		

		beta.sigma2.star= chol2inv(chol(t(x)%*%gamma.inv.nn%*%x + 1/beta.sigma20))
		beta.mu.star=beta.sigma2.star%*%(t(x)%*%gamma.inv.nn%*%y
            	    + (1/beta.sigma20)*beta.mu0)

		#sample beta
		seednum=seednum+10;set.seed(seednum)
		beta.sample <- as.matrix(mvrnorm(1, beta.mu.star, beta.sigma2.star),length(beta.mu.star),1)

		all$beta=beta.sample
		all$eta=x%*%beta.sample	
		

	}

	#option 2: update beta for additive model
	if(global.var$RegModel=='additive'){
		beta.mu0=matrix(global.var$beta.mu0,dim(all$x)[2],1)
		beta.sigma20=global.var$beta.sigma20
		x=all$x
		

		beta.sigma2.star= chol2inv(chol(t(x)%*%gamma.inv.nn%*%x + 1/beta.sigma20))
		beta.mu.star=beta.sigma2.star%*%(t(x)%*%gamma.inv.nn%*%y
            	    + (1/beta.sigma20)*beta.mu0)

		#sample beta
		seednum=seednum+10;set.seed(seednum)
		beta.sample <- as.matrix(mvrnorm(1, beta.mu.star, beta.sigma2.star),length(beta.mu.star),1)

		all$beta=beta.sample
		all$eta=x%*%beta.sample	
	}

	#option 3 : nonparametric Eta
	if(global.var$RegModel=='nonparametric'){
		
		eta.mu0=all$eta.mu0  #eta prior mean
		eta.tau2=all$eta.tau2
		eta.rho=all$eta.rho
		dist0=as.matrix(dist(all$data.s[,-1]))
		eta.Sigma=exp(-eta.rho*dist0)/eta.tau2
		eta.Sigma.inv=chol2inv(chol(eta.Sigma))
		new.Sigma= chol2inv(chol(gamma.inv.nn + eta.Sigma.inv))
		new.mu= new.Sigma%*%gamma.inv.nn%*%(y-eta.mu0) + eta.mu0
		#generate Eta
		seednum=seednum+10;set.seed(seednum)
		all$eta = as.matrix(mvrnorm(1,new.mu,new.Sigma),length(new.mu),1)
	}



	all$epsi.0=y-all$eta
	all$epsi = (all$epsi.0)*sqrt(all$tau2.e)
      all$log.peri = CompLogPerioFunc(all$epsi)
	all$seednum=seednum
	return(all)
}


UpdateEtaPriorFunc = function(all,global.var){
# sub-routine for updating hyperparameters for Eta's Cov

	seednum=all$seednum
	eta=all$eta	
	eta.mu0=all$eta.mu0
	nn=length(eta)

	#update tau2 in Eta
	new.a=global.var$ab.eta.tau[1]+nn/2

	dist0=as.matrix(dist(all$data.s[,-1]))
	eta.Sigma0=exp(-all$eta.rho*dist0)
	eta.Sigma0.inv=chol2inv(chol(eta.Sigma0))
	new.b=(0.5*t(eta-eta.mu0)%*%eta.Sigma0.inv%*%(eta-eta.mu0) + global.var$ab.eta.tau[2])
	seednum=seednum+10;set.seed(seednum)
	all$eta.tau2=rgamma(1,shape=new.a,rate=new.b)		

	#update rho in Eta
	rho.vec=global.var$eta.rho.vec
	rho.dist.inv=all$eta.rho.dist.inv

	logprob=array(0,length(rho.vec))
	for(ri in 1:length(rho.vec)){
		logprob[ri]=-0.5*(all$eta.tau2)*t(eta-eta.mu0)%*%rho.dist.inv[,,ri]%*%(eta-eta.mu0)
		logprob[ri]=logprob[ri]+0.5*determinant(rho.dist.inv[,,ri])$mod
		#'+' since the matrix is already inversed
	}
	logprob1=logprob-max(logprob)
	prob.v=exp(logprob1)/sum(exp(logprob1))
	cumprob=cumsum(c(0,prob.v[1:(length(rho.vec)-1)]))
	seednum=seednum+10;set.seed(seednum)
	uni=runif(1)
	all$eta.rho=rho.vec[sum(cumprob <uni)]
	
	all$seednum=seednum
	return(all)
}


UpdateBetaFunc_Par = function(all,global.var){
# sub-routine for updating beta in the parametric model

	#option 1: update beta
	beta.mu0=global.var$beta.mu0
	beta.sigma20=global.var$beta.sigma20
	x=as.matrix(cbind(1,all$data[,-1]))
	seednum=all$seednum

	n=length(all$data[,1])
	ar.order=global.var$ar.order
	new.epsi=c(all$eps0,all$epsi) #eps -p+1, eps -p+2,.. eps0, eps1,..eps n
	#rewrite data vector with matrix form
	y.star=as.matrix(all$data[,1])
	for(i in 1:ar.order){
		y.star=y.star - all$ar.phi[i]*as.matrix(new.epsi[(ar.order-i+1):(ar.order-i+n)])
	}
	D.inv.nn=(1/all$sigma2)*diag(1,n)

	beta.sigma2.star= chol2inv(chol(t(x)%*%D.inv.nn%*%x + 1/beta.sigma20))
	beta.mu.star=beta.sigma2.star%*%(t(x)%*%D.inv.nn%*%y.star
                + (1/beta.sigma20)*beta.mu0)

	#sample beta
	require(MASS)
	seednum=seednum+10;set.seed(seednum)
	beta.sample <- as.matrix(mvrnorm(1, beta.mu.star, beta.sigma2.star),length(beta.mu.star),1)

	all$beta=beta.sample
	all$seednum=seednum
	all$epsi=as.matrix(all$data[,1])-x%*%beta.sample	

	return(all)
}

UpdateEpsFunc_Par=function(all,global.var){
#sub-routine to update eps 0, eps -1 ,,, for parametric error model
	ar.phi=all$ar.phi
	sigma2=all$sigma2
	beta=all$beta
	x=as.matrix(cbind(1,all$data[,-1]))
	seednum=all$seednum
	n=length(all$data[,1])
	ar.order=global.var$ar.order
	tau2=global.var$tau2
	y=as.matrix(all$data[,1])
	eps0=all$eps0
	eps=all$epsi
	if(ar.order==1){
		#eps 0
		sigma2.eps=1/(ar.phi[1]^2/sigma2 + 1/tau2)
		mu.eps=sigma2.eps*(1/sigma2)*ar.phi[1]*(y[1]-x[1,]%*%beta)
		seednum=seednum+10;set.seed(seednum)
		eps0[1]=rnorm(1,mu.eps, sqrt(sigma2.eps))
	}
	if(ar.order==2){
		#eps-1
		sigma2.eps=1/(ar.phi[2]^2/sigma2 + 1/tau2)
		mu.eps=sigma2.eps*(ar.phi[2]/sigma2)*(y[1]-x[1,]%*%beta -ar.phi[1]*eps0[2])
		seednum=seednum+10;set.seed(seednum)
		eps0[1]=rnorm(1,mu.eps, sqrt(sigma2.eps))
		#eps 0
		sigma2.eps=1/((ar.phi[1]^2+ar.phi[2]^2)/sigma2 + 1/tau2)
		mu.eps=sigma2.eps*(1/sigma2)*(
		ar.phi[1]*(y[1]-x[1,]%*%beta -ar.phi[2]*eps0[1])+
		ar.phi[2]*(y[2]-x[2,]%*%beta -ar.phi[1]*eps[1]))
		seednum=seednum+10;set.seed(seednum)
		eps0[2]=rnorm(1,mu.eps, sqrt(sigma2.eps))
	}
	all$eps0=eps0
	all$seednum=seednum

	return(all)
}

UpdatePhiFunc_Par=function(all,global.var){
#sub-routine to update phi(ar coef) for parametric error model

	beta=all$beta
	x=as.matrix(cbind(1,all$data[,-1]))
	seednum=all$seednum
	n=length(all$data[,1])
	ar.order=global.var$ar.order
	sigma2=all$sigma2
	new.epsi=c(all$eps0,all$epsi) #eps -p+1, eps -p+2,.. eps0, eps1,..eps n
	temp.phi=seq(-0.99,0.99,by=0.01)
	y=as.matrix(all$data[,1])
	
	for(i in 1:ar.order){
		logprob=c()
		ar.phi=all$ar.phi
		for(j in 1:length(temp.phi)){
			y.star=y
			ar.phi[i]=temp.phi[j]
			for(k in 1:ar.order){
				y.star=y.star - ar.phi[k]*as.matrix(new.epsi[(ar.order-k+1):(ar.order-k+n)])
			}
			logprob[j]= (-1/(2*sigma2))*t(y.star-x%*%beta)%*%(y.star-x%*%beta)			
		}
		prob=exp(logprob-max(logprob))/sum(exp(logprob-max(logprob)))
		notAR=TRUE
		while(notAR){
			ar.phi0=all$ar.phi
			seednum=seednum+10;set.seed(seednum)
			ar.phi0[i]=sample(temp.phi,size=1,prob=prob)
			if (sum(abs(polyroot(rev(c(1,-ar.phi0)))) <1)== ar.order){
				notAR=FALSE
			}
		}
		all$ar.phi[i]=ar.phi0[i]
	}
	all$seednum=seednum

	return(all)
}



UpdateSigma2Func_Par=function(all,global.var){
#sub-routine to update sigma2 for parametric error model

	beta=all$beta
	x=as.matrix(cbind(1,all$data[,-1]))
	seednum=all$seednum

	n=length(all$data[,1])
	ar.order=global.var$ar.order
	new.epsi=c(all$eps0,all$epsi) #eps -p+1, eps -p+2,.. eps0, eps1,..eps n
	#rewrite data vector with matrix form
	y.star=as.matrix(all$data[,1])
	for(i in 1:ar.order){
		y.star=y.star - all$ar.phi[i]*as.matrix(new.epsi[(ar.order-i+1):(ar.order-i+n)])
	}
	x=as.matrix(cbind(1,all$data[,-1]))
	
#	new.a=global.var$sigma2.hyp[1]+ n/2
#	new.b=1/(global.var$sigma2.hyp[2]+ (1/2)*t(y.star-x%*%beta)%*%(y.star-x%*%beta))
	new.a=global.var$sigma2.hyp[1]+ n/2
	new.b=(global.var$sigma2.hyp[2]+ (1/2)*t(y.star-x%*%beta)%*%(y.star-x%*%beta))
	
  seednum=seednum+10;set.seed(seednum)
	sigma2.inv=rgamma(1,shape=new.a, rate=new.b)
	all$sigma2=1/sigma2.inv
	all$seednum=seednum

	return(all)
	
}

UpdateTau2eFunc=function(all,global.var)
{
#sub-routine for updating variance of the error for nonparametric error model

	
	seednum=all$seednum
	#y=as.matrix(all$data[,1])
#	gamma.inv.nn0= global.var$Q.mat %*% diag(1/as.vector(all$lambda)) %*% t(global.var$Q.mat)
	cosQ0=cos(global.var$Q0.mat)
	sinQ0=sin(global.var$Q0.mat)
	gamma.inv.nn.1= cosQ0%*%diag(1/as.vector(all$lambda))%*%t(cosQ0)/n
	gamma.inv.nn.2= sinQ0%*%diag(1/as.vector(all$lambda))%*%t(sinQ0)/n
	gamma.inv.nn0=(gamma.inv.nn.1+gamma.inv.nn.2)


	epsi0=all$epsi.0 	#all$epsi.0=y-all$eta
	n.star=length(epsi0)

	new.a=global.var$ab.tau.e[1]+n.star/2
	new.b=(0.5*t(epsi0)%*%gamma.inv.nn0%*%epsi0+global.var$ab.tau.e[2])
	seednum=seednum+10;set.seed(seednum)
	all$tau2.e=rgamma(1,shape=new.a,rate=new.b)	
	all$epsi = (all$epsi.0)*sqrt(all$tau2.e)
      all$log.peri = CompLogPerioFunc(all$epsi)

	all$seednum=seednum

	return(all)
}



UpdateThetaFunc = function(all,global.var)
{
# sub-routine for updating theta

log.peri=all$log.peri
phi.mean=all$phi.mean
phi.var=all$phi.var
#prior for nu.n and tau.nn
nu.n=all$nu.n
tau.nn=all$tau.nn
seednum=all$seednum

	require(MASS)
	
	tau.n.inv=chol2inv(chol(tau.nn))
	phi.var.inv=chol2inv(chol(phi.var))

	#posterior var
	#tau.star <- as.matrix(tau.nn %*% chol2inv(chol(diag(dim(phi.var)[1]) + phi.var.inv%*% tau.nn)))

	tau.star <- chol2inv(chol(phi.var.inv + tau.n.inv))

	#posterior mean
	nu.star <- tau.star%*% phi.var.inv %*% (log.peri - phi.mean - nu.n) + nu.n


	###theta.sample <- rmvnorm(n=length(log.peri), mean= nu.star, sigma = tau.star)
	seednum=seednum+10;set.seed(seednum)
	theta.sample <- mvrnorm(1, nu.star, tau.star)

	all$theta=as.matrix(theta.sample)
	all$lambda = ConvThetaToLambdaFunc(dim(all$data)[1],all$theta)
	#after rescaling lambda, update theta accordingly
	all$theta= ConvLambdaToThetaFunc(length(theta.sample),all$lambda)
	all$seednum=seednum
	

	return(all)
}


UpdateThetaPriorFunc=function(all,global.var)
{
#sub-routine for updating tau2 and rho, hyperparametrs for Theta

	nn=length(all$epsi)
	theta=as.matrix(all$theta)
	nu=as.matrix(all$nu.n)
	tmp.theta.pr=CompThetaPriorFunc(nn,1,all$rho)
	tau0.nn=tmp.theta.pr$tau.nn
	n.star=length(all$theta)
	
	#update tau2
	seednum=all$seednum
	new.a=global.var$ab.tau[1]+n.star/2
	new.b=(0.5*t(theta-nu)%*%chol2inv(chol(tau0.nn))%*%(theta-nu)+global.var$ab.tau[2])
	seednum=seednum+10;set.seed(seednum)
	all$tau2=rgamma(1,shape=new.a,rate=new.b)	
	
	#update rho
	rho.vec=global.var$rho.vec
	rho.dist.inv=all$rho.dist.inv

	
	logprob=array(0,length(rho.vec))
	for(ri in 1:length(rho.vec)){
		logprob[ri]=-0.5*(all$tau2)*t(theta-nu)%*%rho.dist.inv[,,ri]%*%(theta-nu)
		logprob[ri]=logprob[ri]+0.5*determinant(rho.dist.inv[,,ri])$mod
		#'+' since the matrix is already inversed
	}
	logprob1=logprob-max(logprob)
	prob.v=exp(logprob1)/sum(exp(logprob1))
	cumprob=cumsum(c(0,prob.v[1:(length(rho.vec)-1)]))
	seednum=seednum+10;set.seed(seednum)
	uni=runif(1)
	all$rho=rho.vec[sum(cumprob <uni)]
	theta.pr=CompThetaPriorFunc(nn,all$tau2,all$rho)
	#all$nu.n=theta.pr$nu.n
	all$tau.nn=theta.pr$tau.nn
	all$seednum=seednum

	return(all)

}



