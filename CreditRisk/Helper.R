rm(list=ls())
##############################################################
##############################################################
Black.Scholes<-function(type,S,K,r,div,vol,T){
  # type: +1 for call; all other values returns the put
  # S: initial asset value 
  # K: strike price
  # r: continuously compounded yearly risk-free rate
  # div: continuously compounded yearly dividend rate
  # vol: annualized standard deviation of log return
  # T: option maturity
  
  dp<-(log(S/K)+(r-div+0.5*vol^2)*T)/(vol*sqrt(T))
  dm<-dp-vol*sqrt(T)
  if(type==1){
    f<-S*exp(-div*T)*pnorm(dp)-K*exp(-r*T)*pnorm(dm)
  }
  else{
    f<-K*exp(-r*T)*pnorm(-dm)-S*exp(-div*T)*pnorm(-dp)
  }  
}

##############################################################
##############################################################
Delta<-function(type,K,S,r,div,vol,T){
  # type: +1 for call; all other values returns the put
  # K: strike price
  # S: initial asset value 
  # r: continuously compounded yearly risk-free rate
  # div: continuously compounded yearly dividend rate
  # vol: annualized standard deviation of log return
  # T: option maturity
  d<-(log(S/K)+(r-div+0.5*vol^2)*T)/(vol*sqrt(T))
  if(type==1){
    f<-pnorm(d)
  }else{
    f<-pnorm(d)-1
  }
  f<-f*exp(-div*T)
}

##############################################################
##############################################################
Gamma<-function(K,S,r,div,vol,T){
  # K: strike price
  # S: initial asset value 
  # r: continuously compounded yearly risk-free rate
  # div: continuously compounded yearly dividend rate
  # vol: annualized standard deviation of log return
  # T: option maturity
  d<-(log(S/K)+(r-div+0.5*vol^2)*T)/(vol*sqrt(T))
  f<-exp(-div*T)*dnorm(d)/(vol*S*sqrt(T))
}

##############################################################
##############################################################
Theta<-function(type,K,S,r,div,vol,T){
  # type: +1 for call; all other values returns the put
  # K: strike price
  # S: initial asset value 
  # r: continuously compounded yearly risk-free rate
  # div: continuously compounded yearly dividend rate
  # vol: annualized standard deviation of log return
  # T: option maturity
  dp<-(log(S/K)+(r-div+0.5*vol^2)*T)/(vol*sqrt(T))
  dm<-dp-vol*sqrt(T)
  if(type==1){
    f<--exp(-div*T)*vol*S*dnorm(dp)/(2*sqrt(T))-r*K*exp(-r*T)*pnorm(dm)+div*S*exp(-div*T)*pnorm(dp)
  }else{
    f<--exp(-div*T)*vol*S*dnorm(dp)/(2*sqrt(T))+r*K*exp(-r*T)*pnorm(-dm)-div*S*exp(-div*T)*pnorm(-dp)
  }
  f
}


##############################################################
##############################################################
Loss.Quantile<-function(y,Q.CDF, call,put,N.sim,dt,alpha){
  #returns the probability of a loss greater than x, as well as the quantile associated with alpha
  #Portfolio of European calls and puts
  
  #call/put: inputs associated with the call and put options in the portfolio
  #N.sim: number of simulations
  #rho: correlation matrix of the BM driving the underlying asset model
  #dt: horizon
  #N: number of steps to evaluate characteristic function
  #h: step size
  #alpha: tail probability
  
  
  # YOU NEED TO GENERATE THE THETA, DELTA, AND GAMMA OF THE PORTFOLIO 
  # AS WELL AS WRITE A FUNCTION THAT VALUES THE PORTFOLIO

  ##############################################################
  #matrix to save the simulated exact loss, approximate loss, probabilities, and weights
  #column 1: exact loss
  #column 2: delta-gamma approximate loss
  #column 3: complementary CDF
  #column 4: weights
  L.mat<-matrix(0,N.sim,4)
  #run the simulations in a loop
  for(i in 1:N.sim){
    #normally-distributed risk factor vector
    z<-rnorm(m)
    temp<-vol*sqrt(dt)*z
    dS<-S0*temp
    #delta-gamma approximate loss
    L.mat[i,2]<-a+b%*%z+diag(Lambda)%*%(z^2)
    #update underlying values at time dt
    call$S<-S0*(1+temp)
    put$S<-call$S
    #compute exact loss function given the outer simulation scenario
    L.mat[i,1]<-V.init-Portfolio(call,put)
  }
  #sort matrix according to increasing order of exact loss function
  L.mat<-L.mat[order(L.mat[,1],decreasing =FALSE),]
  #determine if delta-gamma loss function exceeds some threshold x specified by user
  L.mat[,3]<-as.numeric(L.mat[,2]>y)
  #sum up exceedances for determining weights
  L.sum.1<-sum(L.mat[,3])
  L.sum.2<-N.sim-L.sum.1
  #determine weights
  w.1<-(1-Q.CDF)/L.sum.1
  w.2<-Q.CDF/L.sum.2
  #   L.mat[,4]<-ifelse(L.mat[,3]>0,w.1,w.2)
  ind.temp<-as.numeric(L.mat[,3])
  L.mat[,4]<-w.1*ind.temp+w.2*(1-ind.temp)
  #sum(weights[j:N.sim]), j=1,2,...,N.sim
  #this is the empirical complementary CDF
  L.mat[,3]<-sapply(X=1:N.sim,FUN=function(x) sum(L.mat[x:N.sim,4]))
  
  #return index of minimum loss >= to the quantile
  i.1<-max(which(L.mat[,3]>alpha))
  i.2<-i.1+1
  #interpolation factor
  a<-(alpha-L.mat[i.2,3])/(L.mat[i.1,3]-L.mat[i.2,3])
  #loss quantile estimate
  x.alpha<-a*L.mat[i.1,1]+(1-a)*L.mat[i.2,1]
  

}

##############################################################
##############################################################
Loss.Probability.IS<-function(y,v,call,put,N.sim,rho,dt,theta.0,tol){
  #returns the probability of a loss greater than x 
  #Portfolio of calls and puts
  
  #x: loss threshold for control variate
  #call/put: inputs associated with the call and put options in the portfolio
  #N.sim: number of simulations
  #rho: correlation matrix of the BM driving the underlying asset model
  #dt: horizon
  #theta.0: starting point for Newton-Raphson algorithm
  #tol: tolerance determining when N-R algorithm stops
  
  
  # YOU NEED TO GENERATE THE THETA, DELTA, AND GAMMA OF THE PORTFOLIO 
  # AS WELL AS WRITE A FUNCTION THAT VALUES THE PORTFOLIO
 
  ##############################################################
  #find IS value of theta
  x = y - a
  theta.max<-max(1/(2*diag(Lambda)))
  phi<-0.97
  
  theta<-uniroot(theta.root, c(tol,phi*theta.max), a,b,v,diag(Lambda),x,tol = tol, maxiter = 1000)$root
  #shifted mean and variance-covariance matrix under P(theta)
  psi<-Q.cumulant(a,b,v,diag(Lambda),theta)
  ##############################################################
  #run simulations under P(theta)
  Y.shape<-v/2
  alpha.theta = - (theta * x/v) + .5 * sum(theta^2 * (b^2/v)/ (1-2*theta*diag(Lambda)))
  Y.scale = 2/(1-2*alpha.theta)
  Y = rgamma(N.sim, shape=Y.shape, scale=Y.scale)
  LR<-rep(0,N.sim)
  V.f<-rep(0,N.sim)
  counter = 0
  for(i in 1:N.sim){
    
    mu.theta<-theta*b*sqrt(Y[i]/v)/(1-2*diag(Lambda)*theta)
    sigma.theta<-sqrt(1/(1-2*diag(Lambda)*theta))
    
    Z<-t(mu.theta+sigma.theta*rnorm(m))
    X = Z/sqrt(Y[i]/v)
    dS<-C %*% X
    Q<- b %*%X + t(X) %*% Lambda %*% X
    #Q<- -deltaV %*% dS -.5 * t(dS) %*% gammaV %*% dS
    Q.x = (Y[i]/v)*(Q-x)
    LR[i]<-exp(-theta*Q.x+psi)
    call$S<-S0+dS
    put$S<-call$S
    if(length(which(call$S < 0))){
      cat("neg S = ", dS, "\n")
      cat("X = ", X, ", Y = " ,Y[i], "\n")
      counter = counter + 1
      next
    }
    call$T = T-dt
    put$T = call$T
 #   call$vol = vol*call$S*sqrt(dt)
 #   put$vol = call$vol
    #compute the portfolio value given the outer simulation scenario
    V.f[i]<-Portfolio(call,put)
  }
  cat("counter = ", counter, "\n")
  #generate indicator function
  V.f<-LR*as.numeric(V.init-V.f>y)
  

}

##############################################################
##############################################################
inversion.CDF<-function(CF,x,h,N){
  # CF: characteristic function of the random variable
  # x: point at which the CDF is evaluated at
  # h: step size for discretizing the inversion integral
  # N: cutoff (integration limits are plus/minus (N-1)*h)
  
  i<-sqrt(as.complex(-1))
  a<-0
  
  for(k in 1:(N-1)){
    CF1<-CF[k]
    temp<-CF1*exp(-i*h*k*x)/(2*pi*i*k)
    a<-a+temp
    CF1<-CF[k+N]
    temp<--CF1*exp(i*h*k*x)/(2*pi*i*k)
    a<-a+temp
  }
  
  f<-0.5+(h*x/(2*pi))-a
  Re(f)
}

##############################################################
##############################################################
Q.CF<-function(a,b,lambda,x){
  #imaginary number "i"
  i<-sqrt(as.complex(-1))
  
  temp<-x^2*b^2/(1-2*i*x*lambda)+log(1-2*i*x*lambda)
  
  f<-exp(i*a*x-0.5*sum(temp))
}

##############################################################
##############################################################
Q.cumulant<-function(a,b,v,lambda,theta){  
  alpha.theta = - (theta * x/v) + .5 * sum(theta^2 * (b^2/v)/ (1-2*theta*lambda))
  temp = 1/sqrt(1-2*theta*lambda)
  f<-(1-2*alpha.theta)^(-v/2) * prod(temp)
  f<-log(f)
  
}

##############################################################
##############################################################
theta.root<-function(theta,a,b,v,lambda,x){
  temp3<-lambda/(1-2*theta*lambda)
  sum_temp3<-sum(temp3)
  dalpha.dtheta = (2*x/v) - (2/v) * sum(theta * b^2 * (1+lambda * theta) / (1-2*theta*lambda))
  alpha.theta = - (theta * x/v) + .5 * sum(theta^2 * (b^2/v)/ (1-2*theta*lambda))
  f = sum_temp3 - (v/2) * (dalpha.dtheta/(1+alpha.theta))
#   temp1<-(theta*b/v)^2 * (1/(1-2*theta*lambda))
#   sum_temp1 = sum(temp1)
#   temp2<-b^2/(1-2*theta*lambda)^2
#   sum_temp2<-sum(temp2)
#   
#   f<- -(v/2) * (1/(1+(2*theta*x/v)-sum_temp1)) * ((2*x/v)-(2*theta/v)*sum_temp2) + sum_temp3
}

##############################################################
