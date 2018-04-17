# +++++++++++++++++++++++++++++++
# + Name: M Masoud +
# +++++++++++++++++++++++++++++++

rm(list=ls())
source("Helper.R")
library(sciplot)

S0 = 100
r = .05
div = 0
vol = .3
T = .5
K=100

n.calls = 10
n.puts = 5
n.assets = 10


CBS = Black.Scholes(type=1,S=S0,K=K,r=r,div=div,vol=vol,T=T)
PBS = Black.Scholes(type=0,S=S0,K=K,r=r,div=div,vol=vol,T=T)
V0.exact = n.assets * (-n.calls * CBS - n.puts *  PBS)

Portfolio<-function(call,put){
  Sc = call$S
  TTMc = call$T
  Kc=call$K
  rc=call$r
  divc=call$div
  volc = call$vol
  
  Sp = put$S
  TTMp = put$T
  Kp=put$K
  rp=put$r
  divp=put$div
  volp = put$vol
  
  
  CBS = Black.Scholes(type=1,S=Sc,K=Kc,r=rc,div=divc,vol=volc,T=TTMc)
  PBS = Black.Scholes(type=0,S=Sp,K=Kp,r=rp,div=divp,vol=volp,T=TTMp)
  V.exact = -n.calls * sum(CBS) - n.puts *  sum(PBS)
  V.exact
}

Portfolio.Theta<-function(call,put){
  Sc = call$S
  TTMc = call$T
  Kc=call$K
  rc=call$r
  divc=call$div
  volc = call$vol
  
  Sp = put$S
  TTMp = put$T
  Kp=put$K
  rp=put$r
  divp=put$div
  volp = put$vol
  
  
  thetaC = Theta(type=1,K=Kc,S=Sc,r=rc,div=divc,vol=volc,T=TTMc)
  thetaP = Theta(type=0,K=Kp,S=Sp,r=rp,div=divp,vol=volp,T=TTMp)
  
  thetaV = -n.calls * sum(thetaC) - n.puts * sum(thetaP)
  thetaV  
}

Portfolio.Delta<-function(call,put){
  Sc = call$S
  TTMc = call$T
  Kc=call$K
  rc=call$r
  divc=call$div
  volc = call$vol
  
  Sp = put$S
  TTMp = put$T
  Kp=put$K
  rp=put$r
  divp=put$div
  volp = put$vol
  
  deltaC = Delta(type=1,S=Sc,K=Kc,r=rc,div=divc,vol=volc,T=TTMc)
  deltaP = Delta(type=0,S=Sp,K=Kp,r=rp,div=divp,vol=volp,T=TTMp)
  
  deltaV = -n.calls* deltaC - n.puts*deltaP
  deltaV  
}

Portfolio.Gamma<-function(call,put){
  S=call$S
  K=call$K
  r=call$r
  div=call$div
  vol=call$vol
  T=call$T
  
  gammaC = Gamma(K=K,S=S,r=r,div=div,vol=vol,T=T)
  gammaP = gammaC
  gammaV = -n.calls*gammaC - n.puts*gammaP
  gammaV = diag(gammaV)
}

cat("##### Q1 ##### \n")

### exact loss simulations
x=311

set.seed(123)
N = 40000
dt = .04
L.p.exact = rep(0,N)
m = n.assets
call = list(S=S0,K=K,S=S0,r=r,div=div,vol=vol,T=T)
put = list(S=S0,K=K,S=S0,r=r,div=div,vol=vol,T=T)

for(i in 1:N){
  z<-rnorm(m)
  S<-S0 * (1 + vol*sqrt(dt)*z)
  call$S = S
  put$S = S
  call$T=T-dt
  put$T=call$T
  L.p.exact[i] = as.numeric(V0.exact-Portfolio(call,put) > x)
}

L.se.exact = se(L.p.exact)
L.p.exact = mean(L.p.exact)

cat("P(L>x) exact = ",L.p.exact,"\n")
cat("P(L>x) exact standard error = ",L.se.exact,"\n")

v=5

sigma.S = diag(rep(vol^2*S0^2*dt,m))
call$S = rep(S0, m)
call$T = T
#call$vol = vol*S0*sqrt(dt)
put = call
gammaV = Portfolio.Gamma(call,put)

C.tilde<-t(chol(sigma.S))
A<--0.5*t(C.tilde)%*%(gammaV)%*%C.tilde
temp<-eigen(A)
U<-temp$vectors 
Lambda<-matrix(0,m,m) 
diag(Lambda)<-temp$values 
C<-C.tilde%*%U


# a, b
thetaV = Portfolio.Theta(call,put)
a = - thetaV * dt

deltaV = Portfolio.Delta(call,put)

b = - t(deltaV) %*% C

call$vol = vol*S0*sqrt(dt)
put$vol = call$vol
V.init<-Portfolio(call,put)



### Important Sampling

max.theta = 1 / (2*max(diag(Lambda))) - .0001
v=5
set.seed(123)
Loss.IS.ind = Loss.Probability.IS(y=x,v=v,call=call,put=put,N.sim=N,rho=NA,dt=dt,theta.0=max.theta,tol=1e-5)

L.p.IS = mean(Loss.IS.ind)
L.se.IS = se(Loss.IS.ind)

cat("P(L>x) Importance Sampling = ",L.p.IS,"\n")
cat("P(L>x) Importance standard error = ",L.se.IS,"\n")


# Q1.res.table = matrix(rep(0,9), ncol=3)
# colnames(Q1.res.table) = c("Loss probability", "Standard error", "VRR")
# rownames(Q1.res.table) = c("Exact Loss", "Control variate", "Importance sampling")
# 
# Q1.res.table[1,] = c(L.p.exact, L.se.exact, 1)
# Q1.res.table[2,] = c(L.p.CV, L.se.CV, L.se.exact/L.se.CV)
# Q1.res.table[3,] = c(L.p.IS, L.se.IS, L.se.exact/L.se.IS)
# 
# write.table(Q1.res.table, file = "Q1.res.csv", sep = ",",qmethod = "double")
# 
# 
# 
