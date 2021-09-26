#### Example 1: A portfolio with a stock and a call option on another stock. 
#### The two stocks are uncorrelated. This example compares different methodologies for estimating 
#### the component VaR of each instrument. 


require("greeks")

kernel<-function(x,h){
  
  ret<-1-abs(x/h)
  idx<-which(ret<0)
  ret[idx]<-0
  return(ret)
}

BlackScholes <- function(S, K, r, T, sig, type){
  
  if(type=="C"){
    d1 <- (log(S/K) + (r + sig^2/2)*T) / (sig*sqrt(T))
    d2 <- d1 - sig*sqrt(T)
    
    value <- S*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
    return(value)}
  
  if(type=="P"){
    d1 <- (log(S/K) + (r + sig^2/2)*T) / (sig*sqrt(T))
    d2 <- d1 - sig*sqrt(T)
    
    value <-  (K*exp(-r*T)*pnorm(-d2) - S*pnorm(-d1))
    return(value)}
}
## Initialize paramters
sigma_1<-0.1 # volatility of returns of S_1
sigma_2<-0.1 # volatility of returns of S_2
dt <- 1/252 # time step size for risk analysis (one day VaR)
T<-7/252 # maturity of option
rho_12<-0.0 # correlation between returns of S_1 and S_2
c<-0.01 # confidence
r<-0.01 # instantaneous annualized short rate
mu_1<-0.0 # drift of returns of S_1
mu_2<-0.0 # drift of returns of S_2
S_1_0<-100 # Spot price of S_1
S_2_0<-100 # Spot price of S_2
MM<-1e7 #Number of samples
# Construct covariance matrix
Sigma<-matrix(c(sigma_1**2*dt, rho_12*sigma_1*sigma_2*dt, rho_12*sigma_1*sigma_2*dt, sigma_2**2*dt), nrow=2, ncol=2)

CVaR_1s.kernel<-rep(0,log10(MM)-3+1)
CVaR_2s.kernel<-rep(0,log10(MM)-3+1)
VaRs.kernel<-rep(0,log10(MM)-3+1)
# Perform Monte-Carlo simulation and using the Kernel method to estimate the Component VaR.
i<-1
for (M in 10**seq(3, log10(MM))){
  Z<-matrix(rep(0,2*M),nrow=2, ncol=M)
  Z[,(M/2+1):M]<-matrix(rnorm(M), nrow=2, ncol=as.integer(M/2))
  Z[,1:(M/2)]<- - Z[,(M/2+1):M]
  C<-chol(Sigma)
  R<-C%*%Z
  R[1,]<-R[1,]+rep(mu_1*dt,M)
  R[2,]<-R[2,]+rep(mu_2*dt,M)
  
  l_p_1<- -R[1,]
  l_p_2<- -R[2,]
  
  Rw <- R
  Rw[1,]<-Rw[1,]*S_1_0
  Rw[2,]<-Rw[2,]*S_2_0
  
  # P=S_1 + V(S_2), where V is call option on stock 2
  
  S_2_0<-100
  K<-100
  # Perform MC sampling with full re-pricing
  V_0<-BS_European_Greeks(initial_price = S_2_0, exercise_price = K,
                     r=0.01, time_to_maturity = T, dividend_yield = 0.0, volatility = sigma_2,
                     greek = c("fair_value"), payoff = "call")
  #V_0<-BlackScholes(S_2_0, K, r, 1, sigma_2, 'C')
  S_2<-S_2_0*(1-l_p_2)
  
  V_1<-BlackScholes(S_2, K, r, T-dt, sigma_2, 'C') # allows vector input
  l_p_2<-(V_0-V_1)
  l_p<-S_1_0*l_p_1+l_p_2
  #l_p_<-S_1_0*l_p_1+S_2_0*res[1]*l_p_2
  l_p_s<-sort(l_p)
  VaR.kernel<--l_p_s[as.integer((1-c)*M)]
  h<-sd(l_p)*2.575*M**(-1/5)
  
  CVaR_1s.kernel[i]<-VaR.kernel*sum(kernel(-l_p-VaR.kernel,h)*-l_p_1*S_1_0)/sum(kernel(-l_p-VaR.kernel,h)*-l_p)
  CVaR_2s.kernel[i]<-VaR.kernel*sum(kernel(-l_p-VaR.kernel,h)*-l_p_2)/sum(kernel(-l_p-VaR.kernel,h)*-l_p)
  VaRs.kernel[i]<-VaR.kernel
  i<-i+1
}
## linear and DG CVaR calculations


# First compute the DG CVaR with CF expansions
res<-BS_European_Greeks(initial_price = S_2_0, exercise_price = K,
                        r, time_to_maturity = T, dividend_yield = 0.0, volatility = sigma_2,
                        greek = c("delta","gamma"), payoff = "call")



w_delta<-matrix(c(S_1_0,0,0,res[1]*S_2_0), nrow=2, ncol=2)
w_gamma<-matrix(c(0,0,0,0,0,0,0,res[2]*S_2_0**2), nrow=2)


Delta<-matrix(c(S_1_0,res[1]*S_2_0), nrow=2, ncol=1)
Gamma<-matrix(c(0,0,0,res[2]*S_2_0**2), nrow=2, ncol=2)
N<- dim(w_delta)[2] # number of risk factors
M<- dim(w_delta)[1] # number of instruments
source("RiskDecomposition.R")

w<-matrix(c(S_1_0,res[1]*S_2_0), nrow=1, ncol=2)
mu<-matrix(c(r*dt,r*dt), nrow=1, ncol=2)

### Sanity check (delta-gamma with gaussian loss distribution: P=S+V)
### Delta-Gamma CVaR (with gaussian loss distribution)
mu_2_x<- 0.5*(S_2_0**2)*res[2]*(sigma_2**2*dt+ mu[2]**2)
mu_dp<-w%*%t(mu) +mu_2_x
sigma_p_2<-(S_1_0)**2*sigma_1**2*dt + (S_2_0*res[1])**2*sigma_2**2*dt + 0.5*(res[2]*S_2_0**2)**2*sigma_2**4*dt**2
sigma_p_1_2<-(S_1_0)**2*sigma_1**2*dt/sqrt(sigma_p_2)
sigma_p_2_2<-((S_2_0*res[1])**2*sigma_2**2*dt + 0.5*(res[2]*S_2_0**2)**2*sigma_2**4*dt**2)/sqrt(sigma_p_2)

CVaR_1.dg<-w[1]*mu[1]- qnorm(1-c)*sigma_p_1_2
CVaR_2.dg<-w[2]*mu[2] + mu_2_x - qnorm(1-c)*sigma_p_2_2
VaR.dg<-mu_dp - qnorm(1-c)*sqrt(sigma_p_2)


### Delta CVaR
mu_dp<-w%*%t(mu)
sigma_p_2<-(S_1_0)**2*sigma_1**2*dt + (S_2_0*res[1])**2*sigma_2**2*dt
sigma_p_1_2<-(S_1_0)**2*sigma_1**2*dt/sqrt(sigma_p_2)
sigma_p_2_2<-((S_2_0*res[1])**2*sigma_2**2*dt)/sqrt(sigma_p_2)

CVaR_1.linear<-w[1]*mu[1]- qnorm(1-c)*sigma_p_1_2
CVaR_2.linear<-w[2]*mu[2]- qnorm(1-c)*sigma_p_2_2
VaR.linear<-mu_dp - qnorm(1-c)*sqrt(sigma_p_2)

precision<-5
print(round(VaR.linear,precision))
print(c(round(CVaR_1.linear,precision),round(CVaR_2.linear,precision)))

print(round(VaR.dg,precision))
print(c(round(CVaR_1.dg,precision),round(CVaR_2.dg,precision) ))
print(round(VaR,precision))
print(round(instrument.VaR,precision))
print(round(VaRs.kernel[length(VaRs.kernel)],precision))
print(c(round(CVaR_1s.kernel[length(CVaR_2s.kernel)],precision),round(CVaR_2s.kernel[length(CVaR_2s.kernel)],precision) ))


plot(seq(3, log10(MM)),VaRs.kernel,type='l', col='black', ylim=c(-1.66,-0.1), xlab='log_10 #simulations', ylab='VaR', cex=1.2)
lines(seq(3, log10(MM)),CVaR_1s.kernel, col='red')
lines(seq(3, log10(MM)),CVaR_2s.kernel, col='green')
lines(seq(3, log10(MM)), rep(VaR,5), col='black', lty='dashed')
lines(seq(3, log10(MM)), rep(instrument.VaR[1],5), col='red', lty='dashed')
lines(seq(3, log10(MM)), rep(instrument.VaR[2],5), col='green', lty='dashed')
legend(5,-0.42, legend=c("MC VaR", "Kernel CVaR 1", "Kernel CVaR 2","DG VaR", "DG CVaR 1", "DG CVaR 2"),
       col=c("black", "red","green","black", "red","green"), lty=c('solid','solid','solid','dashed','dashed','dashed'), cex=1.0,bg='lightblue')


