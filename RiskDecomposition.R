library("matrixcalc")

# calculate moments of the distribution
# inputs: 
# Delta is a vector of Nx1
# Gamma is a N x N matrix
# Sigma is a N x N matrix
# outputs: an array of the first four moments
calcMoments<-function(Delta, Gamma, Sigma){
  mu <- array()
  GS <- Gamma%*%Sigma
  mu[1] <- 0.5*matrix.trace(GS)
  mu[2] <- t(Delta)%*%Sigma%*%Delta + 0.5*matrix.trace(GS)^2
  mu[3] <- 3.0*t(Delta)%*%Sigma%*%GS%*%Delta + matrix.trace(GS)^3
  mu[4] <- 12.0*t(Delta)%*%Sigma%*%GS%*%GS%*%Delta + 3.0*matrix.trace(GS)^4 + 3.0*mu[2]^2
  return(mu)
}

#Calculate the instrument component of the mean of portfolio losses
#inputs: instrument object
calcComponentMean<-function(idx){
  # collect all risk factors associated with the instrument
  mu <- 0.5*matrix.trace(matrix(w_gamma[idx,],N,N)%*%Sigma)
  return(mu);
}

#Calculate the instrument component of the variance of portfolio losses
#inputs: instrument object
calcComponentVariance<-function(idx){
    
    # collect all risk factors associated with the instrument
    # Convexity<- 0.5*matrix.trace(matrix(w_gamma[idx,],N,N)%*%Sigma)*matrix.trace(Gamma%*%Sigma)
    Convexity<- 0.5*matrix.trace(matrix(w_gamma[idx,],N,N)%*%nabla_gamma)
    #print(Convexity)
    sigma <- Convexity;
    for (k in 1:N) {
       sigma <- sigma + 0.5*w_delta[idx,k]*nabla_delta[k]
    }
    return(sigma);
}

z <- qnorm(1-c)
mu <- calcMoments(Delta, Gamma, Sigma)
k <- mu[4]/mu[2]^2 # kurtosis
s <- mu[3]/(mu[2]^1.5) # skew

# calculate the variance of portfolio losses under the parametric loss distribution
nabla_delta <- 2.0*Sigma %*% Delta
nabla_gamma <- matrix.trace(Sigma%*%Gamma)*Sigma
sigma_P <- sqrt(0.5*t(Delta) %*% nabla_delta + 0.5*matrix.trace(Gamma%*%nabla_gamma))

instrument.mu    <- matrix(, nrow = M, ncol = 1)
instrument.sigma2 <- matrix(, nrow = M, ncol = 1)
instrument.VaR <- matrix(, nrow = M, ncol = 1)

for (idx in 1:M){
  instrument.mu[idx]     <-calcComponentMean(idx)  
  instrument.sigma2[idx] <-calcComponentVariance(idx)
  # Calculate instrument component VaR
  instrument.VaR[idx]    <- (instrument.mu[idx] - (z+ 1/6*(z^2-1)*s + 1/24*(z^3-3*z)*(k-3) - 1/36*(2*z^3-5*z)*s^2)*instrument.sigma2[idx]/sigma_P)
}
# Calculate total VaR
VaR <- (mu[1] - (z+ 1/6*(z^2-1)*s + 1/24*(z^3-3*z)*(k-3) - 1/36*(2*z^3-5*z)*s^2)*sigma_P) 
