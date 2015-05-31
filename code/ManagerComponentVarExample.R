source('RiskDecomposition.R')

# This file presents a simple example of how the manager component VaR is calculated

# Arbitrarily assign a number of manager and instruments to those managers for this example
managers <- list()
managers[[1]] <-seq(1,4)
managers[[2]] <-seq(5,7)
managers[[3]] <-seq(8,11)
manager.VaR <- matrix( , nrow = 3, ncol = 1)

for (i in 1:3){
  manager.mu <-0
  manager.sigma2 <-0
  
  for (j in managers[[i]]){
    manager.mu <- manager.mu + com_mu[j]  
    manager.sigma2 <- manager.sigma2 + com_sigma2[j]
  }
  # Calculcate the manager compomnent VaR
  manager.VaR[i] <- -(manager.mu + (z+ 1/6*(z^2-1)*s + 1/24*(z^3-3*z)*(k-3) - 1/36*(2*z^3-5*z)*s^2)*manager.sigma2/sigma_P)
}