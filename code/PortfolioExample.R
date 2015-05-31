#Load example portfolio 

# This is a simple portfolio designed to illustrate the concept of manager component VaR
Portfolio <- read.csv("../data/RiskDecompositionExample.csv", header= TRUE, sep = ",", quote="\"");

# The following files have also been generated since there is no R code to map the instruments to risk factors and calculate the greeks
# The following files have been exported from HedgeFacts Analytics

# w_delta[i,j] is the contribution from instrument i to j^{th} element of the delta vector 
w_delta <- as.matrix(read.csv("../data/w_delta.csv", header= F)) 
# w_gamma[i,j*N + k] is the contribution from instrument i to j,k^{th} element of the gamma matrix
w_gamma <- as.matrix(read.csv("../data/w_gamma.csv", header= F))
# This is the Delta vector - each element represents a risk factor
Delta <- as.matrix(read.csv("../data/Delta.csv", header= F))
# This is the gamma matrix - each row and column correspond to risk factors
Gamma <- as.matrix(read.csv("../data/Gamma.csv", header= F))
# This is the covariance matrix measured from historical returns of risk factors
Sigma <- as.matrix(read.csv("../data/Sigma.csv", header= F)) 

N<- dim(w_delta)[2] # number of risk factors
M<- dim(w_delta)[1] # number of instruments