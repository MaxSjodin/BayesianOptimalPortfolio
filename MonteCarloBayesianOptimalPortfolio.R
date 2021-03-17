library(matlab)
library(dplyr)
library(MASS)

meanVarianceWeights <- function(B = 10000, k = 40, n = 50, alpha = 50, meanlow = -0.01, meanhigh = 0.01, 
                       volatilitylow = 0.002, volatilityhigh = 0.005, corr = 0.6,
                       r_0 = 100, d_0 = 100){

  ones <- t(t(rep(1, k)))
  samples <- c()
  for(i in 1:B){
    # Parameters
    mu <- runif(k, meanlow, meanhigh)
    volatility <- runif(k, volatilitylow, volatilityhigh)
    correlation <- (1-corr)*diag(k)+corr*ones(k,k)
    S_mat <- diag(volatility)%*%correlation%*%diag(volatility)
    # Sample from multivariate normal distribution
    samples <- rbind(samples,mvrnorm(n=1, mu = mu, Sigma = S_mat))
  }
  
  mean_ret <- colMeans(samples)
  
  S_mat <- t(samples - mean_ret)%*%(samples - mean_ret)

  m_0 <- mean_ret + 0.5*runif(k, meanlow, meanhigh)
  S_0 <- S_mat + 0.5*diag(runif(k, volatilitylow, volatilityhigh))
  
  ### Frequentist
  d_n <- 1/(n-1)
  
  # Weights for global minimum variance portfolio
  wts_gmv <- (solve(S_mat)%*%ones)/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 
  
  # Q-matrix used to solve the optimization problem of the quadratic utility funciton
  Q <- solve(S_mat) - (solve(S_mat)%*%ones%*%t(ones)%*%solve(S_mat))/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 

  ### Deriving quantities for diffuse prior
  c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))

  ### Deriving quantities for conjugate prior
  mean_ret_c <- (n*mean_ret+r_0*m_0)/(n+r_0)
  
  S_mat_c <- S_mat+S_0+n*r_0*(m_0-mean_ret_c)%*%t(m_0-mean_ret_c)/(n+r_0)
  
  q_kn <- 1/(n+d_0-2*k-1) + (2*n+r_0+d_0-2*k-1)/((n+r_0)*(n+d_0-2*k-1)*(n+d_0-2*k-2))
  
  Q_c <- solve(S_mat_c) - (solve(S_mat_c)%*%ones%*%t(ones)%*%solve(S_mat_c))/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 
  
  # Weights for global minimum variance portfolio under conjugate prior
  wts_gmv_c <- (solve(S_mat_c)%*%ones)/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 
  
  ###########
  ## Return expected values
  # Sample?
  R_sample <- t(mean_ret)%*%wts_gmv
  # Frequentist
  R_mv_freq <- t(mean_ret)%*%wts_gmv + 1/(alpha*d_n)*t(mean_ret)%*%Q%*%mean_ret
  # Diffuse
  R_mv_diffuse <- t(mean_ret)%*%wts_gmv + 1/(alpha*c_kn)*t(mean_ret)%*%Q%*%mean_ret
  # Conjugate
  R_mv_conjugate <- t(mean_ret_c)%*%wts_gmv_c + 1/(alpha*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c
  
  ## Expected variance
  # Sample?
  V_sample <- d_n/(as.double(t(ones)%*%solve(S_mat)%*%ones))
  # Frequentist
  V_mv_freq <- d_n/(as.double(t(ones)%*%solve(S_mat)%*%ones)) + 1/(alpha^2*d_n)*t(mean_ret)%*%Q%*%mean_ret
  # Diffuse
  V_mv_diffuse <- c_kn/(as.double(t(ones)%*%solve(S_mat)%*%ones)) + 1/(alpha^2*c_kn)*t(mean_ret)%*%Q%*%mean_ret
  # Conjugate
  V_mv_conjugate <- q_kn/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) + 1/(alpha^2*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c
  
  return(cbind(R_sample, R_mv_freq, R_mv_diffuse, R_mv_conjugate, V_sample, V_mv_freq, V_mv_diffuse, V_mv_conjugate))
}

cdratio <- function(k=5,n=50){
  c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
  d_n <- 1/(n-1)
  return(c_kn/d_n)
}


data <- meanVarianceWeights(B=40, k=5, n = 40)

data
abs(data[,1]-data[,2])
abs(data[,1]-data[,3])
abs(data[,1]-data[,4])
abs(data[,5]-data[,6])
abs(data[,5]-data[,7])
abs(data[,5]-data[,8])

data2 <- meanVarianceWeights(B=40, k=40, n = 40)

data2
abs(data2[,1]-data2[,2])
abs(data2[,1]-data2[,3])
abs(data2[,1]-data2[,4])
abs(data2[,5]-data2[,6])
abs(data2[,5]-data2[,7])
abs(data2[,5]-data2[,8])






# Used to look at certain quantities
RMV <- function(B = 10000, k = 5, n = 50, alpha = 50, meanlow = -0.01, meanhigh = 0.01, 
                volatilitylow = 0.002, volatilityhigh = 0.005, corr = 0.6,
                r_0 = 100, d_0 = 100){
  ones <- t(t(rep(1, k)))
  samples <- c()
  for(i in 1:B){
    # Parameters
    mu <- runif(k, meanlow, meanhigh)
    volatility <- runif(k, volatilitylow, volatilityhigh)
    correlation <- (1-corr)*diag(k)+corr*ones(k,k)
    S_mat <- diag(volatility)%*%correlation%*%diag(volatility)
    # Sample from multivariate normal distribution
    samples <- rbind(samples,mvrnorm(n=1, mu = mu, Sigma = S_mat))
  }
  
  mean_ret <- colMeans(samples)
  
  S_mat <- t(samples - mean_ret)%*%(samples - mean_ret)
  ### Frequentist
  d_n <- 1/(n-1)
  
  # Weights for global minimum variance portfolio
  wts_gmv <- (solve(S_mat)%*%ones)/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 
  
  # Q-matrix used to solve the optimization problem of the quadratic utility funciton
  Q <- solve(S_mat) - (solve(S_mat)%*%ones%*%t(ones)%*%solve(S_mat))/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 
  
  ### Deriving quantities for diffuse prior
  c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
  ## Return expected values
  # Frequentist
  R_mv_freq <- t(mean_ret)%*%wts_gmv + 1/(alpha*d_n)*t(mean_ret)%*%Q%*%mean_ret
  # Diffuse
  R_mv_diffuse <- t(mean_ret)%*%wts_gmv + 1/(alpha*c_kn)*t(mean_ret)%*%Q%*%mean_ret
  #return(R_mv_freq, R_mv_diffuse)
  return(Q)
}  



Q<-RMV()
















