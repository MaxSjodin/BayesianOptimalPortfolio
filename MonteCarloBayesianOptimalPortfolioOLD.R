library(matlab)
library(dplyr)

meanVarianceWeights <- function(Samples = 10000, k = 10, n = 50, alpha = 50, meanlow = -0.01, meanhigh = 0.01, 
                       volatilitylow = 0.002, volatilityhigh = 0.005, corr = 0.6,
                       r_0 = 100, d_0 = 100){

  ones <- t(t(rep(1, 10)))
  mean_ret <- runif(k, meanlow, meanhigh)
  volatility <- runif(k, volatilitylow, volatilityhigh)
  correlation <- (1-corr)*diag(k)+corr*ones(k,k)
  S_mat <- diag(volatility)%*%correlation%*%diag(volatility)
  
  m_0 <- mean_ret + 0.5*runif(k, meanlow, meanhigh)
  S_0 <- S_mat + 0.5*runif(k, volatilitylow, volatilityhigh)
  
  ### Frequentist
  d_n <- 1/(n-1)
  
  # Weights for global minimum variance portfolio
  wts_gmv <- (solve(S_mat)%*%ones)/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 
  
  # Q-matrix used to solve the optimization problem of the quadratic utility funciton
  Q <- solve(S_mat) - (solve(S_mat)%*%ones%*%t(ones)%*%solve(S_mat))/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 
  
  # Mean-Variance weights for diffuse prior
  wts_mv_freq <- wts_gmv + 1/(alpha*d_n)*Q%*%mean_ret  
  
  ### Deriving quantities for diffuse prior
  c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
  
  # Mean-Variance weights for diffuse prior
  wts_mv_diffuse <- wts_gmv + 1/(alpha*c_kn)*Q%*%mean_ret
  
  ### Deriving quantities for conjugate prior
  mean_ret_c <- (n*mean_ret+r_0*m_0)/(n+r_0)
  
  S_mat_c <- S_mat+S_0+n*r_0*(m_0-mean_ret_c)%*%t(m_0-mean_ret_c)/(n+r_0)
  
  q_kn <- 1/(n+d_0-2*k-1) + (2*n+r_0+d_0-2*k-1)/((n+r_0)*(n+d_0-2*k-1)*(n+d_0-2*k-2))
  
  Q_c <- solve(S_mat_c) - (solve(S_mat_c)%*%ones%*%t(ones)%*%solve(S_mat_c))/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 
  
  # Weights for global minimum variance portfolio under conjugate prior
  wts_gmv_c <- (solve(S_mat_c)%*%ones)/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 
  
  # Mean-Variance weights for conjugate prior
  wts_mv_conjugate <- wts_gmv_c + 1/(alpha*q_kn)*Q_c%*%mean_ret_c
  
  ## Return weights of mean-variance portfolio, for each prior
  meanVarWeights <- cbind(frequentist = wts_mv_freq, diffuse = wts_mv_diffuse, conjugate = wts_mv_conjugate)
  #return(meanVarWeights)
  
  ### Sample from frequentist perspective
  B <- Samples
  
  sample_freq <- c()
  
  
  ### Sample from posterior distribution of diffuse prior
  sample_diffuse <- c()
  
  for(i in 1:B){
    # Draw samples from t distribution
    t_1 <- rt(1, n-k)
    t_2 <- rt(1, n-k+1)
    
    sample_diffuse[i] <- t(wts_mv_diffuse)%*%mean_ret + sqrt(t(wts_mv_diffuse)%*%S_mat%*%wts_mv_diffuse)*((t_1)/sqrt(n*(n-k))+sqrt(1+(t_1^2)/(n-k))*((t_2)/sqrt(n-k+1)))
  }
  
  ### Sample from the posterior distribution of conjugate prior
  sample_conjugate <- c()
  
  for(i in 1:B){
    # Draw samples from t distribution
    eta_1 <- rt(1, n+d_0-2*k)
    eta_2 <- rt(1, n+d_0-2*k+1)
    
    sample_conjugate[i] <- t(wts_mv_conjugate)%*%mean_ret_c + sqrt(t(wts_mv_conjugate)%*%S_mat_c%*%wts_mv_conjugate)*((eta_1)/sqrt((n+r_0)*(n+d_0-2*k))+sqrt(1+(eta_1^2)/(n+d_0-2*k))*((eta_2)/sqrt(n+d_0-2*k+1)))
  }
  
  ## Return samples from the posterior distribution
  #Samples <- cbind(diffuse = sample_diffuse, conjugate = sample_conjugate)
  #return(Samples)
  
  
  ###########
  ## Return expected values
  # Diffuse
  R_mv_diffuse <- t(mean_ret)%*%wts_gmv + 1/(alpha*c_kn)*t(mean_ret)%*%Q%*%mean_ret
  # Conjugate
  R_mv_conjugate <- t(mean_ret_c)%*%wts_gmv_c + 1/(alpha*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c
  
  ## Expected variance
  # Diffuse
  V_mv_diffuse <- c_kn/(as.double(t(ones)%*%solve(S_mat)%*%ones)) + 1/(alpha^2*c_kn)*t(mean_ret)%*%Q%*%mean_ret
  # Conjugate
  V_mv_conjugate <- q_kn/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) + 1/(alpha^2*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c
  
  expectedMeanVar <- cbind(diffuseMean = R_mv_diffuse, diffuseVar = V_mv_diffuse, conjugateMean = R_mv_conjugate, conjugateVar = V_mv_conjugate)
  
  #return(expectedMeanVar)
  
  ###### Components of R_mv_diffuse and R_mv_conjugate
  return(cbind(mean(sample_diffuse), var(sample_diffuse), R_mv_diffuse, V_mv_diffuse, mean(sample_conjugate), var(sample_conjugate),R_mv_conjugate, V_mv_conjugate))
  
  #Return means
  #return(cbind(mean = mean_ret, meanc = mean_ret_c))
}

meanVarWeights <- meanVarianceWeights(n=100)
meanVarWeights

meanVarWeights[1,1]-meanVarWeights[1,3]
meanVarWeights[1,2]-meanVarWeights[1,4]
meanVarWeights[1,5]-meanVarWeights[1,7]
meanVarWeights[1,6]-meanVarWeights[1,8]




















