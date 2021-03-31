library(matlab)
library(dplyr)
library(MASS)
library(MCMCpack)

meanVarianceWeights <- function(B = 1000, k = 40, n = 50, alpha = 50, meanlow = -0.01, meanhigh = 0.01, 
                       volatilitylow = 0.002, volatilityhigh = 0.005, corr = 0.6,
                       r_0 = 100, d_0 = 100){

 
  ones <- t(t(rep(1, k))) # Vector of ones
  samples <- c()
  samples_c <- c()
  df <- data.frame(matrix(0, ncol = 8, nrow = B))
  colnames(df) <- c("Rpop", "Vpop", "Rfreq", "Vfreq", "Rdif", "Vdif", "Rconj", "Vconj")
  for(i in 1:B){
    # Parameters
    mu_true <- runif(k, meanlow, meanhigh)
    volatility <- runif(k, volatilitylow, volatilityhigh)
    correlation <- (1-corr)*diag(k)+corr*ones(k,k)
    Sigma_true <- diag(volatility)%*%correlation%*%diag(volatility)
    
    m_0 <- mu_true + 0.5*runif(k, meanlow, meanhigh)
    S_0 <- Sigma_true + 0.5*diag(runif(k, volatilitylow, volatilityhigh))
    
    Sigma_true_c <- riwish(d_0,S_0)
    mu_true_c <- mvrnorm(n = 1, mu = m_0, Sigma = Sigma_true_c/r_0)
    
    # Sample from multivariate normal distribution
    # Loop over all, take average of all R_mv
    # set n=140 and use mvrnorm once per R_mv calculated
    samples <- mvrnorm(n = n, mu = mu_true, Sigma = Sigma_true)
    # draw another sample from conjugate prior
    samples_c <- mvrnorm(n = n, mu = mu_true_c, Sigma = Sigma_true_c)
    
    mean_ret <- colMeans(samples)
    
    mean_ret_c <- colMeans(samples_c)
    
    m_mean <- matrix(data=1, nrow=n) %*% mean_ret
    
    m_mean_c <- matrix(data=1, nrow=n) %*% mean_ret_c
    
    m_deviations <- samples - m_mean
    
    m_deviations_c <- samples_c - m_mean_c
    
    S_mat <- t(m_deviations)%*%m_deviations
    
    S_mat_2 <- t(m_deviations_c)%*%m_deviations_c
    #S_mat <- t(samples - mean_ret)%*%(samples - mean_ret) #Calculate covariance matrix properly
    
    #S_mat_2 <- t(samples_c - mean_ret_c)%*%(samples_c - mean_ret_c) # S_mat_c already taken (name)
    
    # Is mu_true and Sigma_true generated for each sample???
    
    #return(samples)
    ### Frequentist
    d_n <- 1/(n-1)
    
    # Weights for global minimum variance portfolio
    wts_gmv <- (solve(S_mat)%*%ones)/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 
    
    # based on population variance matrix
    wts_gmv_true <- (solve(Sigma_true)%*%ones)/(as.double(t(ones)%*%solve(Sigma_true)%*%ones)) 
    
    # Q-matrix used to solve the optimization problem of the quadratic utility function
    Q <- solve(S_mat) - (solve(S_mat)%*%ones%*%t(ones)%*%solve(S_mat))/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 
    
    # based on population variance matrix
    Q_true <- solve(Sigma_true) - (solve(Sigma_true)%*%ones%*%t(ones)%*%solve(Sigma_true))/(as.double(t(ones)%*%solve(Sigma_true)%*%ones)) 
    
    ### Deriving quantities for diffuse prior
    c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
    
    ### Deriving quantities for conjugate prior
    mean_ret_c <- (n*mean_ret_c+r_0*m_0)/(n+r_0)

    S_mat_c <- S_mat_2+S_0+n*r_0*(m_0-mean_ret_c)%*%t(m_0-mean_ret_c)/(n+r_0)
    
    ##############Based on same sample as freq
    mean_ret_c <- (n*mean_ret+r_0*m_0)/(n+r_0)
    S_mat_c <- S_mat+S_0+n*r_0*(m_0-mean_ret_c)%*%t(m_0-mean_ret_c)/(n+r_0)
    ################
    q_kn <- 1/(n+d_0-2*k-1) + (2*n+r_0+d_0-2*k-1)/((n+r_0)*(n+d_0-2*k-1)*(n+d_0-2*k-2))
    
    Q_c <- solve(S_mat_c) - (solve(S_mat_c)%*%ones%*%t(ones)%*%solve(S_mat_c))/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 
    
    # Weights for global minimum variance portfolio under conjugate prior
    wts_gmv_c <- (solve(S_mat_c)%*%ones)/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 
    
    #mean_ret <- t(t(mean_ret))
    #mean_ret_c <- t(t(mean_ret_c))
    ###########
    ## Return expected values
    # Population
    R_mv_pop <- t(mu_true)%*%wts_gmv_true + 1/(alpha)*t(mu_true)%*%Q_true%*%mu_true
    # Frequentist
    R_mv_freq <- t(mean_ret)%*%wts_gmv + 1/(alpha*d_n)*t(mean_ret)%*%Q%*%mean_ret
    # Diffuse
    R_mv_diffuse <- t(mean_ret)%*%wts_gmv + 1/(alpha*c_kn)*t(mean_ret)%*%Q%*%mean_ret
    # Conjugate
    R_mv_conjugate <- t(mean_ret_c)%*%wts_gmv_c + 1/(alpha*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c

    ## Expected variance
    # Sample?
    V_mv_pop <- 1/(as.double(t(ones)%*%solve(Sigma_true)%*%ones)) + 1/(alpha^2)*t(mu_true)%*%Q_true%*%mu_true
    # Frequentist
    V_mv_freq <- d_n/(as.double(t(ones)%*%solve(S_mat)%*%ones)) + 1/(alpha^2*d_n)*t(mean_ret)%*%Q%*%mean_ret
    # Diffuse
    V_mv_diffuse <- c_kn/(as.double(t(ones)%*%solve(S_mat)%*%ones)) + 1/(alpha^2*c_kn)*t(mean_ret)%*%Q%*%mean_ret
    # Conjugate
    V_mv_conjugate <- q_kn/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) + 1/(alpha^2*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c
    # Add to df
    df[i,] <- c(R_mv_pop, V_mv_pop, R_mv_freq, V_mv_freq, R_mv_diffuse, V_mv_diffuse, R_mv_conjugate, V_mv_conjugate)
    #return(cbind(t(mu_true)%*%wts_gmv_true, t(mean_ret)%*%wts_gmv))
    ####return(cbind(Q/c_kn, Q_true))
  }    
  #return(cbind(R_mv_pop, R_mv_freq, R_mv_diffuse, R_mv_conjugate, V_mv_pop, V_mv_freq, V_mv_diffuse, V_mv_conjugate))
  
  #return(cbind(t(mu_true)%*%wts_gmv_true, t(mean_ret)%*%wts_gmv))
  #return(cbind(1/(alpha)*t(mu_true)%*%Q_true%*%mu_true, 1/(alpha*d_n)*t(mean_ret)%*%Q%*%mean_ret))
  
  return(df)
}

# Look at mean for each expected return/variance
data <- meanVarianceWeights() %>% colMeans()
data

# function that fills matrix with values, similar to that of 
# table 1 in "Bayesian mean-variance... under parameter uncertainty"
deviations <- function(k_vec = c(5, 10, 25, 40), n_vec = c(50, 75, 100, 125), df_length = 8){
  n_dev <-df_length/2-1 # number of priors to compute deviations for
  m_ret_dev <- matrix(0, nrow = length(k_vec)*n_dev, ncol = length(n_vec)) 
  m_var_dev <- matrix(0, nrow = length(k_vec)*n_dev, ncol = length(n_vec))
  for (j in 1:length(k_vec)) {
    for (i in 1:length(n_vec)) {
      temp <- meanVarianceWeights(k=k_vec[j],n=n_vec[i]) %>% colMeans()
      
      ret_deviations <- c(abs(temp[1]-temp[3]), abs(temp[1]-temp[5]), abs(temp[1]-temp[7]))
      var_deviations <- c(abs(temp[2]-temp[4]), abs(temp[2]-temp[6]), abs(temp[2]-temp[8]))
      
      m_ret_dev[(j*n_dev-2):(j*n_dev),i] <- ret_deviations
      m_var_dev[(j*n_dev-2):(j*n_dev),i] <- var_deviations
    }
  }  
  return(list(ret = m_ret_dev, var = m_var_dev))
}

data <- deviations(k_vec = c(5, 10, 25, 40), n_vec = c(50, 75, 100, 125))
data


cdratio <- function(k=5,n=50){
  c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
  d_n <- 1/(n-1)
  q_kn <- 1/(n+100-2*k-1) + (2*n+100+100-2*k-1)/((n+100)*(n+100-2*k-1)*(n+100-2*k-2))
  return(c(c_kn,q_kn))
}
cdratio(k=40,n=50)
















