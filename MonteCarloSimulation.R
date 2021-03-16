library(matlab)
library(dplyr)
library(LaplacesDemon)

## Sample from marginal posterior of Sigma for the Objective-based prior
GenerateDataEnvelope <- function(n, acceptMax = 1000000, k = 5, mean_ret = mean_ret, S = S_mat, gamma = 50){
  set.seed(13)
  # Parameters
  S_inv <- solve(S) # invert once to minimize number of computations
  w_ob <- rep(1/k,k) # suitable prior constant
  sigma <- 0.01 # uncertainty of mu
  s2 <- mean(diag(S_mat)) # average of diag of Sigma
  S_ob <- S_mat # prior constant
  v_ob <- k# prior constant
  # Marginalposterior
  f <- function(s2, sigma, gamma, S, S_inv, w_ob, n, mean_ret, v_ob, S_ob){
    r_ob <- (s2/sigma*gamma*S%*%w_ob + n*mean_ret)/(s2/sigma+n)
    
    det(S)^(-(v_ob+n)/2)*exp(-1/2*tr((S_ob + (n-1)*S)*S_inv))*
      exp(-1/2(n*t(mean_ret)%*%S_inv%*%mean_ret+s2/sigma*t(w_ob)%*%S%*%w_ob-(n+s2/sigma)*t(r_ob)%*%S_inv%*%r_ob))
  }
  # Rejection sampling
  accept <- c()
  yesCount <- 0
  sample.x <- c()
  for(i in 1:n){
    U = runif(1, 0, 1)
    sample = rinvwishart(v_ob, S_ob + (n-1)*S)
    sample.x[i] <- sample
    if(U <= f(s2, sigma, gamma, S, S_inv, w_ob, n, mean_ret, v_ob, S_ob)/(2*dinvwishart(sample, v_ob, S_ob + (n-1)*S))) {
      accept[i] = 'Yes'
      
      yesCount <- yesCount + 1
      if(yesCount == acceptMax){
        break
      }
    }
    else {
      accept[i] = 'No'
    }
  }
  T = data.frame(sample.x, accept = factor(accept, levels= c('Yes','No')))
  return(T)
}




















