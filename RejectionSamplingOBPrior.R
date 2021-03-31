library(matlab)
library(dplyr)
library(LaplacesDemon)

## Sample from marginal posterior of Sigma for the Objective-based prior
RejectionSample <- function(acceptMax = 1, attempts = 1000, n = 100, k = 5, mean_ret = bar_x, Sig = Sig, S = S, w = w_ob, s2 = S2, s_ob = s_ob, 
                            SOB = SOB, v_ob = v_ob, gamma = 50){
  set.seed(13)
  # Parameters
  S<-S/(n-1)
  Sig_inv <- solve(Sig) # invert once to minimize number of computations
  # w_ob <- rep(1/k,k) # suitable prior constant
  # s_ob <- 100 # uncertainty of mu
  # s2 <- mean(diag(S_mat)) # average of diag of Sigma
  # SOB <- S_mat # prior constant
  # v_ob <- k# prior constant
  # Marginalposterior
  f <- function(s2, s_ob, gamma, S, Sig_inv, w_ob, n, mean_ret, v_ob, SOB){
    r_ob <- ((s2/s_ob)*gam*Sig%*%w_ob+n*mean_ret)/(s2/s_ob+n)
    
    det(S)^(-(v_ob+n)/2)*exp(-1/2*tr((SOB + (n-1)*S)*Sig_inv))*
      exp(-1/2*(n*(s2/s_ob+n)/(s2/s_ob+n)*t(mean_ret-gam*Sig%*%w_ob)%*%Sig_inv%*%(mean_ret-gam*Sig%*%w_ob)))
  }
  # Rejection sampling
  accept <- c()
  yesCount <- 0
  sample.x <- c()
  for(i in 1:attempts){
    U = runif(1, 0, 1)
    sample = rinvwishart(v_ob+n, SOB + (n-1)*S)
    sample.x <- sample
    if(U <= f(s2, s_ob, gamma, S, Sig_inv, w_ob, n, mean_ret, v_ob, SOB)/(2*sample.x)) {
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
  T = data.frame(sample.x)
  return(T)
}




















