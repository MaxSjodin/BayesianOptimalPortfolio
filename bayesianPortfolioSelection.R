library(tidyquant)
library(plotly)
library(tidyverse)
library(timetk)
# get list of names from a imported yahoo portfolio
source("getNamesFinance.R")
# Gathers stock data in nested list
#source("getData.R")


stockNames <- getNames("Data/quotes.csv")

tick <- stockNames[-1] %>% head(29) %>% str_sort()

#download price data
price_data <- tq_get(tick,
                     from = '2020-05-28',
                     to = '2020-12-29',
                     get = 'stock.prices')

price_data[,-c(1,2)] <- na.approx(price_data[,-c(1,2)])

# Calculate daily returns for assets
ret_data <- price_data %>%
  group_by(symbol) %>%
  tq_transmute(select = adjusted,
               mutate_fun = periodReturn,
               period = "daily",
               col_rename = "ret")

# Wide format
ret_data_wide <- ret_data %>%
  spread(symbol, value = ret) %>%
  tk_xts()

# From the wide format we can perhaps select columns and rewrite the rest of the code as function

### Diffuse prior
# Create the necessary variables from "Bayesian mean variance optimal portfolio selection under parameter..."
#################
n <- nrow(ret_data_wide)
k <- ncol(ret_data_wide)
c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
# Arbitrary value for gamma
gamma <- 10

mean_ret <- colMeans(ret_data_wide)

# multiply to annualize maybe?
cov_mat <- cov(ret_data_wide)*(nrow(ret_data_wide)-1)

# Vector of ones
ones <- t(t(rep(1, length(tick))))

# Weights for global minimum variance portfolio
wts_gmv <- (solve(cov_mat)%*%ones)/(as.double(t(ones)%*%solve(cov_mat)%*%ones)) 

# Q from Eq. 8 and R from Eq. 28
Q <- solve(cov_mat) - (solve(cov_mat)%*%ones%*%t(ones)%*%solve(cov_mat))/(as.double(t(ones)%*%solve(cov_mat)%*%ones)) 




## mean variance portfolio
wts_mv_diffuse <- wts_gmv + 1/(gamma*c_kn)*Q%*%mean_ret

# Same as Eq. 9
R_mv_diffuse <- t(mean_ret)%*%wts_gmv + 1/(gamma*c_kn)*t(mean_ret)%*%Q%*%mean_ret

# Expected annual return
(R_mv_diffuse + 1)^252 - 1

# Expected variance
V_mv_diffuse <- c_kn/(as.double(t(ones)%*%solve(cov_mat)%*%ones)) + 1/(gamma^2*c_kn)*t(mean_ret)%*%Q%*%mean_ret

# Expected annual variance
V_mv_diffuse*252

# Expected standard deviance?
sqrt(V_mv_diffuse)


## Bayesian Efficient frontier

# Return of global minimum variance portfolio
R_gmv <- t(mean_ret)%*%wts_gmv

# Expected annual return
(R_gmv + 1)^252 - 1

# Variance of global minimum variance portfolio
V_gmv <- c_kn/(as.double(t(ones)%*%solve(cov_mat)%*%ones)) 

# Left to do is create efficient frontier using Eq. 13



#Not runnable until good hyperparameter-values of 
### Extended Black-Litterman (Conjugate Prior)

m_0 <- 0
r_0 <- 0
d_0 <- 0
S_0 <- 0


q_kn <- 1/(n+d_0-2*k-1) + (2*n+r_0+d_0-2*k-1)/((n+r_0)*(n+d_0-2*k-1)*(n+d_0-2*k-2))

# Mean variance portfolio

wts_mv_conjugate <- wts_gmv + 1/(gamma*q_kn)*Q%*%mean_ret

# Same as Eq. 9
R_mv_conjugate <- t(mean_ret)%*%wts_gmv + 1/(gamma*q_kn)*t(mean_ret)%*%Q%*%mean_ret

# Expected annual return
(R_mv_conjugate + 1)^252 - 1

# Expected variance
V_mv_conjugate <- q_kn/(as.double(t(ones)%*%solve(cov_mat)%*%ones)) + 1/(gamma^2*q_kn)*t(mean_ret)%*%Q%*%mean_ret

# Expected annual variance
V_mv_conjugate*252

# Expected standard deviance?
sqrt(V_mv_conjugate)





