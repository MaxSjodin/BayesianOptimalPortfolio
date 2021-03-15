library(tidyquant)
library(plotly)
library(tidyverse)
library(timetk)
# get list of names from a imported yahoo portfolio
source("getNamesFinance.R")
# Gathers stock data in nested list
#source("getData.R")



# Stocks from arbitrary "Sweden top 100" portfolio
stockNames <- getNames("Data/quotes.csv")
# All stocks on Swedish market
allstockNames <- getNames("Data/stocks.csv")

# Create data frame containing adjusted price for each stock and time
# Maybe later

# Doing analysis from post "https://www.codingfinance.com/post/2018-04-20-portfolio-stats/" for 
#Select stocks to analyze

tick <- stockNames[-1] %>% head(21) %>% str_sort()

fullNames <- get_stock_name(tick)

#download price data
price_data <- tq_get(tick,
                     from = '2021-01-01',
                     to = '2021-03-02',
                     get = 'stock.prices')


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

# Approximate missing values 
# Should not affect daily data but in case of weekly/monthly data the NA approx 
# should be done before gathering to get weekly/monthly
ret_data_wide <- na.approx(ret_data_wide)

# From the wide format we can perhaps select columns and rewrite the rest of the code as function

## Global variables used in multiple approaches/priors
# function for expected annual return
expectedAnnualReturn <- function(R){
  (R + 1)^252-1
}
n <- nrow(ret_data_wide)
k <- ncol(ret_data_wide)
mean_ret <- colMeans(ret_data_wide)
# Multiply covariance by n-1 to get the S matrix in Eq.4
S_mat <- cov(ret_data_wide)*(nrow(ret_data_wide)-1)

# Vector of ones
ones <- t(t(rep(1, length(tick))))

# Weights for global minimum variance portfolio
wts_gmv <- (solve(S_mat)%*%ones)/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 

### Frequentist
d_n <- 1/(n-1)

# Return of global minimum variance portfolio
R_gmv_freq <- t(mean_ret)%*%wts_gmv


# Variance of global minimum variance portfolio
V_gmv_freq <- d_n/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 



### Diffuse prior
# Create the necessary variables from "Bayesian mean variance optimal portfolio selection under parameter..."
#################
c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))
# Arbitrary value for gamma
gamma <- 5

## Algorithm 1: Simulate from posterior distribution of diffuse prior from minimum variance portfolio
# mean and weights are column vectors
sampleDiffuse <- function(Samples = 10000, mean, wts, S, n, k){
  B <- Samples
  sample_diffuse <- c()
  
  for(i in 1:B){
    # Draw samples from t distribution
    t_1 <- rt(1, n-k)
    t_2 <- rt(1, n-k+1)
    
    sample_diffuse[i] <- t(wts)%*%mean + sqrt(t(wts)%*%S%*%wts)*((t_1)/sqrt(n*(n-k))+sqrt(1+(t_1^2)/(n-k))*((t_2)/sqrt(n-k+1)))
  }
  return(sample_diffuse)
}

sample_diffuse <- sampleDiffuse(Samples = 10000, mean = mean_ret, S = S_mat, wts = wts_gmv, n = n, k = k)

# Does not seem to be correct, compared to Table1 and 2 in article
mean(sample_diffuse)
var(sample_diffuse)

## mean variance portfolio
# Q from Eq. 8 and R from Eq. 28
Q <- solve(S_mat) - (solve(S_mat)%*%ones%*%t(ones)%*%solve(S_mat))/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 

wts_mv_diffuse <- wts_gmv + 1/(gamma*c_kn)*Q%*%mean_ret

# Same as Eq. 9
R_mv_diffuse <- t(mean_ret)%*%wts_gmv + 1/(gamma*c_kn)*t(mean_ret)%*%Q%*%mean_ret

# Expected annual return
#(R_mv_diffuse + 1)^252 - 1

# Expected variance
V_mv_diffuse <- c_kn/(as.double(t(ones)%*%solve(S_mat)%*%ones)) + 1/(gamma^2*c_kn)*t(mean_ret)%*%Q%*%mean_ret

# Expected annual variance
#V_mv_diffuse*252

# Expected standard deviance?
#sqrt(V_mv_diffuse)


## Bayesian Efficient frontier

# Return of global minimum variance portfolio
R_gmv_diffuse <- t(mean_ret)%*%wts_gmv

# Expected annual return
#(R_gmv + 1)^252 - 1

# Variance of global minimum variance portfolio
V_gmv_diffuse <- c_kn/(as.double(t(ones)%*%solve(S_mat)%*%ones)) 

# Left to do is create efficient frontier using Eq. 13
# Given a variance V, return the expected return R. Does not work
efficientFrontierDiffuse <- function(V){
  (R-R_gmv_diffuse)^2 = 1/(c_kn)*t(mean_ret)%*%Q%*%mean_ret*(V-V_gmv_diffuse)
  return(R)
}
# Given a expected return R, return the variance V.
efficientFrontierDiffuse2 <- function(R){
  V <- c_kn/t(mean_ret)%*%Q%*%mean_ret*(R-R_gmv_diffuse)^2 + V_gmv_diffuse
  return(cbind(R,V))
}

base <- ggplot() + xlim(0.0006, 0.01)
base + geom_function(fun = function(x) c_kn/t(mean_ret)%*%Q%*%mean_ret*(x-R_gmv_diffuse)^2 + V_gmv_diffuse) +
        coord_flip()

efficientFrontierDiffuse2(2.5)


#Not runnable until good hyperparameter-values of 
### Extended Black-Litterman (Conjugate Prior)

m_0 <- mean_ret+0.001
r_0 <- 100
d_0 <- 100
S_0 <- S_mat+0.001

mean_ret_c <- (n*mean_ret+r_0*m_0)/(n+r_0)

S_mat_c <- S_mat+S_0+n*r_0*(m_0-mean_ret_c)%*%t(m_0-mean_ret_c)/(n+r_0)

q_kn <- 1/(n+d_0-2*k-1) + (2*n+r_0+d_0-2*k-1)/((n+r_0)*(n+d_0-2*k-1)*(n+d_0-2*k-2))

Q_c <- solve(S_mat_c) - (solve(S_mat_c)%*%ones%*%t(ones)%*%solve(S_mat_c))/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 

wts_gmv_c <- (solve(S_mat_c)%*%ones)/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 

# Mean variance portfolio

wts_mv_conjugate <- wts_gmv_c + 1/(gamma*q_kn)*Q_c%*%mean_ret_c

# Same as Eq. 9
R_mv_conjugate <- t(mean_ret_c)%*%wts_gmv_c + 1/(gamma*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c

# Expected annual return
(R_mv_conjugate + 1)^252 - 1

# Expected variance
V_mv_conjugate <- q_kn/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) + 1/(gamma^2*q_kn)*t(mean_ret_c)%*%Q_c%*%mean_ret_c

# Expected annual variance
#V_mv_conjugate*252

# Expected standard deviance?
#sqrt(V_mv_conjugate)

# Return of global minimum variance portfolio
R_gmv_conjugate <- t(mean_ret_c)%*%wts_gmv_c

# Expected annual return
#(R_gmv + 1)^252 - 1

# Variance of global minimum variance portfolio
V_gmv_conjugate <- q_kn/(as.double(t(ones)%*%solve(S_mat_c)%*%ones)) 

base <- ggplot() + xlim(0.0000005, 1)
base + geom_function(fun = function(x) c_kn/t(mean_ret)%*%Q%*%mean_ret*(x-R_gmv_diffuse)^2 + V_gmv_diffuse, aes(colour = "Diffuse")) +
  geom_function(fun = function(x) d_n/t(mean_ret)%*%Q%*%mean_ret*(x-R_gmv_freq)^2 + V_gmv_freq, aes(colour = "Freq")) +
  geom_function(fun = function(x) q_kn/t(mean_ret_c)%*%Q_c%*%mean_ret_c*(x-R_gmv_conjugate)^2 + V_gmv_conjugate, aes(colour = "Conjugate")) +
  #scale_x_continuous(labels = expectedAnnualReturn, limits = c(0, 0.08)) +
  ylim(0,0.8) +
  scale_color_manual(values = c('red','blue','green')) +
  labs(x="Return", y="Variance") +
  coord_flip()



