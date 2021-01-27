library(tidyquant)
library(plotly)
library(tidyverse)
library(timetk)
# get list of names from a imported yahoo portfolio
source("getNamesFinance.R")
# Gathers stock data in nested list
#source("getData.R")


stockNames <- getNames("Data/quotes.csv")

tick <- stockNames[-1] %>% head(6) %>% str_sort()


#download price data
price_data <- tq_get(tick,
                     from = '2015-12-29',
                     to = '2019-12-29',
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

# Create the necessary variables from "Bayesian mean variance optimal portfolio selection under parameter..."
#################
n <- nrow(ret_data_wide)
k <- ncol(ret_data_wide)
c_kn <- 1/(n-k-1)+(2*n-k-1)/(n*(n-k-1)*(n-k-2))

mean_ret <- colMeans(ret_data_wide)

# multiply to annualize maybe?
cov_mat <- cov(ret_data_wide)

# Vector of ones
ones <- t(t(rep(1, length(tick))))

# Weights
wts <- (solve(cov_mat)%*%ones)/(as.double(t(ones)%*%solve(cov_mat)%*%ones)) 

# Q from Eq. 8 and R from Eq. 28
Q <- solve(cov_mat) - (solve(cov_mat)%*%ones%*%t(ones)%*%solve(cov_mat))/(as.double(t(ones)%*%solve(cov_mat)%*%ones)) 












