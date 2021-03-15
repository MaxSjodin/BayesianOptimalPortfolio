# This script includes functions that....
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

tick <- stockNames[-1] %>% head(8) %>% str_sort()

fullNames <- get_stock_name(tick)

#download price data
price_data <- tq_get(tick,
                     from = '2021-01-01',
                     to = '2021-03-01',
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

#Test one portfolio to look for errors
#################
mean_ret <- colMeans(ret_data_wide)
mean_ret
# multiply to annualize
cov_mat <- cov(ret_data_wide)*252

wts <- runif(n = length(tick))
wts <- wts/sum(wts)

port_returns <- (sum(wts * mean_ret) + 1)^252 - 1
print(port_returns)

port_risk <- sqrt(t(wts) %*% (cov_mat %*% wts))
print(port_risk)

## Theoretical minimum variance portfolio weights

# Column vector of ones
ones <- t(t(rep(1, length(tick))))

denominator <- t(ones)%*%solve(cov_mat)%*%ones

weights <- (solve(cov_mat)%*%ones)/(as.double(denominator)) 

((sum(mean_ret*weights)+1)^252)-1

sqrt(t(weights) %*% (cov_mat %*% weights))

## Below follows code to optimize portfolio and store values

# Number of portfolios
num_port <- 5000

# Creating a matrix to store the weights

all_wts <- matrix(nrow = num_port,
                  ncol = length(tick))

# Creating an empty vector to store
# Portfolio returns

port_returns <- vector('numeric', length = num_port)

# Creating an empty vector to store
# Portfolio Standard deviation

port_risk <- vector('numeric', length = num_port)

# Creating an empty vector to store
# Portfolio Sharpe Ratio

sharpe_ratio <- vector('numeric', length = num_port)

# For loop over random portfolios to 
for (i in seq_along(port_returns)) {
  
  wts <- runif(length(tick))
  wts <- wts/sum(wts)
  wts
  
  # Storing weight in the matrix
  all_wts[i,] <- wts
  
  # Portfolio returns
  
  port_ret <- sum(wts * mean_ret)
  port_ret <- ((port_ret + 1)^252) - 1
  
  # Storing Portfolio Returns values
  port_returns[i] <- port_ret
  
  
  # Creating and storing portfolio risk
  port_sd <- sqrt(t(wts) %*% (cov_mat  %*% wts))
  port_risk[i] <- port_sd
  
  # Creating and storing Portfolio Sharpe Ratios
  # Assuming 0% Risk free rate
  
  sr <- port_ret/port_sd
  sharpe_ratio[i] <- sr
  
}


# Storing the values in the table
portfolio_values <- tibble(Return = port_returns,
                           Risk = port_risk,
                           SharpeRatio = sharpe_ratio)


# Converting matrix to a tibble and changing column names
all_wts <- tk_tbl(all_wts)

colnames(all_wts) <- colnames(ret_data_wide)

# Combing all the values together
portfolio_values <- tk_tbl(cbind(all_wts, portfolio_values))

min_var <- portfolio_values[which.min(portfolio_values$Risk),]
max_sr <- portfolio_values[which.max(portfolio_values$SharpeRatio),]

# Plot Minimum Variance
p <- min_var %>%
  gather(as.factor(tick[1]):as.factor(tick[length(tick)]), key = Asset,
         value = Weights) %>%
  mutate(Asset = as.factor(Asset)) %>%
  ggplot(aes(x = fct_reorder(Asset,Weights), y = Weights, fill = Asset)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(x = 'Assets', y = 'Weights', title = "Minimum Variance Portfolio Weights") +
  scale_y_continuous(labels = scales::percent) 

ggplotly(p)

# Plot Tangent
p <- max_sr %>%
  gather(as.factor(tick[1]):as.factor(tick[length(tick)]), key = Asset,
         value = Weights) %>%
  mutate(Asset = as.factor(Asset)) %>%
  ggplot(aes(x = fct_reorder(Asset,Weights), y = Weights, fill = Asset)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  labs(x = 'Assets', y = 'Weights', title = "Tangency Portfolio Weights") +
  scale_y_continuous(labels = scales::percent) 

ggplotly(p)

# All random portfolios to visualize efficient frontier

p <- portfolio_values %>%
  ggplot(aes(x = Risk, y = Return, color = SharpeRatio)) +
  geom_point() +
  theme_classic() +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = 'Annualized Risk',
       y = 'Annualized Returns',
       title = "Portfolio Optimization & Efficient Frontier") +
  geom_point(aes(x = Risk,
                 y = Return), data = min_var, color = 'red') +
  geom_point(aes(x = Risk,
                  y = Return), data = max_sr, color = 'red') #+
  # annotate('text', x = 0.20, y = 0.42, label = "Tangency Portfolio") +
  # annotate('text', x = 0.18, y = 0.01, label = "Minimum variance portfolio") +
  # annotate(geom = 'segment', x = 0.14, xend = 0.135,  y = 0.01, 
  #          yend = 0.06, color = 'red', arrow = arrow(type = "open")) +
  # annotate(geom = 'segment', x = 0.22, xend = 0.2275,  y = 0.405, 
  #          yend = 0.365, color = 'red', arrow = arrow(type = "open"))


ggplotly(p)
