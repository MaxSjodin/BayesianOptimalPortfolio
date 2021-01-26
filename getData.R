library(quantmod)
library(tidyquant)

#Takes imported yahoo finance portfolio and retrieves the stock names i.e. ABB.ST
source("getNamesFinance.R")

#Creates list, each element contains object of type xts for given stock and timeperiod
getStocks <- function(names = stockNames, startDate = '2011-01-01', endDate = '2021-01-01'){
  stockData <- list()
  for(i in 1:length(stockNames)){
    getSymbols(names[i], from = startDate,
                 to = endDate, src = "yahoo", warnings = FALSE,
                 auto.assign = TRUE)
    stockData[[i]] <- tq_get(names[i], from = startDate,
                             to = endDate, get="stock.prices")
  }
  return(stockData)
}



