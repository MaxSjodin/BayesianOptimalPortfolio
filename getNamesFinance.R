#Retrieving the names from a yahoo finance exported file in csv format
library(httr)
library(jsonlite)


getNames <- function(filepath){
  data <- read.csv(file=filepath)
  return(data[,1])
}

get_stock_name <- function(tick) {
  stock_symbol <- c()
  stock_name <- c()
  
  for (i in 1:length(tick)) {
  
    res = GET(paste("http://d.yimg.com/autoc.finance.yahoo.com/autoc?query=", tick[i], "&region=1&lang=en", sep=""))
    
    data = fromJSON(rawToChar(res$content))
    
    
    stock_symbol[i] = data$ResultSet$Result$symbol[1]
    stock_name[i] = data$ResultSet$Result$name[1]
  }  
  
  stocks <- cbind(symbol = stock_symbol, name = stock_name)
  
  return(stocks)
  
}

