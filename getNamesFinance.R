#Retrieving the names from a yahoo finance exported file in csv format
getNames <- function(filepath){
  data <- read.csv(file=filepath)
  return(data[,1])
}

