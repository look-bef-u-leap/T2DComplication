# Description: This file is used for categorical conversion of trait variables

library(dplyr)

traRfVif2ndEd <- "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traRfVif2ndEd.rds" %>%
  readRDS(.)

trait6 <- traRfVif2ndEd

traCla <- trait6 %>%
  lapply(., class) %>%
  as.data.frame(.) %>%
  t(.) %>%
  as.data.frame(.)
traCla$name <- rownames(traCla)

num <- which((traCla[, 1] == "integer") | (traCla[, 1] == "numeric"))
traNum <- trait6[, num]

# Check whether the column meets the quartile split condition
n <- as.numeric()
for (i in 1:ncol(traNum)) {
  anyD <- anyDuplicated(quantile(traNum[, i], probs = seq(0, 1, 0.25), na.rm = TRUE))
  if (anyD == 0){
    
  }else{
    print(i)
    n <- c(n, i)
  }
}
# The distribution of these variables is difficult to classify and is no longer included
# Remove variables based on n value
if (length(n) == 0){
  traNum2 <- traNum
}else{
  traNum2 <- traNum[, -n]
}

convertLab <- function(column) {
  quartiles <- quantile(column, probs = seq(0, 1, 0.25), na.rm = TRUE)
  categorical_column <- cut(column, breaks = quartiles, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
  return(categorical_column)
}
convertSec <- function(column) {
  quartiles <- quantile(column, probs = seq(0, 1, 0.25), na.rm = TRUE)
  categorical_column <- cut(column, breaks = quartiles, include.lowest = TRUE)
  return(categorical_column)
}
# Quickly get the format of the data frame
traNum2Sec <- traNum2
traNum2Lab <- traNum2

# The trait variables of numerical type are classified in two ways, section and label, 
# and the two ways are corresponding
n <- (ncol(traNum2) - 1)
for (i in 1:n) {
  traNum2Sec[, i+1] <- convertSec(traNum2[, i+1])
  traNum2Lab[, i+1] <- convertLab(traNum2[, i+1])
}
traLab2ndEd <- cbind(traNum2Lab, trait6[, -num])
traSec2ndEd <- cbind(traNum2Sec, trait6[, -num])

write.csv(
  traLab2ndEd, 
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traLab2ndEd.csv",
  row.names = TRUE
)
saveRDS(
  traLab2ndEd,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traLab2ndEd.rds"
)
data.table::fwrite(
  traSec2ndEd, 
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traSec2ndEd.csv",
  row.names = TRUE
)
saveRDS(
  traSec2ndEd,
  file = "D:\\work2\\t2dAndComplicationsUkb\\data\\processed\\integComp2ndEd\\traSec2ndEd.rds"
)
