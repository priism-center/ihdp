############################################
# Test script

############################################
args <- commandArgs(trainlingOnly=TRUE)

fileConn <- file('test.txt')
writeLines(c("Hello World", args), fileconn)
close(fileConn)
