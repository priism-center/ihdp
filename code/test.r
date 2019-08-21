############################################
# Test script

############################################
args <- commandArgs(trailingOnly=TRUE)

fileConn <- file('test.txt')
writeLines(c("Hello World", args), fileconn)
close(fileConn)
