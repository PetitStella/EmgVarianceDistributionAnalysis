
source("EM_algorithm.R")

library(MASS)

data <- scan("emg_data.dat")

result <- estimate.var.dist(data)

print(result)