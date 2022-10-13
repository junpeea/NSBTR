

library(rstan)
getwd()

rt <- stanc(file="210601_test_stan.stan", model_name="8schools")

schools.data <- list(
  n = 8,
  y = c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
)

fit1 <- stan(
  file = "210601_test_stan.stan",  # Stan program
  data = schools.data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  refresh = 1000          # show progress every 'refresh' iterations
)

remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")