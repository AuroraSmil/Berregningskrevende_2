## SINGLE SITE MCMC 

MH_single_alg <- function(n,data, t_0,t_2, t, lambda_0,lambda_1, beta, sigma){

  # Find number of deaths in each interval given initial t
  y_0 = sum(data[date <t]$event)
  y_1 = sum(data[date>= t]$event)
  
  # Construct result matrix
  results <- matrix(nrow = n, ncol = 5)
  results[1, 1:5] <- c(t, lambda_0, lambda_1, beta, 1)
  for(i in 2:n){
    print(i)
    # Collect old values for t and beta and compute corresponding y_0 and y_1 values
    t_old <- results[i-1, 1]
    
    y_0_o = sum(data[date <t_old]$event)
    y_1_o = sum(data[date>= t_old]$event)
    
    beta <- results[i-1, 4]
    
    # Preform Gibbs step
    lambda_0 <- rgamma(1, y_0 + 2, scale = 1/(t_old - t_0 + 1/beta))
    lambda_1 <- rgamma(1, y_1 + 2, scale = 1/(t_2 - t_old + 1/beta))

    beta_tilde <- rgamma(1, 6, scale = 1/(lambda_0 + lambda_1 + 1)) #is this okay to d0?
    beta <- 1/beta_tilde
    
    results[i, 2:4] <- c( lambda_0, lambda_1, beta)
    
    # Preform MH step with the random walk of order one as the proposal density
    t_new <- rnorm(1, t_old, sigma)
    
    ### Compute corresponding values for y_0 and y_1
    y_0_n = sum(data[date <t_new]$event)
    y_1_n = sum(data[date>= t_new]$event)
    
    ### Compute the acceptance probability in log scale and compare with u unif(0,1)
    a_log <- (((y_0_n +1)*log(lambda_0) + (y_1_n +1)*log(lambda_1) - t_new *(lambda_0 - lambda_1 )) - 
                ((y_0_o +1)*log(lambda_0) + (y_1_o +1)*log(lambda_1) - t_old *(lambda_0 - lambda_1 )))
    u <- runif(1)
    results[i, 5] <-  exp(a_log)
    if (u< exp(a_log)){
      results[i, 1] <-  t_new
    } else{
      results[i,1] <- t_old
      }
  }
  return (results)
}

# Load libraries and data
library(boot)
library(data.table)
library(tidyr)
library(ggplot2)

# Create data table
data <- as.data.table(coal)
data[, event:= 1]
data[, cum_event := cumsum(event)]
data[, date_shifted := date- data[1, date]]


q <- ggplot(data = data, aes(x = date, y = cum_event))
q <- q + geom_line()
q

q <- ggplot(data = data, aes(x = date_shifted, y = cum_event))
q <- q + geom_line()
q


# Initial values for t, lambda_0, lambda_1, beta
t = data[1, date] + 50
lambda_0 <- 10 
lambda_1 <- 5 
beta = 3

# Set n, t_0, t_2 and tuning parameter sigma
n = 10000
t_0 = data[1, date]
t_2 = data[nrow(data), date] #maybe should be 1963? 
sigma = 3

# Run algorithm
sim_MH_single <- MH_single_alg(n,data, t_0, t_2, t, lambda_0,lambda_1, beta, sigma)

#Format and plot the reuslts
sim_MH_single <- as.data.table(sim_MH_single)
setnames(sim_MH_single, c("t", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_single[, itteration := seq(1:n)]

q <- ggplot(data = sim_MH_single, aes(x = itteration) )
q <- q + geom_line(aes(y = lambda_0, colour = "lambda_0"))
q
q <- ggplot(data = sim_MH_single, aes(x = itteration) )
q <- q + geom_line(aes(y = lambda_1, colour = "lambda_1"))
q
q <- ggplot(data = sim_MH_single, aes(x = itteration) )
q <- q + geom_line(aes(y = beta, colour = "beta"))
q
q <- ggplot(data = sim_MH_single, aes(x = itteration) )
q <- q + geom_line(aes(y = t, colour = "t"))
q

# Check if values seem reasonable
summary(sim_MH_single)
3.655*(1891- t_0) + 0.62*(t_2-1891)

q <- ggplot(data = data, aes(x = date, y = cum_event))
q <- q + geom_line()
q

##Looking at the acf of t

burnin <- 1000

sim_t <- sim_MH_single$t[burnin:n]

acf(sim_t)
