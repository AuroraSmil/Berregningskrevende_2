
library(boot)
library(data.table)
library(tidyr)
library(ggplot2)

data <- as.data.table(coal)
data[, event:= 1]
data[, cum_event := cumsum(event)]
data[, date := date- 1851]

MH_blokking <- function(n, t_0, t_2, t, lambda_0, lambda_1, beta, sigma_t, sigma_beta, data){
  
  # n = 100
  # t <- 45
  # lambda_0 <- 3
  # lambda_1 <- 1.5
  # beta <- 1
  # sigma_t <- 3
  # sigma_beta <- 0.5
  # t_0 <- 0
  # t_2 <- 113
  
  results <- matrix(nrow = n, ncol = 5)
  results[1, 1:5] <- c(t, lambda_0, lambda_1, beta, 1)

  
  
  for (i in 2:n){
    print(i)
    if(i %% 2 == 0){
      t_old <- results[i-1, 1]
      t_new <- rnorm(1, t_old, sigma_t)
      
      y_0_n <-  sum(data[date <t_new]$event)
      y_1_n <- sum(data[date>= t_new]$event)
      y_0_o <-  sum(data[date <t_old]$event)
      y_1_o <- sum(data[date>= t_old]$event)
      beta <- results[i-1, 4]
      lambda_0_n <- rgamma(1, y_0_n + 2, scale = 1/(t_new - t_0 + 1/beta))
      lambda_1_n <- rgamma(1, y_1_n + 2, scale = 1/(t_2 - t_new + 1/beta))
      lambda_0_o <- results[i-1, 2]
      lambda_1_o <- results[i-1, 3]
      #Acceptance probability 
      
      ## Proposal density
      
      prop_nominator <- dnorm(t_old, mean = t_new, sigma_t, log = TRUE) + 
        dgamma(1, y_0_o + 2, scale = 1/(t_old - t_0 + 1/beta), log = TRUE) + 
        dgamma(1, y_1_o + 2, scale = 1/(t_2 - t_old + 1/beta), log = TRUE)
      
      prop_denominator <- dnorm(t_new, mean = t_old, sigma_t, log = TRUE) + 
        dgamma(1, y_0_n + 2, scale = 1/(t_new - t_0 + 1/beta), log = TRUE) + 
        dgamma(1, y_1_n + 2, scale = 1/(t_2 - t_new + 1/beta), log = TRUE)
      
      target_nominator <- lambda_0_n*(t_new- t_0) - lambda_1_n*(t_2- t_new) - 1/beta *(lambda_0_n + lambda_1_n) + (y_0_n + 1)*log(lambda_0_n) + (y_1_n + 1)*log(lambda_1_n)
      
      exp(-lambda_0_n*(t_new-t_0) - lambda_1_n*(t_2-t_new)- beta*(lambda_0 + lambda_1))*lambda_0^{y_0+ 1}*lambda_1^{y_1+1}
      
      target_denominator <- lambda_0_o*(t_old- t_0) - lambda_1_o*(t_2- t_old) - 1/beta *(lambda_0_o + lambda_1_o) + (y_0_o + 1)*log(lambda_0_o) + (y_1_o + 1)*log(lambda_1_o)
      
      acceptance_log <- prop_nominator + target_nominator - prop_denominator - target_denominator
      
      a_log <- min(0, acceptance_log)
      
      u <- runif(1)
      a <- exp(a_log)
      a
      if (u < a){
        results[i, 1:5] <- c(t_new, lambda_0_n, lambda_1_n, beta, a)  
      }else{
        results[i, 1:5] <- c(t_old, lambda_0_o, lambda_1_o, beta, a)  
      }
      
    }else{
        
      t <- results[i-1, 1]
      y_0 <-  sum(data[date <t]$event)
      y_1 <- sum(data[date>= t]$event)
      
      beta_old <- results[i-1, 4]
      beta_new <- rnorm(1, mean = beta_old, sigma_beta)

      lambda_0_o <- results[i-1, 2]
      lambda_1_o <- results[i-1, 3]
      
      
      if(beta_new<=0){
        results[i, 1:5] <- c(t, lambda_0_o, lambda_1_o, beta_old, 0)
      }else{
        lambda_0_n <- rgamma(1, y_0 + 2, scale = 1/(t - t_0 + 1/beta_new))
        lambda_1_n <- rgamma(1, y_1 + 2, scale = 1/(t_2 - t + 1/beta_old))
        
        #Acceptance probability 
        
        ## Proposal density
        
        prop_nominator <- dnorm(beta_old, mean = beta_new, sigma_beta, log = TRUE) + 
          dgamma(1, y_0 + 2, scale = 1/(t - t_0 + 1/beta_old), log = TRUE) + 
          dgamma(1, y_1 + 2, scale = 1/(t_2 - t + 1/beta_old), log = TRUE)
        
        prop_denominator <- dnorm(beta_new, mean = beta_old, sigma_beta, log = TRUE) + 
          dgamma(1, y_0 + 2, scale = 1/(t - t_0 + 1/beta_new), log = TRUE) + 
          dgamma(1, y_1 + 2, scale = 1/(t_2 - t + 1/beta_new), log = TRUE)
        
        target_nominator <- lambda_0_n*(t- t_0) - lambda_1_n*(t_2- t) - 1/beta_new *(lambda_0_n + lambda_1_n + 1) + 
          (y_0 + 1)*log(lambda_0_n) + (y_1 + 1)*log(lambda_1_n) - 5 *log(beta_new)
        
        target_denominator <- lambda_0_o*(t- t_0) - lambda_1_o*(t_2- t) - 1/beta_old *(lambda_0_o + lambda_1_o + 1) + 
          (y_0 + 1)*log(lambda_0_o) + (y_1 + 1)*log(lambda_1_o) - 5 *log(beta_old)
        
        acceptance_log <- prop_nominator + target_nominator - prop_denominator - target_denominator
        
        a_log <- min(0, acceptance_log)
        
        u <- runif(1)
        a <- exp(a_log)
        
        if (u < a){
          results[i, 1:5] <- c(t_new, lambda_0_n, lambda_1_n, beta_new, a)  
        }else{
          results[i, 1:5] <- c(t_old, lambda_0_o, lambda_1_o, beta_old, a)  
        }
      }
    }
  }
  return(results)
}


n = 100
t <- 45
lambda_0 <- 3
lambda_1 <- 1.5
beta <- 1
sigma_t <- 3
sigma_beta <- 0.1
t_0 <- 0
t_2 <- 113

sim_MH <- MH_blokking(n, t_0, t_2, t, lambda_0, lambda_1, beta, sigma_t, sigma_beta, data)

sim_MH <- as.data.table(sim_MH)

setnames(sim_MH, c("t", "lambda_0", "lambda_1", "beta", "a"))

sim_MH[, itteration := seq(1:n)]
sim_MH

q <- ggplot(data = sim_MH, aes(x = itteration) )
q <- q + geom_line(aes(y = lambda_0, colour = "lambda_0"))
q
q <- ggplot(data = sim_MH, aes(x = itteration) )
q <- q + geom_line(aes(y = lambda_1, colour = "lambda_1"))
q
q <- ggplot(data = sim_MH, aes(x = itteration) )
q <- q + geom_line(aes(y = beta, colour = "beta"))
q
q <- ggplot(data = sim_MH, aes(x = itteration) )
q <- q + geom_line(aes(y = t, colour = "t"))
q
