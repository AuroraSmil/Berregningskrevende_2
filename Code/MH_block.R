#MH block

library(boot)
library(data.table)
library(tidyr)
library(ggplot2)

data <- as.data.table(coal)
data[, event:= 1]
data[, cum_event := cumsum(event)]
data[, date := date- 1851]

f_cond_beta <- function(beta, lambda_0, lambda_1){
  return(exp(-1/beta * (lambda_0 + lambda_1 + 1))*(1/beta)^5)
}



f_target_1 <- function(lambda_0, lambda_1, t, t_0, t_2, beta, y_0, y_1){
  return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1))*lambda_0^y_0*lambda_1^y_1)
}

f_prop_1 <- function(lambda_0, lambda_1, t, t_0, t_2, beta, y_0, y_1){
  return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1))*lambda_0^y_0*lambda_1^y_1)
}

f_target_2 <- function(lambda_0, lambda_1, beta, t_0, t_2, t, y_0, y_1){
  return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1+1))*lambda_0^y_0*lambda_1^y_1)*(1/beta^5)
}

f_prop_2 <- function(lambda_0, lambda_1, t, t_0, t_2, beta, y_0, y_1){
  return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1))*lambda_0^y_0*lambda_1^y_1)
}

MH_alg <- function(n,data, t_0,t_2, t, lambda_0,lambda_1, beta, sigma){
  
  
  y_0 = sum(data[date <t]$event)
  y_1 = sum(data[date>= t]$event)
  y_0
  y_1
  results <- matrix(nrow = n, ncol = 4)
  results[1, 1:4] <- c(t, lambda_0, lambda_1, beta)
  for(i in 2:(n)){
    print(i)
    if(i %% 2 == 0){
      t_old <- results[i-1, 1]
      t_old
      t_new <- rnorm(1, t_old, sigma)
      t_new
      beta_old <- results[i-1, 4]
      #print("t")
      #print(results[i-1,1])
      y_0 = sum(data[date <results[i-1,1]]$event)
      y_1 = sum(data[date>= results[i-1,1]]$event)
      #print(i)
      #gibs step
      #t <- rexp(1, 1/(lambda_0 - lambda_1)) #must be wrong! 
      lambda_0 <- rgamma(1, (y_0 + 1), scale= (1/(t_new - t_0 + 1/beta_old)))
      print("scale")
      print(1/(t_new - t_0 + 1/beta_old))
      print("t")
      print(t_new)
      print("beta")
      print(beta_old)
      lambda_1 <- rgamma(1, (y_1 + 1), scale=(1/(t_2 - t_new + 1/beta_old)))
      print(lambda_0)
      print(lambda_1)
      
      target_ratio <-  f_target_2(lambda_0, lambda_1, beta_new, t_0, t_2, t_old, y_0, y_1)/f_target_2(results[i-1, 2], results[i-1, 3], beta_old, t_0, t_2, t_old, y_0, y_1)
      target_ratio
      prop_ratio <- f_prop_1(lambda_0, lambda_1, t_new, t_0, t_2, beta_old, y_0, y_1)/f_prop_1(results[i-1, 2], results[i-1, 3], t_old, t_0, t_2, beta_old, y_0, y_1)
      prop_ratio
      a <- target_ratio*prop_ratio
      a
      u <- runif(1)
      u
      if (u< a){
        results[i, 1:4] <- c(t_new, lambda_0, lambda_1, beta_old)
      } else{
        results[i,1:4] <- c(t_old, results[i-1, 2], results[i-1, 3], beta_old)
      }
    }
    else{
      t_old <- results[i-1, 1]
      t_old
      beta_old <- results[i-1, 4]
      beta_new <- rnorm(1, beta_old, sigma)
      
      lambda_0 <- rgamma(1, (y_0 + 1), scale= (1/(t_old - t_0 + 1/beta_new)))
      lambda_1 <- rgamma(1, (y_1 + 1), scale=(1/(t_2 - t_old + 1/beta_new)))
      lambda_0
      lambda_1
      print(lambda_0)
      print(lambda_1)
      
      target_ratio <-  f_target_2(lambda_0, lambda_1, t_new, t_0, t_2, beta_old, y_0, y_1)/f_target_1(results[i-1, 2], results[i-1, 3], t_old, t_0, t_2, beta_old, y_0, y_1)
      target_ratio
      prop_ratio <- f_prop_2(lambda_0, lambda_1, t_old, t_0, t_2, beta_new, y_0, y_1)/f_prop_2(results[i-1, 2], results[i-1, 3], t_old, t_0, t_2, beta_old, y_0, y_1)
      prop_ratio
      a <- target_ratio*prop_ratio
      a
      u <- runif(1)
      u
      if (u< a){
        results[i, 1:4] <- c(t_new, lambda_0, lambda_1, beta_old)
      } else{
        results[i,1:4] <- c(t_old, results[i-1, 2], results[i-1, 3], beta_old)
      }
      
    }
  }
  return (results)
}




t = 45

lambda_0 <- 135/50
lambda_1 <- 56/50
beta = 10

n = 10000
t_0 = 0
t_2 = 112 #maybe should be 1963? 
sigma = 3

sim_MH <- MH_alg(n,data, t_0, t_2, t, lambda_0,lambda_1, beta, sigma)

sim_MH <- as.data.table(sim_MH)
sim_MH

setnames(sim_MH, c("t", "lambda_0", "lambda_1", "beta"))

sim_MH[, itteration := seq(1:n)]
sim_MH

q <- ggplot(data = sim_MH, aes(x = itteration) )
q <- q + geom_line(aes(y = lambda_0, colour = "lambda_0"))
q
q <- q + geom_line(aes(y = lambda_1, colour = "lambda_1"))
q
q <- q + geom_line(aes(y = lambda_1, colour = "beta"))
q

#only beta has a burn in period! kind of

summary(sim_MH)
2.92*45 + 0.92*(112-45)