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
  if(t<t_0){
    return (0)
  }
  if(t> t_2){
    return (0)
  }
  return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1))*lambda_0^y_0*lambda_1^y_1)
}

f_prop_1 <- function(lambda_0, lambda_1, t_given, t_cur, t_0, t_2, beta, y_0, y_1){
  if(t<t_0){
    retval <- 0
  }
  if(t> t_2){
   retval <- 0
  }
  #return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1))*lambda_0^y_0*lambda_1^y_1)
  
  
  retval <- dgamma(lambda_0, (y_0 ), scale= (1/(t_cur - t_0 + 1/beta)))*dgamma(lambda_1, (y_1 ), scale= (1/(t_2 - t_cur + 1/beta)))*dnorm(t, t_given,1)
  print("retval_prop")
  print(retval)
}

f_target_2 <- function(lambda_0, lambda_1, beta, t_0, t_2, t, y_0, y_1){
  if(beta<0){
    return (0)
  }
  print("beta_target")
  print(beta)
  return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1+1))*lambda_0^y_0*lambda_1^y_1)*(1/beta^5)
}

f_prop_2 <- function(lambda_0, lambda_1, t, t_0, t_2, beta_cur, beta_given, y_0, y_1){
  if(beta_cur<0){
    retval <- 0
  } else if (beta_given<0){
    retval <- 0
  } else{
    print("beta_prop")
    print(beta_cur)
    #return(exp(-lambda_0*(t-t_0) - lambda_1*(t_2-t)- beta*(lambda_0 + lambda_1))*lambda_0^y_0*lambda_1^y_1)
    # print((1/(t - t_0 + 1/beta_cur)))
    # print(1/(t_2 - t + 1/beta_cur))
    retval <- dgamma(lambda_0, (y_0 ), scale= (1/(t - t_0 + 1/beta_cur)))*dgamma(lambda_1, (y_1 ), scale= (1/(t_2 - t + 1/beta_cur)))*dnorm(beta_cur, beta_given,1)
  }

  print("retval prop_ratio")
  print(retval)
  return (retval)
}

MH_alg <- function(n,data, t_0,t_2, t, lambda_0,lambda_1, beta, sigma_t, sigma_beta){
  
  
  y_0 = sum(data[date <t]$event)
  y_1 = sum(data[date>= t]$event)
  y_0
  y_1
  results <- matrix(nrow = n, ncol = 4)
  results[1, 1:4] <- c(t, lambda_0, lambda_1, beta)
  for(i in 2:(n)){
    print(i)
    if(i %% 2 == 0){
      beta <- results[i-1, 4]
      beta
      t_old <- results[i-1, 1]
      t_old
      t_new <- rnorm(1, t_old, sigma_t)
      print("t_new")
      print(t_new)
      if(t_new <t_0 | t_new > t_2){
        results[i,1:4] <- c(t_old, lambda_0_old, lambda_1_old, beta)
      }else{
     
        #print("t")
        #print(results[i-1,1])
        y_0_old = sum(data[date <results[i-1,1]]$event)
        y_1_old = sum(data[date>= results[i-1,1]]$event)
        y_0_n = sum(data[date <t_new]$event)
        y_1_n = sum(data[date>= t_new]$event)
        #print(i)
        #gibs step
        #t <- rexp(1, 1/(lambda_0 - lambda_1)) #must be wrong! 
        lambda_0_old <- results[i-1, 2]
        lambda_1_old <- results[i-1, 3]
        lambda_0 <- rgamma(1, (y_0 ), scale= (1/(t_new - t_0 + 1/beta_old)))
         print("scale")
         print(1/(t_new - t_0 + 1/beta_old))
        # print("t")
        # print(t_new)
         print("beta")
         print(beta_old)
        lambda_1 <- rgamma(1, (y_1), scale=(1/(t_2 - t_new + 1/beta_old)))
        # print(lambda_0)
        # print(lambda_1)
        
        target_ratio <-  f_target_1(lambda_0, lambda_1, t_new, t_0, t_2, beta, y_0_n, y_1_n)/
          f_target_1(lambda_0_old, lambda_1_old, t_old, t_0, t_2, beta_old, y_0_old, y_1_old)
        target_ratio
        prop_ratio <- f_prop_1(lambda_0_old, lambda_1_old, t_cur = t_old, t_given = t_new,  t_0, t_2, beta, y_0_old, y_1_old)/
          f_prop_1(lambda_0, lambda_1, t_new, t_given = t_old, t_0, t_2, beta, y_0, y_1)
        prop_ratio
        a <- target_ratio*prop_ratio
        a
        a <- min(1, a)
        a
        print("a")
        print(a)
        u <- runif(1)
        u
        if (u< a){
          results[i, 1:4] <- c(t_new, lambda_0, lambda_1, beta)
        } else{
          results[i,1:4] <- c(t_old, lambda_0_old, lambda_1_old, beta)
        }
      }
    }
    else{
      t <- results[i-1, 1]
      t
      beta_old <- results[i-1, 4]
      beta_new <-  rgamma(1, 6, scale = 1/(lambda_0 + lambda_1 + 1)) #rnorm(1, beta_old, sigma_beta)
      y_0 = sum(data[date <results[i-1,1]]$event)
      y_1 = sum(data[date>= results[i-1,1]]$event)
      if(beta_new< 0){
        results[i,1:4] <- c(t, lambda_0_old, lambda_1_old, beta_old)
      }else {
        print("t")
        print(t)
        print("beta_n")
        print(beta_new)
        print("beta_o")
        print(beta_old)
        
        lambda_0_n <- rgamma(1, (y_0 ), scale= (1/(t_old - t_0 + 1/beta_new)))
        lambda_1_n <- rgamma(1, (y_1), scale=(1/(t_2 - t_old + 1/beta_new)))
        lambda_0_n
        lambda_1_n
        # print(lambda_0)
        # print(lambda_1)
        
        lambda_0_old <- results[i-1, 2]
        lambda_1_old <- results[i-1, 3]
        
        target_ratio <-  f_target_2(lambda_0, lambda_1, beta_new, t_0, t_2, t, y_0, y_1)/
          f_target_2(lambda_0_old, lambda_1_old, beta_old, t_0, t_2, t, y_0, y_1)
        target_ratio
        # prop_ratio <- f_prop_2(lambda_0_old, lambda_1_old, beta_cur = beta_old, beta_given = beta_new, t = t, t_0 = t_0, t_2 = t_2, y_0, y_1)/
        #   f_prop_2(lambda_0, lambda_1, beta_cur = beta_new, beta_given = beta_old, t_0 = t_0, t_2 = t_2, t = t, y_0, y_1)
        # prop_ratio
        
        prop_ratio <- (dgamma(lambda_0_old, (y_0 ), scale= (1/(t - t_0 + 1/beta_old)))*
          dgamma(lambda_1_old, (y_1), scale= (1/(t_2 - t + 1/beta_old)))*
          dgamma(beta_old,  6, scale = 1/(lambda_0_old + lambda_1_old + 1)))/
          (dgamma(lambda_0_n, (y_0 ), scale= (1/(t - t_0 + 1/beta_new)))*
          dgamma(lambda_1_n, (y_1 ), scale= (1/(t_2 - t+ 1/beta_new)))*
          dgamma(beta_new,  6, scale = 1/(lambda_0_n + lambda_1_n + 1)))
        
        
        #print(target_ratio)
        print("Prop_ratio")
        print(prop_ratio)
        a <- min(1, target_ratio*prop_ratio)
        a
        print("a")
        print(a)
        u <- runif(1)
        u
        if (u< a){
          results[i, 1:4] <- c(t, lambda_0, lambda_1, beta_new)
        } else{
          results[i,1:4] <- c(t, lambda_0_old, lambda_1_old, beta_old)
        }
      }
    }
  }
  return (results)
}




t = 45

lambda_0 <- 3#135/50
lambda_1 <- 3# 56/50
beta = 1

n = 10000
t_0 = 0
t_2 = 112 #maybe should be 1963? 
sigma_t = 0.5
sigma_beta = 0.5

sim_MH <- MH_alg(n,data, t_0, t_2, t, lambda_0,lambda_1, beta, sigma_t = sigma_t, sigma_beta = sigma_beta)

sim_MH <- as.data.table(sim_MH)

setnames(sim_MH, c("t", "lambda_0", "lambda_1", "beta"))

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

#only beta has a burn in period! kind of

summary(sim_MH)
2.92*45 + 0.92*(112-45)