#implemetn MCMC


# Want Q to be a random walk


# f_cond_beta <- function(beta, lambda_0, lambda_1){
#   return(exp(-1/beta * (lambda_0 + lambda_1 + 1))*(1/beta)^5)
# }

MH_alg <- function(n,data, t_0,t_2, t, lambda_0,lambda_1, beta, sigma){

  y_0 = sum(data[date <t]$event)
  y_1 = sum(data[date>= t]$event)
  y_0
  y_1
  
  t = t
  t
  beta = beta
  results <- matrix(nrow = n, ncol = 5)
  results[1, 1:5] <- c(t, lambda_0, lambda_1, beta, 1)
  for(i in 2:n){
    print(i)
    results[2:10, 1:4]
    #print("t")
    #print(results[i-1,1])
    t_old <- results[i-1, 1]
    
    y_0_o = sum(data[date <t_old]$event)
    y_1_o = sum(data[date>= t_old]$event)
    #print(i)
    #gibs step
    
    beta <- results[i-1, 4]
    beta
    
    lambda_0 <- rgamma(1, y_0 + 2, scale = 1/(t_old - t_0 + 1/beta))
    lambda_1 <- rgamma(1, y_1 + 2, scale = 1/(t_2 - t_old + 1/beta))
    
    
    lambda_0
    lambda_1
    beta_tilde <- rgamma(1, 6, scale = 1/(lambda_0 + lambda_1 + 1)) #is this okay to d0?
    beta_tilde
    beta <- 1/beta_tilde
    beta
    results[i, 2:4] <- c( lambda_0, lambda_1, beta)
    
    #MH step
    
    t_old <- results[i-1, 1]
    t_new <- rnorm(1, t_old, sigma)
    t_new
    print(t_new)
    
    y_0_n = sum(data[date <t_new]$event)
    y_1_n = sum(data[date>= t_new]$event)
    #print(c("T_old", t_old, "T_new", t_new))
    ## acceptannce prob
    
    a_log <- (y_0_n +1)*log(lambda_0) + (y_1_n +1)*log(lambda_1) - t_new *(lambda_0 - lambda_1 ) - ((y_0_o +1)*log(lambda_0) + (y_1_o +1)*log(lambda_1) - t_old *(lambda_0 - lambda_1 ))
    
    a_log
    
#    a <-  (lambda_0^{y_0_n+1} *lambda_1^{y_1_n+1}*exp(-t_new *(lambda_0 - lambda_1 )))/(lambda_0^{y_0_o+1} *lambda_1^{y_1_o+1}*exp(-t_old *(lambda_0 - lambda_1 )))

 #   print(a)
    
    exp(a_log)
    
    
    u <- runif(1)
    u
    results[i, 5] <-  a
    if (u< exp(a_log)){
      results[i, 1] <-  t_new
    } else{
      results[i,1] <- t_old
      }
  }
  return (results)
}


library(boot)
library(data.table)
library(tidyr)
library(ggplot2)

data <- as.data.table(coal)
data[, event:= 1]
data[, cum_event := cumsum(event)]
data[, date := date- 1851]
q <- ggplot(data = data, aes(x = date, y = cum_event))
q <- q + geom_line()
q

t = 45

lambda_0 <- 10 #135/50
lambda_1 <- 5 #56/50
beta = 3

n = 10000
t_0 = 0
t_2 = 112 #maybe should be 1963? 

sigma = 3

sim_MH <- MH_alg(n,data, t_0, t_2, t, lambda_0,lambda_1, beta, sigma)

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
#only beta has a burn in period! kind of

summary(sim_MH)
2.92*45 + 0.92*(112-45)
#seems very reasonable!! 


##Looking at the acf of t

burnin <- 1000

sim_t <- sim_MH$t[burnin:n]

acf(sim_t)
