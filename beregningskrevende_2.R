#implemetn MCMC


# Want Q to be a random walk


f_cond_beta <- function(beta, lambda_0, lambda_1){
  return(exp(-1/beta * (lambda_0 + lambda_1 + 1))*(1/beta)^5)
}

MH_alg <- function(n,data, t_0,t_2, t, lambda_0,lambda_1, beta, sigma){
  
   t = 45
  # y_0 = sum(data[date <t]$event)
  # y_1 = sum(data[date>= t]$event)
  # lambda_0 <- y_0/t
  # lambda_1 <- y_1/(190-t)
  # beta = 1
  # 
  # n = 100
  # t_0 = 0 #1851
  # t_2 = 1963- 1851 #maybe should be 1963? 
  # sigma = 1
  
  
  y_0 = sum(data[date <t]$event)
  y_1 = sum(data[date>= t]$event)
  y_0
  y_1
  results <- matrix(nrow = n, ncol = 4)
  results[1, 1:4] <- c(t, lambda_0, lambda_1, beta)
  for(i in 2:n){
    results[2:10, 1:4]
    #print("t")
    #print(results[i-1,1])
    y_0 = sum(data[date <results[i-1,1]]$event)
    y_1 = sum(data[date>= results[i-1,1]]$event)
    #print(i)
    #gibs step
    #t <- rexp(1, 1/(lambda_0 - lambda_1)) #must be wrong! 
    lambda_0 <- rgamma(1, (y_0 + 1), scale= (1/(t - t_0 + 1/beta)))
    lambda_1 <- rgamma(1, (y_1 + 1), scale=(1/(t_2 - t + 1/beta)))
    lambda_0
    lambda_1

    results[i, 1:3] <- c(t, lambda_0, lambda_1)
    
    #MH step
    beta_old <- results[i-1, 4]
    beta_new <- rnorm(1, beta_old, sigma)
    beta_old
    beta_new
    ## acceptannce prob
    
    a <- min(1, f_cond_beta(beta_new, lambda_0 , lambda_1)/f_cond_beta(beta_old, lambda_0 , lambda_1))
    a
    u <- runif(1)
    u
    if (u< a){
      results[i, 4] <-  beta_new
    } else{
      results[i,4] <- beta_old
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
# data_aggregated <- data.table(
#   y_0 = sum(data[date <t]$event),
#   y_1 = sum(data[date>= t]$event)
# )

lambda_0 <- 10 #135/50
lambda_1 <- 5 #56/50
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


par(mfrow=c(3,1))
plot(x = seq(1:n), y = sim_MH[1:nrow(sim_MH) ,3], type="l", xlab="Iteration", ylab="x")
plot(trace2, type="l", xlab="Iteration", ylab="x")
plot(trace3, type="l", xlab="Iteration", ylab="x")

## show histograms

plot_dist <- function(trace, until, limits){
  
  if(until > length(trace)){
    stop("the trace is not long enough!")
  }
  
  my_x = 0:20
  my_y = dpois(my_x, lambda=10)
  tmp <- factor(trace[1:until], levels=0:20)
  tab <- table(tmp)
  plot(my_x, my_y, type="b", col=2, xlab="x", ylab="Density", 
       main=paste("Iteration", until), ylim=limits)
  lines(0:20, tab/length(tmp), type="h", cex=1.2)
}

trace1 <- simMC_pois10_mcmc(5000, x_start = 0)