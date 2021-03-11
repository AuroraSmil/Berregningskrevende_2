library(ggplot2)
library(data.table)
library(Matrix)
library(MASS)
library(mvtnorm)

#importing data
data <- as.data.table(read.delim("code/Gaussiandata.txt", header = FALSE))
data
data[, t := seq(1:nrow(data))]
data

setnames(data, "V1", "value")
data

q <- ggplot(data, aes(t, value))
q <- q + geom_point()
q

q <- ggplot(data, aes(x = value))
q <- q + geom_histogram(aes(y = ..density..), bins = 50 )
q

T = 20

make_Q <- function(T){
  if (T <3){
    return("Dimention to low")
  }
  M <- bandSparse(T, T, #dimensions
                  (-2):2, #band, diagonal is number 0
                  list(rep(1, T-2), 
                       rep(-4, T-1),
                       rep(6, T), 
                       rep(-4, T-1),
                       rep(1, T-2)))
  M[1,1]  <- 1
  M[2,2]  <- 5
  M[T,T]  <- 1
  M[T-1,T-1]  <- 5
  M[1, 2] <- -2
  M[2, 1] <- -2
  M[T, T-1] <- -2
  M[T-1, T] <- -2
  return(M)
}

# pi_theta_cond <- function(T){
#   return ( rgamma(T/2, scale= 2/3))
# }

pi_theta_cond <- function(T, eta, Q){
  return ( rgamma(1, shape = ( (T/2 )), scale= 1/(1 + (t(eta)%*% Q%*% eta)*0.5)))
} # missing theta and

Q <-as.matrix( make_Q(20))
I <- diag(nrow = 20, ncol = 20)
y <- data$value

pi_eta_full_cond <- function(Q, I ,  y){
  mu <- solve((I + Q))%*%y
  sigma <- solve(I + Q)
  retval <- mvrnorm(n = 1, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  return (retval)
}


Gibbs_alg <- function(n ,data){
  T = 20
  Q <-as.matrix( make_Q(20))
  I <- diag(nrow = 20, ncol = 20)
  y <- data$value
  
  theta <-  0#1 #pi_theta_cond(T)
  eta <-  pi_eta_full_cond(Q*(theta), I, y)
  
  results <- matrix(nrow = n, ncol = 21)
  results[1, 1:21] <- c(theta, eta)
  
  for(i in 2:n){
    theta <- pi_theta_cond(T, results[i-1, 2:21], Q)
    theta
    eta <-  pi_eta_full_cond(Q*(theta), I, y)
    results[i, 1:21] <- c(theta, eta)
  }
  return (results)
}

n = 100000

posterior_dist_MC <- Gibbs_alg(n, data)
summary(posterior_dist_MC)


posterior_dist_MC_usefull <- posterior_dist_MC[25000:nrow(posterior_dist_MC), ]

post_data_MC <- as.data.table(posterior_dist_MC)
post_data_MC[, itteration:= 1:nrow(post_data_MC)]

q <- ggplot(data = post_data_MC, aes(x = itteration) )
q <- q + geom_line(aes(y = V1))
q

q <- ggplot(post_data_MC, aes(x = V1))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + geom_vline(xintercept = 1.75)
q


q <- ggplot(post_data_MC, aes(x = V11))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q

pi_eta_given_theta<- function(Q, I ,  y){
  mu <- as.matrix(rep(0, 20))
  sigma <- solve(Q)
  retval <- mvrnorm(n = 1, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  return (retval)
}


#INLA scheme #################

#just pick a value for eta? 



pi_theta_given_y <- function(theta, Q, I, y, eta){
  #d_theta <- dgamma(theta, 1,1)
  #print(d_theta)

  #sigma_1 <- solve(Q*theta + diag(0.0000001, 20, 20)) #small pertiibaitn to awoid singularity
  #mu_1 <- as.matrix(rep(0, 20))
  #eta_given_theta <- dmvnorm(x = eta, mean = mu_1, sigma = sigma_1)
  
  #eta_theta <- dgamma(theta, 1, 1/(1 + 0.5*t(eta)%*%Q%*%eta))
  eta_theta <- theta^((T-2)/2)*exp(-theta - 0.5*theta*t(eta)%*%Q%*%eta)
  #print(eta_given_theta)
  likelihood <- dmvnorm(x = y, mean = eta, sigma = diag(1, nrow= T, ncol= T))
  print(likelihood)
  mu_2 <- solve((Q*theta+ I))%*%y
  sigma_2 <- solve(Q*theta+ I)
  pi_eta_full <- dmvnorm(eta, mu_2, sigma_2)
  print(pi_eta_full)
  #retval <- c(((d_theta*eta_given_theta*likelihood)/pi_eta_full), pi_eta_full)
  retval <- c(eta_theta*likelihood/pi_eta_full, pi_eta_full)
  return(retval)
}

mu <- solve((I + Q))%*%y
sigma <- solve(I + Q)
retval <- mvrnorm(n = 1, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

Q <-as.matrix( make_Q(20))
I <- diag(nrow = 20, ncol = 20)
y <- data$value
eta <- pi_eta_full_cond(Q*1, I, y)
n = 100
theta_grid <- (seq(0, 6, length.out = n))
#theta =  pi_theta_given_y(theta_grid[1], Q, I, y, eta)

post_theta <- (lapply(seq(0.0000, 6, length.out = n), pi_theta_given_y, Q, I, y, eta))


require(purrr)
data_post <- map_df(post_theta, ~as.data.frame(t(.)))
data_post <- as.data.table(data_post)
setnames(data_post, c("V1", "V2"), c("post_theta", "eta_full"))

data_post[, x := theta_grid]
#data_post[,   theta_pri := dgamma(x,1,1)]


q <- ggplot(data = as.data.table(data_post)[10:nrow(data_post)])
q <- q + geom_point(aes(x =x, y =post_theta ))
q <- q + geom_vline(xintercept = 1.75)
q

###### last inla step

theta_grid
eta_grid <- seq(-2.5,2.5, length.out  = 100)
Q <- make_Q(20)
I <- diag(1, 20, 20)
results <- vector("list", length = length(theta_grid))
j = 10
plot(1, type="n", xlab="", ylab="", xlim=c(-5, 5), ylim=c(0, 2))
for (i in seq_along(results)){
  sigma <- solve(Q*theta_grid[i] + I)
  sigma_i <- sigma[j,j]
  mu <- sigma%*%y
  mu_i <- mu[j]
  #print(mu_i)
  print(sigma_i)
  density <- dnorm(eta_grid, mu_i, sigma_i)
  points(eta_grid, density)
  results[[i]] <- as.data.table(density)*post_theta[i] #her er feilen den ganger inn to tall
}

(as.data.table(results))

eta <- as.data.table(results)[, density_total := rowSums(.SD)][]

eta[, grid:= eta_grid]
plot(eta_grid, eta$density_total)

q <- ggplot(data = eta, aes(x = eta_grid, y = density_total))
q <- q + geom_point()
q



##### INLA
library(INLA)
data <- as.data.table(read.delim("code/Gaussiandata.txt", header = FALSE))
data
data[, t := seq(1:nrow(data))]
data

setnames(data, "V1", "value")
data

formula  = value ~ -1 + f(t, model = "rw2", hyper = list(prec = list(prior = "loggamma", param = c(1,1))))
result  = inla(formula = formula, family = "gaussian",
               data = data, 
               verbose = TRUE)


precision <- result$marginals.hyperpar$`Precision for t`




q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V1))
q <- q + geom_histogram(aes(y = ..density.., colour = "MCMC"), bins = 100)
q <- q + geom_vline(xintercept = 1.75)
q <- q + geom_line(data = as.data.table(precision), aes(x =x, y = y , colour = "INLA"))
q <- q + geom_point(data = data_post, aes(x =x, y = post_theta *0.25*10^ 9, colour = "INLA_malually"))
q <- q + geom_vline(xintercept = 0.66) + xlim(0, 6)
q

