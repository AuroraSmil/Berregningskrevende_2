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
  M[T,T]  <- 1
  M[1, 2] <- -2
  M[2, 1] <- -2
  M[T, T-1] <- -2
  M[T-1, T] <- -2
  return(M)
}

pi_theta_cond <- function(T){
  return ( rgamma(T/2, scale= 2/3))
}

pi_theta_cond <- function(T, eta, Q){
  return ( rgamma(1, shape = ( (T/2-2)), scale= 1/(t(eta)%*% Q%*% eta)))
}

Q <-as.matrix( make_Q(20))
I <- diag(nrow = 20, ncol = 20)
y <- data$value

pi_eta_full_cond <- function(Q, I ,  y){
  mu <- t(y) %*% solve((Q+ I))
  sigma <- solve(Q+ I)
  retval <- mvrnorm(n = 1, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  return (retval)
}


Gibbs_alg <- function(n ,data){
  Q <-as.matrix( make_Q(20))
  I <- diag(nrow = 20, ncol = 20)
  y <- data$value
  
  theta_0 <- 1# pi_theta_cond(T)
  eta_0 <- pi_eta_full_cond(Q*theta_0, I, y)
  
  results <- matrix(nrow = n, ncol = 21)
  results[1, 1:21] <- c(theta_0, eta_0)
  
  for(i in 2:n){
    theta <- pi_theta_cond(T, eta, Q)
    eta <- pi_eta_full_cond(Q, I, y)
    results[i, 1:21] <- c(theta, eta)
  }
  return (results)
}

n = 10000

posterior_dist_MC <- Gibbs_alg(n, data)
summary(posterior_dist_MC)

q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V1))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + geom_vline(xintercept = 0.66)
q

q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V11))
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
  print("next")
  print(theta)
  d_theta <- dgamma(theta, 1,1)
  print(d_theta)
  
  sigma_1 <- solve(Q*theta/2)
  mu_1 <- as.matrix(rep(0, 20))
  eta_given_theta <- dmvnorm(x = eta, mean = mu_1, sigma = sigma_1)
  print(eta_given_theta)
  likelihood <- dmvnorm(x = t(y), mean = eta, sigma = diag(1, nrow= T, ncol= T))
  print(likelihood)
  mu_2 <- t(y) %*% solve((Q*theta/2+ I))
  sigma_2 <- solve(Q*theta/2+ I)
  pi_eta_full <- dmvnorm(eta, mu_2, sigma_2)
  print(pi_eta_full)
  retval <- c(((d_theta*eta_given_theta*likelihood)/pi_eta_full), pi_eta_full)
  return(retval)
}


Q <-as.matrix( make_Q(20))
I <- diag(nrow = 20, ncol = 20)
y <- data$value
eta <- pi_eta_full_cond(Q*1, I, y)
n = 100
theta_grid <- (seq(0.00001, 2, length.out = n))
#theta =  pi_theta_given_y(theta_grid[1], Q, I, y, eta)

post_theta <- (lapply(seq(0.00001, 6, length.out = n), pi_theta_given_y, Q, I, y, eta))


require(purrr)
data_post <- map_df(post_theta, ~as.data.frame(t(.)))
data_post <- as.data.table(data_post)
setnames(data_post, c("V1", "V2"), c("post_theta", "eta_full"))

data_post[, x := theta_grid]
#data_post[,   theta_pri := dgamma(x,1,1)]


q <- ggplot(data = as.data.table(data_post))
q <- q + geom_point(aes(x =x, y = post_theta ))
q

###### last inla step


##### INLA
library(INLA)
data <- as.data.table(read.delim("code/Gaussiandata.txt", header = FALSE))
data
data[, t := seq(1:nrow(data))]
data

setnames(data, "V1", "value")
data

formula  = value ~ 1 + f(t, model = "rw2", hyper = list(prec = list(prior = "loggamma", param = c(1,1))))
result  = inla(formula = formula, family = "gaussian",
               data = data, 
               verbose = TRUE)

precision <- result$marginals.hyperpar$`Precision for the Gaussian observations`
theta <- inla.tmarginal(function(x) 1/sqrt((x)), precision)


q <- ggplot(data = as.data.table(theta))
q <- q + geom_line(aes(x =x, y = y ))
q
