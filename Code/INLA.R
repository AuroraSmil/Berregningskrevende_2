library(ggplot2)
library(data.table)
library(Matrix)
library(MASS)
library(mvtnorm)
library(latex2exp)
library(tikzDevice)

# Importing data
data <- as.data.table(read.delim("code/Gaussiandata.txt", header = FALSE))
data[, t := seq(1:nrow(data))]
setnames(data, "V1", "value")
head(data)


# Plot the data
q <- ggplot(data, aes(t, value))
q <- q + geom_point()
q <- q + labs(title = "Gaussian data", x = "t",  y = "Value")
q

ggsave(
  "gaussian_data.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

#Define T = number of events
T = 20

#Make precision matrix
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

# Define full conditional for theta
pi_theta_cond <- function(T, eta, Q){
  return ( rgamma(1, shape = ( (T/2 )), scale= 1/(1 + (t(eta)%*% Q%*% eta)*0.5)))
}

#Define full conditional for eta
pi_eta_full_cond <- function(Q, I ,  y){
  mu <- solve((I + Q))%*%y
  sigma <- solve(I + Q)
  retval <- mvrnorm(n = 1, mu, sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  return (retval)
}

# Define the gibbs algorithm
Gibbs_alg <- function(n ,data, theta_0, eta_0){
  T = 20
  Q <-as.matrix( make_Q(20))
  I <- diag(nrow = 20, ncol = 20)
  y <- data$value
  
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

#Run algorithm
n = 100000

#Define initial values for theta and eta
theta <-  0
## Collect Q, I and y
Q <-as.matrix( make_Q(20))
I <- diag(nrow = 20, ncol = 20)
y <- data$value
eta <-  pi_eta_full_cond(Q*(theta), I, y)

posterior_dist_MC <- Gibbs_alg(n, data, eta_0 = eta, theta_0 = theta)

q <- ggplot(data = as.data.table(posterior_dist_MC[1:25000,]), aes(x = 1:25000) )
q <- q + geom_line(aes(y = V1))
q <- q + ylab(unname(TeX(c("$\\theta$")))) +
  ggtitle(unname(TeX("Traceplot for MCMC samples of theta $\\theta$"))) +
  xlab("Iteration")
q

ggsave(
  "trace_theta_mcmc.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

# Discard burn in period
posterior_dist_MC_usefull <- posterior_dist_MC[10000:nrow(posterior_dist_MC), ]
posterior_dist_MC_usefull <- as.data.table(posterior_dist_MC_usefull)
posterior_dist_MC_usefull[, itteration:= 1:nrow(posterior_dist_MC_usefull)]

summary(posterior_dist_MC_usefull)

q <- ggplot(posterior_dist_MC_usefull, aes(x = V1))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + geom_vline(xintercept = 1.797) + geom_text(x=2.2, y=0.5, label="Mean")
q <- q + geom_vline(xintercept = 1.6136) + geom_text(x=1.8, y=0.6, label="Mode")
q <- q + xlab(unname(TeX(c("$\\theta$"))))
q <- q + ylab("Density")
q <- q + ggtitle(unname(TeX("Posterior density of $\\theta$")))
q


ggsave(
  "post_theta_mcmc.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)



posterior_dist_MC_usefull.long = melt(posterior_dist_MC_usefull, id.vars = c("V1", "itteration"),
             measure.vars = c(2:21))



post_eta_mcmc_aggregated <- posterior_dist_MC_usefull.long[, .(
  mean = mean(value),
  q_lower = quantile(value, probs = 0.025), 
  q_upper = quantile(value, probs = 0.975)),
  keyby = .(variable)]

post_eta_mcmc_aggregated[,variable:= seq(1,20,1)]

q <- ggplot(post_eta_mcmc_aggregated, aes(x = as.factor(variable), y = mean))
q <- q + geom_point()+
  geom_errorbar(aes(ymin=q_lower, ymax=q_upper), width=.2,
                position=position_dodge(0.05))
q <- q + xlab(unname(TeX(c("$\\eta$"))))
q <- q + ylab("value")
q <- q + ggtitle("Mean and 95 % credible interval for the smoothing parameter")
q

ggsave(
  "post_eta_mcmc.pdf",
  path ="/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

#INLA scheme #################

## First inla step

# pi(theta|y)
pi_theta_given_y <- function(theta, Q, I, y, eta){
  eta_theta <- theta^((T-2)/2)*exp(-theta - 0.5*theta*t(eta)%*%Q%*%eta)
  likelihood <- dmvnorm(x = y, mean = eta, sigma = diag(1, nrow= T, ncol= T))
  mu_2 <- solve((Q*theta+ I))%*%y
  sigma_2 <- solve(Q*theta+ I)
  pi_eta_full <- dmvnorm(eta, mu_2, sigma_2)
  retval <- c(eta_theta*likelihood/pi_eta_full, pi_eta_full)
  return(retval)
}

Q <-as.matrix( make_Q(20))
I <- diag(nrow = 20, ncol = 20)
y <- data$value

#Pick random eta
eta <- pi_eta_full_cond(Q*1, I, y)
n = 100
theta_grid <- (seq(0, 6, length.out = n))

post_theta <- (lapply(seq(0.0000, 6, length.out = n), pi_theta_given_y, Q, I, y, eta))

require(purrr)
data_post <- map_df(post_theta, ~as.data.frame(t(.)))
data_post <- as.data.table(data_post)
setnames(data_post, c("V1", "V2"), c("post_theta", "eta_full"))

data_post[, x := theta_grid]

q <- ggplot(data = as.data.table(data_post)[1:nrow(data_post)])
q <- q + geom_point(aes(x =x, y =post_theta ))
q <- q + xlab(unname(TeX(c("$\\theta$"))))
q <- q + ylab("Density") + ggtitle(unname(TeX(c("Posterior density for $\\theta$"))))
q


ggsave(
  "post_theta_inla.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

###### last inla step
theta_grid <- (seq(0, 6, length.out = n))
eta_grid <- seq(-3.5,3.5, length.out  = 100)
Q <- make_Q(20)
I <- diag(1, 20, 20)
results <- vector("list", length = length(theta_grid))
j = 10

pdf(file="/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images/eta_10.pdf")
plot(1, type="n", xlab=unname(TeX(c("$\\eta_{10}$"))), ylab="", xlim=c(-3.5, 3.5), ylim=c(0, 2))

for (i in seq_along(results)){
  sigma <- solve(Q*theta_grid[i] + I)
  sigma_i <- sigma[j,j]
  mu <- sigma%*%y
  mu_i <- mu[j]
  density <- dnorm(eta_grid, mu_i, sigma_i)
  points(eta_grid, density)
  results[[i]] <- as.data.table(density)* post_theta[[i]][1] 
}

dev.off()

eta <- as.data.table(results)[, density_total := rowSums(.SD)][]
eta[, grid:= eta_grid]
plot(eta_grid, eta$density_total)

q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V11))
q <- q + geom_histogram(aes(y = ..density.., colour = "MCMC"), bins = 100)
q <- q + geom_point(data = eta, aes(x = eta_grid, y = density_total*0.9*10^7, colour = "INLA_manualy"))
q <- q + xlab(unname(TeX(c("$\\eta_{10}$"))))
q <- q + ylab("Density") + ggtitle(unname(TeX(c("Posterior density for $\\eta_{10}$"))))
q

ggsave(
  "post_eta_inla.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = eta, aes(x = eta_grid, y = density_total))
q <- q + geom_point()
q


##### INLA
library(INLA)
data <- as.data.table(read.delim("code/Gaussiandata.txt", header = FALSE))
data[, t := seq(1:nrow(data))]
setnames(data, "V1", "value")

formula  = value ~ -1 + f(t, model = "rw2", hyper = list(prec = list(prior = "loggamma", param = c(1,1))))
result  = inla(formula = formula, family = "gaussian",
               data = data, 
               verbose = TRUE)

precision <- result$marginals.hyperpar$`Precision for t`
eta_10_inla<- result$marginals.random$t$index.10


q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V1))
q <- q + geom_line(data = as.data.table(precision), aes(x =x, y = y , colour = "INLA"))
q <- q + xlim(0, 6) + xlab(unname(TeX(c("$\\eta_{10}$")))) + ylab("Density")
q


ggsave(
  "R_inla_theta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V1))
q <- q + geom_histogram(aes(y = ..density.., colour = "MCMC"), bins = 100)
q <- q + geom_vline(xintercept = 1.75)
q <- q + geom_line(data = as.data.table(precision), aes(x =x, y = y , colour = "INLA"))
q <- q + geom_point(data = data_post, aes(x =x, y = post_theta *0.25*10^ 9, colour = "INLA_manually"))
q <- q + geom_vline(xintercept = 0.66) + xlim(0, 6) + xlab(unname(TeX(c("$\\theta$")))) + ylab("Density")
q

ggsave(
  "theta_comparison.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V11))
q <- q + geom_point(data = eta, aes(x = eta_grid, y = density_total*0.9*10^7, colour = "INLA")) + 
   xlab("x") + ylab("eta")
q


ggsave(
  "R_inla_eta.pdf",
  path ="/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V11))
q <- q + geom_histogram(aes(y = ..density.., colour = "MCMC"), bins = 100)
q <- q + geom_point(data = eta, aes(x = eta_grid, y = density_total*0.9*10^7, colour = "INLA_manually"))
q <- q + geom_line(data = as.data.table(eta_10_inla), aes(x =x, y = y , colour = "INLA")) +
  xlab(unname(TeX(c("$\\eta_{10}$")))) + ylab("Density")
q

ggsave(
  "smoothing_comparison.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


q <- ggplot(as.data.table(posterior_dist_MC), aes(x = V11))
q <- q + geom_histogram(aes(y = ..density.., colour = "MCMC"), bins = 100)
q <- q + geom_point(data = eta, aes(x = eta_grid, y = density_total*0.9*10^7, colour = "INLA_manually"))
q <- q + geom_line(data = as.data.table(eta_10_inla), aes(x =x + 0.65, y = y , colour = "INLA")) +
  xlab(unname(TeX(c("$\\eta_{10}$")))) + ylab("Density")
q

ggsave(
  "smoothing_comparison_shifted.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

post_eta_mcmc_aggregated[, group := "MCMC"]
post_eta_inla_summary <- as.data.table(result$summary.random)[, .(t.ID, t.mean, t.0.025quant, t.0.975quant)]
post_eta_inla_summary[, group := "INLA"]

post_eta_total <- rbindlist(list(post_eta_mcmc_aggregated, post_eta_inla_summary))

q <- ggplot(post_eta_total, aes(x = as.factor(variable), y = mean, group = group, colour = group))
q <- q + geom_point()+
  geom_errorbar(aes(ymin=q_lower, ymax=q_upper), width=.2,
                position=position_dodge(0.3))
q <- q + xlab(unname(TeX(c("$\\eta$"))))
q <- q + ylab("value")
q <- q + ggtitle("Mean and 95 % credible interval for the smoothing parameter")
q

ggsave(
  "smoothing_comparison_all eta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)
