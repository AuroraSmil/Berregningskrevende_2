
library(ggplot2)
library(data.table)
library(Matrix)
library(MASS)
library(mvtnorm)
library(latex2exp)
library(tikzDevice)

## SINGLE SITE MCMC 

MH_single_alg <- function(n,data, t_0,t_2, t, lambda_0,lambda_1, beta, sigma){

  # Find number of deaths in each interval given initial t
  # y_0 = sum(data[date <t]$event)
  # y_1 = sum(data[date>= t]$event)
  y_0 = sum(data[date <t]$event) -1
  y_1 = sum(data[date>= t]$event) -1
  
  # Construct result matrix
  results <- matrix(nrow = n, ncol = 5)
  results[1, 1:5] <- c(t, lambda_0, lambda_1, beta, 1)
  for(i in 2:n){
    #print(i)
    # Collect old values for t and beta and compute corresponding 
    # y_0 and y_1 values
    t_old <- results[i-1, 1]
    
    y_0_o = sum(data[date <t_old]$event) -1
    y_1_o = sum(data[date>= t_old]$event) -1
    
    beta <- results[i-1, 4]
    
    # Preform Gibbs step
    lambda_0 <- rgamma(1, y_0_o + 2, scale = 1/(t_old - t_0 + 1/beta))
    lambda_1 <- rgamma(1, y_1_o + 2, scale = 1/(t_2 - t_old + 1/beta))

    beta_tilde <- rgamma(1, 6, scale = 1/(lambda_0 + lambda_1 + 1)) 
    beta <- 1/beta_tilde
    
    results[i, 2:4] <- c( lambda_0, lambda_1, beta)
    
    # Preform MH step with the random walk of order one 
    # as the proposal density
    t_new <- rnorm(1, t_old, sigma)
    if(t_new <= t_0){
      results[i,1] <- t_old
    } else if (t_new >= t_2){
      results[i,1] <- t_old
    } else{
    
      ### Compute corresponding values for y_0 and y_1
      y_0_n = sum(data[date <t_new]$event) -1
      y_1_n = sum(data[date>= t_new]$event) -1
      
      ### Compute the acceptance probability in log scale and
      #compare with u unif(0,1)
      a_log <- (((y_0_n +1)*log(lambda_0)
                 + (y_1_n +1)*log(lambda_1)
                 - t_new *(lambda_0 - lambda_1 )) -
                  ((y_0_o +1)*log(lambda_0)
                   + (y_1_o +1)*log(lambda_1)
                   - t_old *(lambda_0 - lambda_1 )))
      
      
      u <- runif(1)
      results[i, 5] <-  exp(a_log)
      if (u< exp(a_log)){
        results[i, 1] <-  t_new
      } else{
        results[i,1] <- t_old
      }
    }
  }
  return (results)
}

# Load libraries and data
library(boot)
library(data.table)
library(tidyr)
library(ggplot2)
library(latex2exp)

# Create data table
data <- as.data.table(coal)
data[, event:= 1]
data[, cum_event := cumsum(event)-1]
data[, date_shifted := date- data[1, date]]


q <- ggplot(data = data[1:(nrow(data)-1),], aes(x = date, y = cum_event))
q <- q + geom_line() + ggtitle("Cumulative plot of no. of accidents between 1851 and 1962") + 
  xlab("Year") + ylab("Number of accidents with 10 or more casualties")
q

ggsave(
  "cumulative_plot_data.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

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

t_2

# Run algorithm
sim_MH_single <- MH_single_alg(n,data, t_0, t_2, t, lambda_0,lambda_1, beta, sigma)

#Format and plot the results
sim_MH_single <- as.data.table(sim_MH_single)
setnames(sim_MH_single, c("t", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_single[, iteration := seq(1:n)]

q <- ggplot(data = sim_MH_single, aes(x = iteration) )
q <- q + geom_line(aes(y = lambda_0, colour = "lambda_0")) 
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for $\\lambda_0$"))))
q

ggsave(
  "sim_lambda0.pdf",
<<<<<<< HEAD
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
=======
  path = "C:\\Users\\sara_\\OneDrive\\Documents\\NTNU\\10.Semester\\Beregningskrevende\\Prosjekt1\\Berregningskrevende_2\\Images",
>>>>>>> f59f2de7abbf5f76d1e12d14c265b18b6518d5ee
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single, aes(x = iteration) )
q <- q + geom_line(aes(y = lambda_1, colour = "lambda_1"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_1$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for $\\lambda_1$"))))
q

ggsave(
  "sim_lambda1.pdf",
<<<<<<< HEAD
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
=======
  path = "C:\\Users\\sara_\\OneDrive\\Documents\\NTNU\\10.Semester\\Beregningskrevende\\Prosjekt1\\Berregningskrevende_2\\Images",
>>>>>>> f59f2de7abbf5f76d1e12d14c265b18b6518d5ee
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single, aes(x = iteration) )
q <- q + geom_line(aes(y = beta, colour = "beta"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for $\\beta$"))))
q

ggsave(
  "sim_beta.pdf",
<<<<<<< HEAD
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
=======
  path = "C:\\Users\\sara_\\OneDrive\\Documents\\NTNU\\10.Semester\\Beregningskrevende\\Prosjekt1\\Berregningskrevende_2\\Images",
>>>>>>> f59f2de7abbf5f76d1e12d14c265b18b6518d5ee
  width = 17,
  height = 10,
  units = "cm"
)
<<<<<<< HEAD
=======

>>>>>>> f59f2de7abbf5f76d1e12d14c265b18b6518d5ee
q <- ggplot(data = sim_MH_single, aes(x = iteration) )
q <- q + geom_line(aes(y = t, colour = "t"))
q <- q +  theme(legend.position = "none")
q <- q + ylab("t")
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for t"))))
q

ggsave(
  "sim_t.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


#mixing of the MH-algorithm
t_alt_1 = data[1, date] + 70
t_alt_2 = data[1, date] + 30

sim_MH_single_alt_1 <- MH_single_alg(n,data, t_0, t_2, t_alt_1, lambda_0,lambda_1, beta, sigma)
sim_MH_single_alt_2 <- MH_single_alg(n,data, t_0, t_2, t_alt_2, lambda_0,lambda_1, beta, sigma)

sim_MH_single_alt_1 <- as.data.table(sim_MH_single_alt_1)
sim_MH_single_alt_2 <- as.data.table(sim_MH_single_alt_2)
setnames(sim_MH_single_alt_1, c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_single_alt_2, c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_single_alt_1[, iteration := seq(1:n)]
sim_MH_single_alt_2[, iteration := seq(1:n)]

sim_MH_single[, run := 1]
sim_MH_single_alt_1[, run := 2]
sim_MH_single_alt_2[, run := 3]
sim_MH_single_compare <- rbindlist(list(sim_MH_single, sim_MH_single_alt_1, sim_MH_single_alt_2))

#[iteration <2600 & iteration >2400 ]
q <- ggplot(data = as.data.table(sim_MH_single_compare), aes(x = iteration, group = as.factor(run), colour = as.factor(run)))
q <- q + geom_line(aes(y = t))
q <- q +  theme(legend.position = "none")
q <- q + ylab("t")
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for t"))))
q

ggsave(
  "mixing_t.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = as.data.table(sim_MH_single_compare), aes(x = iteration, group = as.factor(run), colour = as.factor(run)))
q <- q + geom_line(aes(y = lambda_0))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for $\\lambda_0$"))))
q

ggsave(
  "mixing_lambda_0.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = as.data.table(sim_MH_single_compare), aes(x = iteration, group = as.factor(run), colour = as.factor(run)))
q <- q + geom_line(aes(y = lambda_0))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for $\\lambda_1$"))))
q

ggsave(
  "mixing_lambda_1.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)
q <- ggplot(data = as.data.table(sim_MH_single_compare), aes(x = iteration, group = as.factor(run), colour = as.factor(run)))
q <- q + geom_line(aes(y = beta))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for $\\beta$"))))
q

ggsave(
  "mixing_lambda_1.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

# Check if values seem reasonable
summary(sim_MH_single)
3.395*(1891- 1851.632) + 0.797*(1960.489-1891)

q <- ggplot(data = data, aes(x = date, y = cum_event))
q <- q + geom_line()
q

##Looking at the acf of t

burnin <- 500



sim_t <- sim_MH_single$t[burnin:n]

acf(sim_t)

no_burnin_sim_MH <- sim_MH_single[burnin:n]

summary(no_burnin_sim_MH)
