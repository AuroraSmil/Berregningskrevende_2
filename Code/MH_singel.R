

library(ggplot2)
library(data.table)
library(Matrix)
library(MASS)
library(mvtnorm)
library(latex2exp)
library(tikzDevice)

## SINGLE SITE MCMC

MH_single_alg <-
  function(n,
           data,
           t_0,
           t_2,
           t,
           lambda_0,
           lambda_1,
           beta,
           sigma) {
    # Find number of deaths in each interval given initial t substract one to
    # adjust for h data stat and end date{}
    y_0 = sum(data[date < t]$event) - 1
    y_1 = sum(data[date >= t]$event) - 1
    
    # Construct result matrix
    results <- matrix(nrow = n, ncol = 5)
    results[1, 1:5] <- c(t, lambda_0, lambda_1, beta, 1)
    for (i in 2:n) {
      #print(i)
      # Collect old values for t and beta and compute corresponding
      # y_0 and y_1 values
      t_old <- results[i - 1, 1]
      
      y_0_o = sum(data[date < t_old]$event) - 1
      y_1_o = sum(data[date >= t_old]$event) - 1
      
      beta <- results[i - 1, 4]
      
      # Preform Gibbs step
      lambda_0 <-
        rgamma(1, y_0_o + 2, scale = 1 / (t_old - t_0 + 1 / beta))
      lambda_1 <-
        rgamma(1, y_1_o + 2, scale = 1 / (t_2 - t_old + 1 / beta))
      
      beta_tilde <- rgamma(1, 6, scale = 1 / (lambda_0 + lambda_1 + 1))
      beta <- 1 / beta_tilde
      
      results[i, 2:4] <- c(lambda_0, lambda_1, beta)
      
      # Preform MH step with the random walk of order one
      # as the proposal density
      t_new <- rnorm(1, t_old, sigma)
      if (t_new <= t_0) {
        results[i, 1] <- t_old
      } else if (t_new >= t_2) {
        results[i, 1] <- t_old
      } else{
        ### Compute corresponding values for y_0 and y_1
        y_0_n = sum(data[date < t_new]$event) - 1
        y_1_n = sum(data[date >= t_new]$event) - 1
        
        ### Compute the acceptance probability in log scale and
        #compare with u unif(0,1)
        a_log <- (((y_0_n + 1) * log(lambda_0)
                   + (y_1_n + 1) * log(lambda_1)
                   - t_new * (lambda_0 - lambda_1)
        ) -
          ((y_0_o + 1) * log(lambda_0)
           + (y_1_o + 1) * log(lambda_1)
           - t_old * (lambda_0 - lambda_1)
          ))
        
        
        u <- runif(1)
        results[i, 5] <-  exp(a_log)
        if (u < exp(a_log)) {
          results[i, 1] <-  t_new
        } else{
          results[i, 1] <- t_old
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

# Create data table ----
data <- as.data.table(coal)
data[, event := 1]
data[, cum_event := cumsum(event) - 1]
data[, date_shifted := date - data[1, date]]


q <-
  ggplot(data = data[1:(nrow(data) - 1), ], aes(x = date, y = cum_event))
q <-
  q + geom_line() + ggtitle("Cumulative plot of no. of accidents between 1851 and 1962") +
  xlab("Year") + ylab("Number of accidents with 10 or more casualties")
q

ggsave(
  "cumulative_plot_data.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

# Initial values for t, lambda_0, lambda_1, beta ----
t_1_1 = data[1, date] + 50
t_1_2 = data[1, date] + 10
t_1_3 = data[1, date] + 80
lambda_0 <- 10
lambda_1 <- 5
beta = 3


# Set n, t_0, t_2 and tuning parameter sigma
n = 100000
t_0 = data[1, date]
t_2 = data[nrow(data), date] #maybe should be 1963?
sigma = 10


# Run algorithm
sim_MH_single <-
  MH_single_alg(n, data, t_0, t_2, t_1_1, lambda_0, lambda_1, beta, sigma)

#Format and plot the results -----
sim_MH_single <- as.data.table(sim_MH_single)
setnames(sim_MH_single, c("t", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_single[, iteration := seq(1:n)]

q <- ggplot(data = sim_MH_single, aes(x = iteration))
q <- q + geom_line(aes(y = lambda_0, colour = "lambda_0"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_0$"
))))
q

ggsave(
  "sim_lambda0.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single, aes(x = iteration))
q <- q + geom_line(aes(y = lambda_1, colour = "lambda_1"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_1$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_1$"
))))
q

ggsave(
  "sim_lambda1.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single, aes(x = iteration))
q <- q + geom_line(aes(y = beta, colour = "beta"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\beta$"
))))
q

ggsave(
  "sim_beta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single, aes(x = iteration))
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


#mixing of the MH-algorithm ------

n = 10000
sigma = 10
sim_MH_single_1 <-
  MH_single_alg(n, data, t_0, t_2, t_1_1, lambda_0, lambda_1, beta, sigma)
sim_MH_single_2 <-
  MH_single_alg(n, data, t_0, t_2, t_1_2, lambda_0, lambda_1, beta, sigma)
sim_MH_single_3 <-
  MH_single_alg(n, data, t_0, t_2, t_1_3, lambda_0, lambda_1, beta, sigma)

sim_MH_single_1 <- as.data.table(sim_MH_single_1)
sim_MH_single_2 <- as.data.table(sim_MH_single_2)
sim_MH_single_3 <- as.data.table(sim_MH_single_3)
setnames(sim_MH_single_1,
         c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_single_2,
         c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_single_3,
         c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_single_1[, iteration := seq(1:n)]
sim_MH_single_2[, iteration := seq(1:n)]
sim_MH_single_3[, iteration := seq(1:n)]

sim_MH_single_1[, run := 1]
sim_MH_single_2[, run := 2]
sim_MH_single_3[, run := 3]
sim_MH_single_compare <-
  rbindlist(list(sim_MH_single_1, sim_MH_single_2, sim_MH_single_3))

#[iteration <2600 & iteration >2400 ]
q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = t_alt))
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

q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_0$"
))))
q

ggsave(
  "mixing_lambda_0.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_1$"
))))
q

ggsave(
  "mixing_lambda_1.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)
q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = beta))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\beta$"
))))
q

ggsave(
  "mixing_beta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

# Check if values seem reasonable
summary(sim_MH_single)
3.395 * (1891 - 1851.632) + 0.797 * (1960.489 - 1891)

q <- ggplot(data = data, aes(x = date, y = cum_event))
q <- q + geom_line()
q

##Looking at the acf of t

burnin <- 500

sim_t <- sim_MH_single$t[burnin:n]

acf(sim_t)

no_burnin_sim_MH <- sim_MH_single[burnin:n]

summary(no_burnin_sim_MH)




## TUNING #######

n = 100000
sigma_1 = 10
sigma_2 = 1.5
sigma_3  = 20
sim_MH_single_1 <-
  MH_single_alg(n, data, t_0, t_2, t_1_3, lambda_0, lambda_1, beta, sigma_1)
sim_MH_single_2 <-
  MH_single_alg(n, data, t_0, t_2, t_1_3, lambda_0, lambda_1, beta, sigma_2)
sim_MH_single_3 <-
  MH_single_alg(n, data, t_0, t_2, t_1_3, lambda_0, lambda_1, beta, sigma_3)

sim_MH_single_1 <- as.data.table(sim_MH_single_1)
sim_MH_single_2 <- as.data.table(sim_MH_single_2)
sim_MH_single_3 <- as.data.table(sim_MH_single_3)
setnames(sim_MH_single_1,
         c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_single_2,
         c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_single_3,
         c("t_alt", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_single_1[, iteration := seq(1:n)]
sim_MH_single_2[, iteration := seq(1:n)]
sim_MH_single_3[, iteration := seq(1:n)]

sim_MH_single_1[, run := 1]
sim_MH_single_2[, run := 2]
sim_MH_single_3[, run := 3]
sim_MH_single_compare <-
  rbindlist(list(sim_MH_single_1, sim_MH_single_2, sim_MH_single_3))

#[iteration <2600 & iteration >2400 ]
q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = t_alt))
#q <- q +  theme(legend.position = "none")
q <- q + ylab("t")
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for t"))))
q

ggsave(
  "tuning_t.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_0$"
))))
q

ggsave(
  "tuning_lambda_0.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_1$"
))))
q

ggsave(
  "mixing_lambda_1.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)
q <-
  ggplot(
    data = as.data.table(sim_MH_single_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = beta))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\beta$"
))))
q

ggsave(
  "tuning_beta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)



##################### Density plots ##################################################
sim_MH_single
mean<- summary(sim_MH_single)[4,]
cov <- cov(sim_MH_single[2000:nrow(sim_MH_single),lambda_1], sim_MH_single[2000:nrow(sim_MH_single),lambda_0])

sum[4,]

q <- ggplot(data = sim_MH_single[2000:nrow(sim_MH_single),], aes(x = t))
q <- q + geom_histogram(aes(y = ..density..), bins = 150) 
q <- q + scale_x_continuous( limits = c(1880,1900))
q <- q + scale_y_continuous( limits = c(0,0.6))
q <- q + xlab(unname(TeX(c("$t$"))))
q <- q + ylab("Density")
q <- q + ggtitle(unname(TeX("Posterior density of $t$")))
q


ggsave(
  "post_t_single.pdf",
  path = "C:\\Users\\sara_\\OneDrive\\Documents\\NTNU\\10.Semester\\Beregningskrevende\\Prosjekt1\\Berregningskrevende_2\\Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single[2000:nrow(sim_MH_single),], aes(x = beta))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + xlab(unname(TeX(c("$\\beta$"))))
q <- q + scale_x_continuous( limits = c(0,4))
q <- q + scale_y_continuous( limits = c(0,1.5))
q <- q + ylab("Density")
q <- q + ggtitle(unname(TeX("Posterior density of $\\beta$")))
q

ggsave(
  "post_beta_single.pdf",
  path = "C:\\Users\\sara_\\OneDrive\\Documents\\NTNU\\10.Semester\\Beregningskrevende\\Prosjekt1\\Berregningskrevende_2\\Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single[2000:nrow(sim_MH_single),], aes(x = lambda_0))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + xlab(unname(TeX(c("$\\lambda_0$"))))
q <- q + ylab("Density")
q <- q + scale_x_continuous( limits = c(1,5))
#q <- q + scale_y_continuous( limits = c(0,1.3))
q <- q + ggtitle(unname(TeX("Posterior density of $\\lambda_0$")))
q
ggsave(
  "post_lambda_0_single.pdf",
  path = "C:\\Users\\sara_\\OneDrive\\Documents\\NTNU\\10.Semester\\Beregningskrevende\\Prosjekt1\\Berregningskrevende_2\\Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = sim_MH_single[2000:nrow(sim_MH_single),], aes(x = lambda_1))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + xlab(unname(TeX(c("$\\lambda_1$"))))
q <- q + ylab("Density")
q <- q + scale_x_continuous( limits = c(0,2))
#q <- q + scale_y_continuous( limits = c(0,2))
q <- q + ggtitle(unname(TeX("Posterior density of $\\lambda_1$")))
q

ggsave(
  "post_lambda_1_single.pdf",
  path = "C:\\Users\\sara_\\OneDrive\\Documents\\NTNU\\10.Semester\\Beregningskrevende\\Prosjekt1\\Berregningskrevende_2\\Images",
  width = 17,
  height = 10,
  units = "cm"
)

