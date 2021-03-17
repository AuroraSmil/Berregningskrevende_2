#MH block

library(boot)
library(data.table)
library(tidyr)
library(ggplot2)
library(gridExtra)

# data <- as.data.table(coal)
# data[, event:= 1]
# data[, cum_event := cumsum(event)]
# data[, date := date- 1851]

data <- as.data.table(coal)
data[, event := 1]
data[, cum_event := cumsum(event) - 1]
data[, date_shifted := date - data[1, date]]

f_target_1 <-
  function(lambda_0,
           lambda_1,
           t,
           t_0,
           t_2,
           beta,
           y_0,
           y_1) {
    if (t < t_0) {
      return (0)
    }
    if (t > t_2) {
      return (0)
    }
    return(exp(
      -lambda_0 * (t - t_0) - lambda_1 * (t_2 - t) - (1 / beta) * (lambda_0 + lambda_1)
    ) * lambda_0 ^ {
      y_0 + 1
    } * lambda_1 ^ {
      y_1 + 1
    })
  }

f_prop_1 <-
  function(lambda_0,
           lambda_1,
           t_given,
           t_cur,
           t_0,
           t_2,
           beta,
           y_0,
           y_1,
           sigma_t) {
    if (t < t_0) {
      retval <- 0
    }
    if (t > t_2) {
      retval <- 0
    }
    retval <-
      (
        dgamma(lambda_0, (y_0 + 2), scale = (1 / (t_cur - t_0 + 1 / beta))) *
          dgamma(lambda_1, (y_1 + 2), scale = (1 / (t_2 - t_cur + 1 / beta))) *
          dnorm(t, t_given, sigma_t)
      )
    return(retval)
  }

f_target_2 <-
  function(lambda_0,
           lambda_1,
           beta,
           t_0,
           t_2,
           t,
           y_0,
           y_1) {
    if (beta < 0) {
      return (0)
    }
    
    return(exp(
      -lambda_0 * (t - t_0) - lambda_1 * (t_2 - t) - (1 / beta) *
        (lambda_0 + lambda_1 + 1)
    ) * lambda_0 ^ {
      y_0 + 1
    } *
      lambda_1 ^ {
        y_1 + 1
      } * (1 / beta ^ 5))
  }

f_prop_2 <-
  function(lambda_0,
           lambda_1,
           t,
           t_0,
           t_2,
           beta_cur,
           beta_given,
           y_0,
           y_1,
           sigma_beta) {
    retval <-
      dgamma(lambda_0, (y_0 + 2), scale = (1 / (t - t_0 + 1 / beta_cur))) *
      dgamma(lambda_1, (y_1 + 2), scale = (1 / (t_2 - t + 1 / beta_cur))) *
      dnorm(beta_cur, beta_given, sigma_beta)
    
    return (retval)
  }

MH_alg <-
  function(n,
           data,
           t_0,
           t_2,
           t,
           lambda_0,
           lambda_1,
           beta,
           sigma_t,
           sigma_beta) {
    y_0 = sum(data[date < t]$event) - 1
    y_1 = sum(data[date >= t]$event) - 1
    
    results <- matrix(nrow = n, ncol = 5)
    results[1, 1:5] <- c(t, lambda_0, lambda_1, beta, 1)
    for (i in 2:(n)) {
      #print(i)
      if (i %% 2 == 0) {
        beta <- results[i - 1, 4]
        t_old <- results[i - 1, 1]
        t_new <- rnorm(1, t_old, sigma_t)
        if (t_new < t_0 | t_new > t_2) {
          results[i, 1:5] <- c(t_old, lambda_0_old, lambda_1_old, beta, 0)
        } else{
          y_0_old = sum(data[date < t_old]$event) - 1
          y_1_old = sum(data[date >= t_old]$event) - 1
          y_0_n = sum(data[date < t_new]$event) - 1
          y_1_n = sum(data[date >= t_new]$event) - 1
          
          #gibbs step
          lambda_0_old <- results[i - 1, 2]
          lambda_1_old <- results[i - 1, 3]
          lambda_0 <-
            rgamma(1, (y_0_n + 2), scale = (1 / (t_new - t_0 + 1 / beta)))
          lambda_1 <-
            rgamma(1, (y_1_n + 2), scale = (1 / (t_2 - t_new + 1 / beta)))
          
          target_ratio <-
            f_target_1(lambda_0, lambda_1, t_new, t_0, t_2, beta, y_0_n, y_1_n) /
            f_target_1(lambda_0_old,
                       lambda_1_old,
                       t_old,
                       t_0,
                       t_2,
                       beta,
                       y_0_old,
                       y_1_old)
          
          prop_ratio <-
            f_prop_1(
              lambda_0_old,
              lambda_1_old,
              t_cur = t_old,
              t_given = t_new,
              t_0,
              t_2,
              beta,
              y_0_old,
              y_1_old,
              sigma_t = sigma_t
            ) /
            f_prop_1(
              lambda_0,
              lambda_1,
              t_new,
              t_given = t_old,
              t_0,
              t_2,
              beta,
              y_0_n,
              y_1_n,
              sigma_t = sigma_t
            )
          
          a <- target_ratio * prop_ratio
          a <- min(1, a)
          
          u <- runif(1)
          
          if (u < a) {
            results[i, 1:5] <- c(t_new, lambda_0, lambda_1, beta, a)
          } else{
            results[i, 1:5] <- c(t_old, lambda_0_old, lambda_1_old, beta, a)
          }
        }
      }
      else{
        t <- results[i - 1, 1]
        t
        beta_old <- results[i - 1, 4]
        
        beta_new <-
          rgamma(1, 6, scale = 1 / (lambda_0 + lambda_1 + 1))
        
        
        y_0 = sum(data[date < results[i - 1, 1]]$event) - 1
        y_1 = sum(data[date >= results[i - 1, 1]]$event) - 1
        if (beta_new < 0) {
          results[i, 1:5] <- c(t, lambda_0_old, lambda_1_old, beta_old, 0)
        } else {
          lambda_0_n <-
            rgamma(1, (y_0 + 2), scale = (1 / (t_old - t_0 + 1 / beta_new)))
          lambda_1_n <-
            rgamma(1, (y_1 + 2), scale = (1 / (t_2 - t_old + 1 / beta_new)))
          lambda_0_n
          lambda_1_n
          
          
          lambda_0_old <- results[i - 1, 2]
          lambda_1_old <- results[i - 1, 3]
          
          target_ratio <-
            f_target_2(lambda_0, lambda_1, beta_new, t_0, t_2, t, y_0, y_1) /
            f_target_2(lambda_0_old,
                       lambda_1_old,
                       beta_old,
                       t_0,
                       t_2,
                       t,
                       y_0,
                       y_1)
          
          
          
          prop_ratio <-
            f_prop_2(
              lambda_0_old,
              lambda_1_old,
              beta_cur = beta_old,
              beta_given = beta_new,
              t = t,
              t_0 = t_0,
              t_2 = t_2,
              y_0,
              y_1,
              sigma_beta = sigma_beta
            ) /
            f_prop_2(
              lambda_0,
              lambda_1,
              beta_cur = beta_new,
              beta_given = beta_old,
              t_0 = t_0,
              t_2 = t_2,
              t = t,
              y_0,
              y_1,
              sigma_beta = sigma_beta
            )
          
          
          a <- min(1, target_ratio * prop_ratio)
          u <- runif(1)
          if (u < a) {
            results[i, 1:5] <- c(t, lambda_0, lambda_1, beta_new, a)
          } else{
            results[i, 1:5] <- c(t, lambda_0_old, lambda_1_old, beta_old, a)
          }
        }
      }
    }
    return (results)
  }

t = data[1, date] + 50

lambda_0 <- 10
lambda_1 <-  5
beta = 3

n = 10000
t_0 = data[1, date]
t_2 = data[nrow(data), date] #maybe should be 1963?
sigma_t = 10
sigma_beta = 1

######## Fist run ####################
sim_MH <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t,
    lambda_0,
    lambda_1,
    beta,
    sigma_t = sigma_t,
    sigma_beta = sigma_beta
  )

sim_MH <- as.data.table(sim_MH)

setnames(sim_MH, c("t", "lambda_0", "lambda_1", "beta", "a"))

sim_MH[, itteration := seq(1:n)]

summary(sim_MH)

q <- ggplot(data = sim_MH, aes(x = itteration))
q <- q + geom_line(aes(y = t, colour = "t"))
q <- q +  theme(legend.position = "none")
q <- q + ylab("t")
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for t"))))
q

ggsave(
  "block_sim_t_sigma03.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


q <- ggplot(data = sim_MH, aes(x = itteration))
q <- q + geom_line(aes(y = lambda_0, colour = "lambda_0"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_0$"
))))
q

ggsave(
  "block_sim_lambda0.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


q <- ggplot(data = sim_MH, aes(x = itteration))
q <- q + geom_line(aes(y = lambda_1, colour = "lambda_1"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\lambda_1$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_1$"
))))
q

ggsave(
  "block_sim_lambda1.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


q <- ggplot(data = sim_MH, aes(x = itteration))
q <- q + geom_line(aes(y = beta, colour = "beta"))
q <- q +  theme(legend.position = "none")
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\beta$"
))))
q

ggsave(
  "block_sim_beta_sigma3.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


######## Eval burn in ################

t_1_1 = data[1, date] + 50
t_1_2 = data[1, date] + 10
t_1_3 = data[1, date] + 80
lambda_0 <- 10
lambda_1 <- 5
beta_1 = 3
beta_2 = 10
beta_3 = 0.1
n = 10000
sigma = 10
sim_MH_1 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_1,
    lambda_0,
    lambda_1,
    beta_1,
    sigma_t = sigma_t,
    sigma_beta = sigma_beta
  )
sim_MH_2 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_2,
    lambda_0,
    lambda_1,
    beta_2,
    sigma_t = sigma_t,
    sigma_beta = sigma_beta
  )
sim_MH_3 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_3,
    lambda_0,
    lambda_1,
    beta_3,
    sigma_t = sigma_t,
    sigma_beta = sigma_beta
  )

sim_MH_1 <- as.data.table(sim_MH_1)
sim_MH_2 <- as.data.table(sim_MH_2)
sim_MH_3 <- as.data.table(sim_MH_3)
setnames(sim_MH_1,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_2,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_3,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_1[, iteration := seq(1:n)]
sim_MH_2[, iteration := seq(1:n)]
sim_MH_3[, iteration := seq(1:n)]

sim_MH_1[, run := 1]
sim_MH_2[, run := 2]
sim_MH_3[, run := 3]
sim_MH_block_compare <-
  rbindlist(list(sim_MH_1, sim_MH_2, sim_MH_3))

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(run),
      colour = as.factor(run)
    )
  )
q <- q + geom_line(aes(y = t))
q <- q +  theme(legend.position = "none")
q <- q + ylab("t")
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for t"))))
q

ggsave(
  "mixing_t_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
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
  "mixing_lambda_0_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
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
  "mixing_lambda_1_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)
q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
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
  "mixing_beta_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


########### TUNING  sigma_t ######################

t_1_3 = data[1, date] + 80
lambda_0 <- 10
lambda_1 <- 5
beta_2 = 10

n = 10000
sigma_t_1 = 3
sigma_t_2 = 10
sigma_t_3 = 30
sim_MH_1 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_3,
    lambda_0,
    lambda_1,
    beta_2,
    sigma_t = sigma_t_1,
    sigma_beta = sigma_beta
  )
sim_MH_2 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_3,
    lambda_0,
    lambda_1,
    beta_2,
    sigma_t = sigma_t_2,
    sigma_beta = sigma_beta
  )
sim_MH_3 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_3,
    lambda_0,
    lambda_1,
    beta_2,
    sigma_t = sigma_t_2,
    sigma_beta = sigma_beta
  )

sim_MH_1 <- as.data.table(sim_MH_1)
sim_MH_2 <- as.data.table(sim_MH_2)
sim_MH_3 <- as.data.table(sim_MH_3)
setnames(sim_MH_1,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_2,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_3,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_1[, iteration := seq(1:n)]
sim_MH_2[, iteration := seq(1:n)]
sim_MH_3[, iteration := seq(1:n)]

sim_MH_1[, sigma_t := sigma_t_1]
sim_MH_2[, sigma_t := sigma_t_2]
sim_MH_3[, sigma_t := sigma_t_3]
sim_MH_block_compare <-
  rbindlist(list(sim_MH_1, sim_MH_2, sim_MH_3))

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_t),
      colour = as.factor(sigma_t)
    )
  )
q <- q + geom_line(aes(y = t))
q <- q + ylab("t")
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for t"))))
q

ggsave(
  "tuning_t_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_t),
      colour = as.factor(sigma_t)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_0$"
))))
q

ggsave(
  "tuning_lambda_0_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_t),
      colour = as.factor(sigma_t)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_1$"
))))
q

ggsave(
  "tuning_lambda_1_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_t),
      colour = as.factor(sigma_t)
    )
  )
q <- q + geom_line(aes(y = beta))
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\beta$"
))))
q

ggsave(
  "tuning_beta_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)


########### TUNING  beta ######################

t_1_3 = data[1, date] + 80
lambda_0 <- 10
lambda_1 <- 5
beta_2 = 10

n = 10000
sigma_beta_1 = 0.3
sigma_beta_2 = 1
sigma_beta_3 = 10
sim_MH_1 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_3,
    lambda_0,
    lambda_1,
    beta_2,
    sigma_t = sigma_t_2,
    sigma_beta = sigma_beta_1
  )
sim_MH_2 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_3,
    lambda_0,
    lambda_1,
    beta_2,
    sigma_t = sigma_t_2,
    sigma_beta = sigma_beta_1
  )
sim_MH_3 <-
  MH_alg(
    n,
    data,
    t_0,
    t_2,
    t_1_3,
    lambda_0,
    lambda_1,
    beta_2,
    sigma_t = sigma_t_2,
    sigma_beta = sigma_beta_1
  )

sim_MH_1 <- as.data.table(sim_MH_1)
sim_MH_2 <- as.data.table(sim_MH_2)
sim_MH_3 <- as.data.table(sim_MH_3)
setnames(sim_MH_1,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_2,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
setnames(sim_MH_3,
         c("t", "lambda_0", "lambda_1", "beta", "a"))
sim_MH_1[, iteration := seq(1:n)]
sim_MH_2[, iteration := seq(1:n)]
sim_MH_3[, iteration := seq(1:n)]

sim_MH_1[, sigma_beta := sigma_beta_1]
sim_MH_2[, sigma_beta := sigma_beta_2]
sim_MH_3[, sigma_beta := sigma_beta_3]
sim_MH_block_compare <-
  rbindlist(list(sim_MH_1, sim_MH_2, sim_MH_3))

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_beta),
      colour = as.factor(sigma_beta)
    )
  )
q <- q + geom_line(aes(y = t))
q <- q + ylab("t")
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c("Traceplot for t"))))
q

ggsave(
  "tuning_t_block_beta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_beta),
      colour = as.factor(sigma_beta)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_0$"
))))
q

ggsave(
  "tuning_lambda_0_block_beta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_beta),
      colour = as.factor(sigma_beta)
    )
  )
q <- q + geom_line(aes(y = lambda_0))
q <- q + ylab(unname(TeX(c("$\\lambda_0$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\lambda_1$"
))))
q

ggsave(
  "tuning_lambda_1_block_beta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <-
  ggplot(
    data = as.data.table(sim_MH_block_compare),
    aes(
      x = iteration,
      group = as.factor(sigma_beta),
      colour = as.factor(sigma_beta)
    )
  )
q <- q + geom_line(aes(y = beta))
q <- q + ylab(unname(TeX(c("$\\beta$"))))
q <- q + xlab("Iteration")
q <- q + ggtitle(unname(TeX(c(
  "Traceplot for $\\beta$"
))))
q

ggsave(
  "tuning_beta_block_beta.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)





#grid.arrange(q1,q2,q3,q4, ncol = 2)

#only beta has a burn in period! kind of

summary(sim_MH)
3 * 41 + 0.92 * (112 - 41)

burnin <- 1000

# sim_t <- sim_MH$t[burnin:n]
# acf(sim_t)

no_burnin_sim_MH <- sim_MH[burnin:n]

summary(no_burnin_sim_MH)




######### Posterior ####################

mean <- summary(no_burnin_sim_MH)[4, ]
cov <- cov(no_burnin_sim_MH[, lambda_1], no_burnin_sim_MH[, lambda_0])

mean
cov

q <- ggplot(data = no_burnin_sim_MH, aes(x = t))
q <- q + geom_histogram(aes(y = ..density..), bins = 150)
q <- q + xlab(unname(TeX(c("$t$"))))
q <- q + ylab("Density")
q <- q + ggtitle(unname(TeX("Posterior density of $t$")))
q


ggsave(
  "post_t_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = no_burnin_sim_MH, aes(x = beta))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + xlab(unname(TeX(c("$\\beta$"))))
q <- q + ylab("Density")
q <- q + ggtitle(unname(TeX("Posterior density of $\\beta$")))
q

ggsave(
  "post_beta_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = no_burnin_sim_MH, aes(x = lambda_0))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + xlab(unname(TeX(c("$\\lambda_0$"))))
q <- q + ylab("Density")
q <- q + ggtitle(unname(TeX("Posterior density of $\\lambda_0$")))
q
ggsave(
  "post_lambda_0_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)

q <- ggplot(data = no_burnin_sim_MH, aes(x = lambda_0))
q <- q + geom_histogram(aes(y = ..density..), bins = 100)
q <- q + xlab(unname(TeX(c("$\\lambda_1$"))))
q <- q + ylab("Density")
q <- q + ggtitle(unname(TeX("Posterior density of $\\lambda_1$")))
q

ggsave(
  "post_lambda_1_block.pdf",
  path = "/Users/aurorahofman/Documents/NTNU/5 klasse/Beregningskrevende statistikk/Berregningskrevende_2/Images",
  width = 17,
  height = 10,
  units = "cm"
)
