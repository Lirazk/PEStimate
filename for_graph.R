# Not sure why I used ggplot for the first one
# and tinyplot for the rest.
library(ggplot2)
library(tinyplot)
library(tidyr)
library(dplyr)

source("R/EmbryoSelection.R")

# Figure 1
# Random number of embryos plot
r2 <- 0.1
K <- 0.01
N <- 1:10
p <- 0.4
# For binomial, we solve:
# n_0 * p /(1-(1-p)^n_0) = N
n0 <- sapply(N, function(N) 
  optimise(function(n0) (n0*p/(1-(1-p)^n0) - N)^2, 
           c(0.001, 100))$minimum)
# For poisson, we solve
# lambda * p /(1-exp(-lambda*p)) = N
lambda <- sapply(N, function(N) 
  optimise(function(lambda) (lambda*p/(1-exp(-lambda*p)) - N)^2, 
         c(0.001, 100))$minimum)
result <- matrix(nrow=length(N), ncol=3)
for (n in N) {
  result[n, 1] <- 100*risk_reduction_lowest(r2, K, n)
  result[n, 2] <- 100*risk_reduction_lowest_bin(r2, K, n0[n], p)
  result[n, 3] <- 100*risk_reduction_lowest_pois(r2, K, lambda[n] * p)
}

result <- data.frame(result) 
names(result) <- c("Fixed live births", "Binomial", "Poisson")
result$n <- N
result_long <- pivot_longer(result, cols = 1:3)
result_long <- result_long |>
  subset(n > 1)

result_long$name <- factor(result_long$name, 
                           levels = c("Fixed live births", "Binomial", "Poisson"))

ggplot(data=result_long, aes(n, value, col=name)) +
  geom_point(size=3) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = "bottom") +
  labs(x = "Expected number of live births",
       y = "Relative risk reduction(%)") +
  ylim(c(0, max(result_long$value)))

ggsave("~/Desktop/figure.png",
       width=8, height=5)

# Figure S1
# Shared environment
shared_env <- function(iter, h2, r2, K, n, pi, case=1) {
  # Cases -
  # 1: Unconditional
  # 2: Condition on sib
  # 3: Condition on two sibs
  # 4: Condition on a parent
  # 5: Condition on a grandparent
  stopifnot(case %in% 1:5)
  
  # Shared genetic components.
  c <- rnorm(iter, sd=sqrt(r2/2))
  v <- rnorm(iter, sd=sqrt((h2-r2)/2))
  e_shared <- rnorm(iter, sd=sqrt(pi*(1-h2)))
  
  zk <- qnorm(1-K)
  index <- 1:iter
  
  cond_list <- list()
  
  if (case == 2) {
    cond_list$sib_self <- 1
    
    ys <- c+v+e_shared+rnorm(iter, sd=sqrt(1-h2/2-pi*(1-h2)))
    index <- which(ys > zk)
  }
  if (case == 3) {
    cond_list$sib_self <- c(1, 1)
    
    ys1 <- c+v+e_shared+rnorm(iter, sd=sqrt(1-h2/2-pi*(1-h2)))
    ys2 <- c+v+e_shared+rnorm(iter, sd=sqrt(1-h2/2-pi*(1-h2)))
    index <- which(ys1 > zk & ys2 > zk)
  }
  if (case == 4) {
    cond_list$p1 <- 1
    
    sf <- rnorm(iter, sd=sqrt(r2))
    sm <- rnorm(iter, sd=sqrt(r2))
    gf <- rnorm(iter, sd=sqrt(h2-r2))
    gm <- rnorm(iter, sd=sqrt(h2-r2))
    
    c <- (sf+sm)/2
    v <- (gf+gm)/2
    
    ys <- sf+gf+e_shared+rnorm(iter, sd=sqrt(1-h2-pi*(1-h2)))
    index <- which(ys > zk)
  }
  if (case == 5) {
    cond_list$gp1a <- 1
    # How pretty.
    gp1sf <- rnorm(iter, sd=sqrt(r2))
    gp1sm <- rnorm(iter, sd=sqrt(r2))
    gp1gf <- rnorm(iter, sd=sqrt(h2-r2))
    gp1gm <- rnorm(iter, sd=sqrt(h2-r2))
    sf <- (gp1sf + gp1sm)/2+rnorm(iter, sd=sqrt(r2/2))
    gf <- (gp1gf + gp1gm)/2+rnorm(iter, sd=sqrt((h2-r2)/2))
    
    gp2sf <- rnorm(iter, sd=sqrt(r2))
    gp2sm <- rnorm(iter, sd=sqrt(r2))
    gp2gf <- rnorm(iter, sd=sqrt(h2-r2))
    gp2gm <- rnorm(iter, sd=sqrt(h2-r2))
    sm <- (gp2sf + gp2sm)/2+rnorm(iter, sd=sqrt(r2/2))
    gm <- (gp2gf + gp2gm)/2+rnorm(iter, sd=sqrt((h2-r2)/2))
    
    c <- (sf+sm)/2
    v <- (gf+gm)/2
    
    ys <- gp1sf+gp1gf+e_shared+rnorm(iter, sd=sqrt(1-h2-pi*(1-h2)))
    index <- which(ys > zk)
  }
  
  adjust_h2 <- case_when(
    case %in% 1:4 ~ (h2*(1-2*pi) + 2*pi),
    case == 5 ~ (h2+4*pi*(1-h2))
  )
  
  # y <- replicate(iter, min(rnorm(n, sd=sqrt(r2/2))))+ 
  #   c +
  #   v +
  #   e_shared +
  #   # Specific genetic component
  #   rnorm(iter, sd=sqrt((h2-r2)/2)) +
  #   # Specific environment
  #   rnorm(iter, sd=sqrt((1-pi)*(1-h2)))
  # 
  # baseline <- mean(c[index]+
  #                    v[index]+
  #                    e_shared[index]+
  #                    rnorm(length(index), sd=sqrt(r2/2)) +
  #                    rnorm(length(index), sd=sqrt((h2-r2)/2)) +
  #                    rnorm(length(index), sd=sqrt((1-pi)*(1-h2))) > zk)
  # 
  # risk <- mean(y[index] > zk)
  y <- replicate(iter, min(rnorm(n, sd=sqrt(r2/2))))+ 
    c +
    v +
    e_shared
  
  baseline <- mean(pnorm(c[index]+
                     v[index]+
                     e_shared[index] - zk,
                   sd = sqrt(h2/2 + (1-pi)*(1-h2))))
  
  risk <- mean(pnorm(y[index] - zk, 
                     sd = sqrt((h2-r2)/2 + (1-pi)*(1-h2))))
  
  adjusted <- risk_prediction_exact(iter, n, r2, 
                                    adjust_h2, 
                                    K, cond_list)
  
  c(baseline=baseline, selection=risk, 
    relative_reduction = (baseline-risk)/baseline,
    adjusted[1], adjusted[2], 
    adjusted[3])
}

plot_shared_env <- function(results, pi) {
  par(mfrow = c(1, 2))
  plt(results[3, ] ~ pi, 
      xlab = expression(pi), ylab = "rrr",
      ylim = c(min(results[3, ], results[6, ]),
               max(results[3, ], results[6, ])), pch=19, 
      theme="classic")
  points(results[6, ] ~ pi, col = "red", pch=19)
  
  arr1 <- results[1, ] - results[2, ]
  arr2 <- results[4, ] - results[5, ]
  plt(arr1 ~ pi, xlab = expression(pi), ylab = "arr",
      ylim = c(min(arr1, arr2),
               max(arr1, arr2)), pch=19,
      theme="classic")
  points(arr2 ~ pi, col = "red", pch=19)
}
# Parameters. Somewhat large K so the conditional case
# variance would be lower. Although there are many
# other ways to reduce the variance.
iter <- 1000000
r2 <- 0.1
h2 <- 0.3
K <- 0.1
n <- 5
pi <- seq(0, 0.5, 0.05)

# Unconditional
results <- sapply(pi, function(pi) 
  shared_env(iter, h2, r2, K, n, pi, 1))
plot_shared_env(results, pi)

# One sibling
results <- sapply(pi, function(pi) 
  shared_env(iter, h2, r2, K, n, pi, 2))
plot_shared_env(results, pi)

# Two siblings
results <- sapply(pi, function(pi) 
  shared_env(iter, h2, r2, K, n, pi, 3))
plot_shared_env(results, pi)

# Parent
results <- sapply(pi, function(pi) 
  shared_env(iter, h2, r2, K, n, pi, 4))
plot_shared_env(results, pi)

# Grandparent
pi <- seq(0, 0.25, 0.05)
results <- sapply(pi, function(pi) 
  shared_env(iter, h2, r2, K, n, pi, 5))
plot_shared_env(results, pi)

# Figure S2 - not used
# Correlation between g and s
correlated_genetic <- function(iter, var_env, r2, K, n, rho,
                               case=1) {
  # Cases -
  # 1: Unconditional
  # 2: Condition on a sibling
  # 3: Condition on two siblings
  # 4: Condition on a parent
  # 5: Condition on a parent+sibling
  stopifnot(case %in% 1:5)
  
  zk <- qnorm(1-K)
  # Find the h2 which would result in a specific
  # environment variance.
  # Can also be solved explicitly.
  h2 <- optimise(function(h2) abs(h2+2*rho * 
                                    sqrt(r2 * (h2-r2))-
                                    (1-var_env)), 
                 c(r2+1e-8, 1))$minimum
  
  cov <- rho * sqrt(r2 * (h2-r2))
  # For each parent c and v
  joint_cv1 <- mvnfast::rmvn(iter, mu = c(0, 0), 
                             sigma = rbind(c(r2, cov),
                                           c(cov, h2-r2)))
  joint_cv2 <- mvnfast::rmvn(iter, mu = c(0, 0), 
                             sigma = rbind(c(r2, cov),
                                           c(cov, h2-r2)))
  
  c <- (joint_cv1[, 1] + joint_cv2[, 1])/2
  v <- (joint_cv1[, 2] + joint_cv2[, 2])/2
  
  cond_list <- list()
  index <- 1:iter
  if (case == 2) {
    cond_list <- list(sib_self = 1)
    sib <- mvnfast::rmvn(iter, mu = c(0, 0), 
                         sigma = rbind(c(r2/2, cov/2),
                                       c(cov/2, (h2-r2)/2)))
    ys <- sib[, 1] + sib[, 2] + c + v + 
      rnorm(iter, sd=sqrt(var_env))
    index <- which(ys > zk)
  }
  if (case == 3) {
    cond_list <- list(sib_self = c(1, 1))
    sib <- mvnfast::rmvn(iter, mu = c(0, 0), 
                         sigma = rbind(c(r2/2, cov/2),
                                       c(cov/2, (h2-r2)/2)))
    ys1 <- sib[, 1] + sib[, 2] + c + v + 
      rnorm(iter, sd=sqrt(var_env))
    
    sib <- mvnfast::rmvn(iter, mu = c(0, 0), 
                         sigma = rbind(c(r2/2, cov/2),
                                       c(cov/2, (h2-r2)/2)))
    ys2 <- sib[, 1] + sib[, 2] + c + v + 
      rnorm(iter, sd=sqrt(var_env))
    index <- which(ys1 > zk & ys2 > zk)
  }
  if (case == 4) {
    cond_list <- list(p1 = 1)
    
    ys <- joint_cv1[, 1] + joint_cv1[, 2] +
      rnorm(iter, sd=sqrt(1-h2-2*cov))
    index <- which(ys > zk)
  }
  if (case == 5) {
    cond_list <- list(sib_self = 1, p1 = 1)
    sib <- mvnfast::rmvn(iter, mu = c(0, 0), 
                         sigma = rbind(c(r2/2, cov/2),
                                       c(cov/2, (h2-r2)/2)))
    ys1 <- sib[, 1] + sib[, 2] + c + v + 
      rnorm(iter, sd=sqrt(var_env))
    
    ys2 <- joint_cv1[, 1] + joint_cv1[, 2] +
      rnorm(iter, sd=sqrt(1-h2-2*cov))
    index <- which(ys1 > zk & ys2 > zk)
  }
  
  # Next sample the specific x_i, v_i
  temp <- replicate(n, {
    mvnfast::rmvn(iter, mu = c(0, 0), 
                  sigma = rbind(c(r2/2, cov/2),
                                c(cov/2, (h2-r2)/2)))
  })
  
  # Find the minimum embryo and add the terms for y_min.
  idx <- apply(temp[, 1, ], 1, which.min)
  # y <- sapply(1:iter, function(i)
  #   temp[i, 1, idx[i]]+temp[i, 2, idx[i]]
  # ) + c + v + rnorm(iter, sd=sqrt(var_env))
  # 
  # y_baseline <- temp[, 1, 1]+temp[, 2, 1] + c + v + 
  #   rnorm(iter, sd=sqrt(var_env))
  # 
  # baseline <- mean(y_baseline[index] > zk)
  # risk <- mean(y[index] > zk)
  
  # For some variance reduction
  y <- sapply(1:iter, function(i)
    temp[i, 1, idx[i]]+temp[i, 2, idx[i]]
  ) + c + v
  
  y_baseline <- temp[, 1, 1]+temp[, 2, 1] + c + v #+ 

  baseline <- mean(pnorm(y_baseline[index]-zk, sd=sqrt(var_env)))
  risk <- mean(pnorm(y[index]-zk, sd=sqrt(var_env)))
  
  if (case == 1) {
    risk_original_model <- (risk_prediction_analytical(r2+2*rho*sqrt(r2*(h2-r2))+
                                                         rho^2*(h2-r2), K, n) |>
                              unlist())[1:2]
  }
  else {
    risk_original_model <- risk_prediction_exact(iter, n, 
                                                 r2+2*rho*sqrt(r2*(h2-r2))+rho^2*(h2-r2), 
                                                 h2+2*rho*sqrt(r2*(h2-r2)), 
                                                 K, cond_list)[1:2]
  }
  
  c(baseline_sim = baseline, 
    selection_sim = risk,
    risk_original_model)
}

plot_correlated_genetic <- function(result) {
  rho <- result$rho
  arr1 <- result$baseline-result$selection
  rrr1 <- arr1/result$baseline
  arr2 <- result$baseline_sim-result$selection_sim
  rrr2 <- arr2/result$baseline_sim
  
  par(mfrow = c(1, 2))
  plt(rrr1 ~ rho, xlab = expression(rho), 
      ylab = "rrr",
      ylim = c(min(rrr1, rrr2),
               max(rrr1, rrr2)), pch=19, 
      theme="classic")
  points(rrr2 ~ rho, col = "red", pch=19)
  
  plt(arr1 ~ rho, xlab = expression(rho), ylab = "arr",
      ylim = c(min(arr1, arr2),
               max(arr1, arr2)), pch=19,
      theme="classic")
  points(arr2 ~ rho, col = "red", pch=19)
  # par()
}

# Parameters
rho <- seq(0, 0.9, length=10)
n <- 2
K <- 0.1
var_env <- 0.7

# Unconditional
result <- sapply(rho, function(rho) 
  correlated_genetic(iter, var_env, r2, K, n, rho))
result <- data.frame(cbind(rho, t(result)))

plot_correlated_genetic(result)

# One sibling
result <- sapply(rho, function(rho) 
  correlated_genetic(iter, var_env, r2, K, n, rho, 2))
result <- data.frame(cbind(rho, t(result)))

plot_correlated_genetic(result)

# Two siblings
result <- sapply(rho, function(rho) 
  correlated_genetic(iter, var_env, r2, K, n, rho, 3))
result <- data.frame(cbind(rho, t(result)))

plot_correlated_genetic(result)

# Parent
result <- sapply(rho, function(rho) 
  correlated_genetic(iter, var_env, r2, K, n, rho, 4))
result <- data.frame(cbind(rho, t(result)))

plot_correlated_genetic(result)

# Parent + siblings
result <- sapply(rho, function(rho) 
  correlated_genetic(iter, var_env, r2, K, n, rho, 5))
result <- data.frame(cbind(rho, t(result)))

plot_correlated_genetic(result)

# Figure S3
# Correlation between genetic component and environment

# env_gen <- function(iter, n, h2, r2, K, w, rho1, rho2, cond=F) {
env_gen <- function(iter, n, h2, r2, K, rho1, rho2, cond=F) {
  zk <- qnorm(1-K)
  # s, g, e. Ended up removing w, and just dividing rho1 and rho2
  # equally between score and non-score genetic part.
  # cor(s_1, s_2) = r2/2
  # cor(g_1, g_2) = (h2-r2)/2
  # cor(s_1, e_1) = w rho1
  # cor(g_1, e_1) = (1-w) rho1
  # cor(s_1, e_2) = w rho2
  # cor(g_1, e_2) = (1-w) rho2
  var_s <- r2
  var_g <- h2-r2
  
  # Again can be calculated explicitly, but that is simpler
  var_e <- optimize(function(x)
    abs(1-(h2+x+2*rho1*sqrt(r2*x)+2*rho1*sqrt((h2-r2)*x))), 
    c(0, 1))
  var_e <- var_e$minimum
  
  # Pretty ugly, but define the covariance matrix
  # for siblings with the correlated g, s, e
  sigma <- diag(c(rep(var_s, n+cond),
                  rep(var_g, n+cond),
                  rep(var_e, n+cond)))
  s_index <- 1:(n+cond)
  g_index <- (n+cond+1):(2*(n+cond))
  e_index <- (2*(n+cond)+1):(3*(n+cond))
  
  count <- 1
  for (i in 1:(n+cond)) {
    sigma[i, s_index] <- var_s/2
    sigma[i, i] <- var_s
    
    # sigma[i, e_index] <- w * rho2 * sqrt(var_s * var_e)
    # sigma[i, e_index[count]] <- w*rho1 * sqrt(var_s * var_e)
    sigma[i, e_index] <- rho2 * sqrt(var_s * var_e)
    sigma[i, e_index[count]] <- rho1 * sqrt(var_s * var_e)
    
    # sigma[e_index, i] <- w * rho2 * sqrt(var_s * var_e)
    # sigma[e_index[count], i] <- w*rho1 * sqrt(var_s * var_e)
    sigma[e_index, i] <- rho2 * sqrt(var_s * var_e)
    sigma[e_index[count], i] <- rho1 * sqrt(var_s * var_e)
    
    count <- count + 1
  }
  
  count <- 1
  for (i in (n+cond+1):(2*(n+cond))) {
    sigma[i, g_index] <- (h2-r2)/2
    sigma[i, i] <- h2-r2
    
    # sigma[i, e_index] <- (1-w) * rho2 * sqrt(var_g * var_e)
    # sigma[i, e_index[count]] <- (1-w)*rho1 * sqrt(var_g * var_e)
    sigma[i, e_index] <- rho2 * sqrt(var_g * var_e)
    sigma[i, e_index[count]] <- rho1 * sqrt(var_g * var_e)
    
    # sigma[e_index, i] <- (1-w) * rho2 * sqrt(var_g * var_e)
    # sigma[e_index[count], i] <- (1-w)*rho1 * sqrt(var_g * var_e)
    sigma[e_index, i] <- rho2 * sqrt(var_g * var_e)
    sigma[e_index[count], i] <- rho1 * sqrt(var_g * var_e)
    
    count <- count + 1
  }
  
  temp <- mvnfast::rmvn(iter, numeric(3*(n+cond)), sigma)
  s <- temp[, s_index]
  g <- temp[, g_index]
  e <- temp[, e_index]
  y <- s + g + e
  # mean(y[, 1] > zk)
  if (cond) {
    index <- which(y[, 1] > zk)
    
    idx <- apply(s[, -1], 1, which.min)
    y_min <- sapply(1:nrow(y), function(i) y[i, 1+idx[i]])
    baseline <- mean(y[index, 2] > zk)
    risk <- mean(y_min[index] > zk)
  }
  else {
    baseline <- mean(y > zk)
    idx <- apply(s, 1, which.min)
    y_min <- sapply(1:nrow(y), function(i) y[i, idx[i]])
    risk <- mean(y_min > zk)
  }
  return(c(baseline=baseline, risk=risk, 
           arr=baseline-risk,
           rrr=(baseline-risk)/baseline))
}

rho1_l <- seq(-0.4, 0.4, 0.1)
# rho2_l <- seq(-0.5, 0, 0.1)
rho2 <- 0
result <- NULL
for (rho1 in rho1_l) {
  # for (rho2 in rho2_l) {
  #   rho2 <- ifelse(sign(rho2) != sign(rho1), -rho2, rho2)
  #   result <- rbind(result,
  #                 c(rho1, rho2,
  #                   env_gen(iter, n, h2, r2, K,
  #                           rho1/2, rho2/2, F)))
  # }
  
  result <- rbind(result,
                  c(rho1, rho2, 
                    env_gen(iter, n, h2, r2, K, 
                            rho1, rho2, F)))
}

result <- data.frame(result)
names(result)[1:2] <- c("rho1", "rho2")

# par(mfrow = c(1, 1))
# plt(rrr ~ rho1 | rho2, data=result,
#     # col = rho2,
#     legend=legend(title = ""),
#     xlab = expression(rho[1]), ylab = "rrr",
#     # ylim = c(min(results[3, ], results[6, ]),
#     # max(results[3, ], results[6, ])), pch=19, 
#     theme="classic")
# abline(h = 0, lty=2)
# 
# plt(arr ~ rho1 | rho2, data=result,
#     # col = rho2,
#     legend=legend(title = ""),
#     xlab = expression(rho[1]), ylab = "arr",
#     # ylim = c(min(results[3, ], results[6, ]),
#     # max(results[3, ], results[6, ])), pch=19, 
#     theme="classic")

par(mfrow=c(1, 2))
plt(rrr ~ rho1, data=result,
    # col = rho2,
    # legend=legend(title = ""),
    xlab = expression(rho), ylab = "rrr",
    # ylim = c(min(results[3, ], results[6, ]),
    # max(results[3, ], results[6, ])), pch=19, 
    theme="classic")
abline(h = 0, lty=2)

plt(arr ~ rho1, data=result,
    # col = rho2,
    # legend=legend(title = ""),
    xlab = expression(rho), ylab = "arr",
    # ylim = c(min(results[3, ], results[6, ]),
    # max(results[3, ], results[6, ])), pch=19, 
    theme="classic")
abline(h = 0, lty=2)

# Conditional
result <- NULL
for (rho1 in rho1_l) {
  # for (rho2 in rho2_l) {
  #   rho2 <- ifelse(sign(rho2) != sign(rho1), -rho2, rho2)
  #   result <- rbind(result,
  #                   c(rho1, rho2, 
  #                     env_gen(iter, n, h2, r2, K, 
  #                             rho1/2, rho2/2, T)))
  # }
  result <- rbind(result,
                  c(rho1, rho2, 
                    env_gen(iter, n, h2, r2, K, 
                            rho1, rho2, T)))
}

# result <- data.frame(result)
# names(result)[1:2] <- c("rho1", "rho2")
# 
# tpar(mfrow = c(1, 2))
# plt(rrr ~ rho1 | rho2, data=result, 
#     # col = rho2,
#     legend=legend(title = ""),
#     xlab = expression(rho[1]), ylab = "rrr",
#     # ylim = c(min(results[3, ], results[6, ]),
#     # max(results[3, ], results[6, ])), pch=19, 
#     theme="classic")
# abline(h = 0, lty=2)
# 
# plt(arr ~ rho1 | rho2, data=result,
#     # col = rho2,
#     legend=legend(title = ""),
#     xlab = expression(rho[1]), ylab = "arr",
#     # ylim = c(min(results[3, ], results[6, ]),
#     # max(results[3, ], results[6, ])), pch=19, 
#     theme="classic")
# abline(h = 0, lty=2)

result <- data.frame(result)
names(result)[1:2] <- c("rho1", "rho2")

tpar(mfrow = c(1, 2))
plt(rrr ~ rho1, data=result, 
    # col = rho2,
    # legend=legend(title = ""),
    xlab = expression(rho), ylab = "rrr",
    # ylim = c(min(results[3, ], results[6, ]),
    # max(results[3, ], results[6, ])), pch=19, 
    theme="classic")
abline(h = 0, lty=2)

plt(arr ~ rho1, data=result,
    # col = rho2,
    xlab = expression(rho), ylab = "arr",
    # ylim = c(min(results[3, ], results[6, ]),
    # max(results[3, ], results[6, ])), pch=19, 
    theme="classic")
abline(h = 0, lty=2)