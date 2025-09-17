# library(tmvtnorm)

#' Risk reduction using the lowest PRS strategy
#' 
#' Calculate the relative risk reduction of the lowest PRS embryo based on the 
#' PRS accuracy ($r^2$), disease prevalence (K) and number of embryos (n).
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param n The number of live births. Any integer >= 1, though with 1 of course
#' the risk reduction is 0 (or close to it due to numerical reasons I didn't bother
#' handling separately).
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' risk_reduction_lowest(0.05, 0.01, 5)
#'
#' @export
risk_reduction_lowest = function(r2,K,n)
{
  stopifnot(r2 >=0, r2 <=1, K > 0,  K < 1, n >=1)
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  integrand_lowest = function(t)
  {
    arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  # risk/K, so if we want overall error < 1e-5, we need risk/K < 1e-5, so the integral error < K * 1e-5
  # Let's do 4 digits?
  risk = integrate(integrand_lowest,-Inf,Inf, abs.tol = K * 1e-4)$value
  reduction = (K-risk)/K
  # If abs risk: K-risk
  return(reduction)
}

# risk_reduction_exclude_old = function(r2,K,q,n)
# {
#   stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, q > 0, q < 1, n >=1)
#   r = sqrt(r2)
#   zk = qnorm(K, lower.tail=F)
#   zq = qnorm(q, lower.tail=F)
#   integrand_t = function(t,u)
#   {
#     y = dnorm(t)*pnorm((zk-r/sqrt(2)*(u+t))/sqrt(1-r2),lower.tail=F)
#     return(y)
#   }
#   # Replacement for integrand_t, should be numerically better.
#   inner_integral <- function(x, a, b) {
#     if (is.infinite(x)) {
#       OwenT(x, a/(x * sqrt(1+b^2))) + OwenT(a/sqrt(1+b^2), (x * sqrt(1+b^2)) / a) -
#         # OwenT(x, b) - OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
#         OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
#     }
#     else {
#       OwenT(x, a/(x * sqrt(1+b^2))) + OwenT(a/sqrt(1+b^2), (x * sqrt(1+b^2)) / a) -
#         OwenT(x, (a+b*x)/x) - OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
#     }
#   }
#   
#   integrand_u = function(us)
#   {
#     y = numeric(length(us))
#     beta_vec = zq*sqrt(2)-us
#     denom = pnorm(beta_vec)
#     dnorm_u <- dnorm(us)
#     denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
#     
#     a <- (zk - r/sqrt(2) * us) / sqrt(1-r2)
#     b <- -(r/sqrt(2)) / sqrt(1-r2)
#     
#     for (i in seq_along(us))
#     {
#       u = us[i]
#       beta <- beta_vec[i]
#       # internal_int = integrate(integrand_t,-Inf,beta,u, abs.tol = 1e-10)$value
#       
#       internal_int = ifelse(r2 != 1, pnorm(beta) - inner_integral(beta, a[i], b) + inner_integral(-Inf, a[i], b),
#                             ifelse(beta < sqrt(2) / r * zk - u, 0, pnorm(beta) - pnorm(sqrt(2) / r * zk - u)))
#       numer = dnorm_u[i]*(1-(1-denom[i])^n) * internal_int
#       term1 = numer/denom[i]
#       
#       # Probably doesn't matter, but also numerically more stable
#       # term1 <- exp(dnorm(u, log = T) + VGAM::log1mexp(-n*pnorm(beta, lower.tail = F, log.p = T)) -
#       #                pnorm(beta, log.p = T)) * internal_int
#       
#       # internal_int = integrate(integrand_t,beta,Inf,u, abs.tol = 1e-10)$value
#       internal_int <- ifelse(r2 != 1, 1 - pnorm(beta) - inner_integral(Inf, a[i], b) + inner_integral(beta, a[i], b),
#                              ifelse(beta < sqrt(2) / r * zk - u, 1-pnorm(sqrt(2) / r * zk - u), 1-pnorm(beta)))
#       term2 = dnorm_u[i]*(1-denom[i])^(n-1) * internal_int
#       y[i] = term1 + term2
#       # print(sprintf("%f", y[i]))
#     }
#     return(y)
#   }
#   risk = integrate(integrand_u,-Inf,Inf, abs.tol = K * 1e-4)$value
#   reduction = (K-risk)/K
#   # K-risk
#   return(reduction)
# }

#' Risk reduction using the exclude strategy
#' 
#' Calculate the relative risk reduction of a randomly selected embryo out
#' of embryos with PRS lower than a specific quantile (q), as a function of
#' PRS accuracy ($r^2$), disease prevalence (K) and number of embryos (n).
#' *In case there are no embryos to select from, a random one is selected.*
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param q The threshold quantile, in range of 0-1. When possible, embryos with PRS greater than
#' then 1-q quantile are excluded from the selection, and a random embryo
#' is selected from the non-excluded ones. If there are no embryos after the
#' exclusion, one embryo is selected at random.
#' @param n The number of live births. Any integer >= 1, though with 1 of course
#' the risk reduction is 0 (or close to it due to numerical reasons I didn't bother
#' handling separately).
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' risk_reduction_exclude(0.05, 0.01, 0.05, 5) # Exclude top 5% risk.
#'
#' @export
risk_reduction_exclude <- function(r2, K, q, n) {
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, q > 0, q < 1, n >=1)
  
  r <- sqrt(r2)
  zk <- qnorm(K, lower.tail = FALSE)
  zq <- qnorm(q, lower.tail = FALSE)
  
  integrand_u <- function(us) {
    us <- as.numeric(us)
    beta_vec <- zq * sqrt(2) - us
    denom <- pnorm(beta_vec)
    denom_safe <- ifelse(denom == 0, 1e-300, denom)
    dnorm_u <- dnorm(us)
    
    a <- (zk - r / sqrt(2) * us) / sqrt(1 - r2)
    b <- -(r / sqrt(2)) / sqrt(1 - r2)
    
    if (r2 != 1) {
      inner_beta <- inner_integral_vec_rcpp(beta_vec, a, b) 
      p_beta <- pnorm(beta_vec)
      
      # precompute inner integrals at ±Inf once per element
      inner_inf_neg <- inner_integral_vec_rcpp(rep(-Inf, length(a)), a, b)
      inner_inf_pos <- inner_integral_vec_rcpp(rep(Inf, length(a)), a, b)
      
      internal_int1 <- p_beta - inner_beta + inner_inf_neg
      internal_int2 <- 1 - p_beta - inner_inf_pos + inner_beta
      
      numer1 <- dnorm_u * (1 - (1 - denom)^n) * internal_int1
      term1 <- numer1 / denom_safe
      
      term2 <- dnorm_u * (1 - denom)^(n - 1) * internal_int2
      
      y <- term1 + term2
    } else {
      # r2 == 1 special case vectorised
      thr <- sqrt(2) / r * zk - us
      p_thr <- pnorm(thr)
      
      # internal_int1: if beta < thr then 0 else p_beta - p_thr
      p_beta <- pnorm(beta_vec)
      internal_int1 <- ifelse(beta_vec < thr, 0, p_beta - p_thr)
      internal_int2 <- ifelse(beta_vec < thr, 1 - p_thr, 1 - p_beta)
      
      numer1 <- dnorm_u * (1 - (1 - denom)^n) * internal_int1
      term1 <- numer1 / denom_safe
      term2 <- dnorm_u * (1 - denom)^(n - 1) * internal_int2
      y <- term1 + term2
    }
    
    y
  }
  
  risk <- integrate(integrand_u, -Inf, Inf, abs.tol = K * 1e-4)$value
  reduction <- (K - risk) / K
  reduction
}

#' Risk reduction using the lowest PRS strategy given the parents' PRS
#' 
#' Calculate the relative risk reduction of the lowest PRS embryo, 
#' when the parents' PRS is known (qf and qm), based on the 
#' PRS accuracy ($r^2$), disease prevalence (K) and number of embryos (n).
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param n The number of live births. Any integer >= 1, though with 1 of course
#' the risk reduction is 0 (or close to it due to numerical reasons I didn't bother
#' handling separately).
#' @param qf The father's PRS quantile, in range of 0-1. 
#' @param qm The mother's PRS quantile, in range of 0-1.
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' # The risk when one parent is in the top 5% of PRS, and the second in the median.
#' risk_reduction_lowest_conditional(0.05, 0.01, 5, 0.95, 0.5) 
#'
#' @export
risk_reduction_lowest_conditional = function(r2,K,n,qf,qm)
{
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, n >=1,
            qf > 0, qf < 1,
            qm > 0, qf < 1)
  
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  c = (zqf+zqm)/2 * r
  
  baseline = pnorm((zk-c)/sqrt(1-r2/2),lower.tail=F)
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  # risk = integrate(integrand_lowest_cond,-Inf,Inf,rel.tol=1e-9)$value
  risk = integrate(integrand_lowest_cond,-Inf,Inf, abs.tol = baseline * 1e-4)$value
  reduction = (baseline-risk)/baseline
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

#' Risk reduction using the exclude strategy, conditional on the parents' PRS.
#' 
#' Calculate the relative risk reduction of a randomly selected embryo out
#' of embryos with PRS lower than a specific quantile (q), when the parents' PRS is known (qf and qm), 
#' as a function of PRS accuracy ($r^2$), disease prevalence (K) and number of embryos (n).
#' *In case there are no embryos to select from, a random one is selected.*
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param q The threshold quantile, in range of 0-1. When possible, embryos with PRS greater than
#' then 1-q quantile are excluded from the selection, and a random embryo
#' is selected from the non-excluded ones. If there are no embryos after the
#' exclusion, one embryo is selected at random.
#' @param qf The father's PRS quantile, in range of 0-1. 
#' @param qm The mother's PRS quantile, in range of 0-1.
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' # The risk when one parent is in the top 5% of PRS, and the second in the median.
#' risk_reduction_lowest_conditional(0.05, 0.01, 5, 0.95, 0.5) 
#'
#' @export
risk_reduction_exclude_conditional = function(r2,K,q,n,qf,qm)
{
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, q > 0, q < 1, n >=1,
            qf > 0, qf < 1,
            qm > 0, qm < 1)
  
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  c = (zqf+zqm)/2 * r
  baseline= pnorm((zk-c)/sqrt(1-r^2/2),lower.tail=F)
  gamma = zq*sqrt(2) - c/(r/sqrt(2))
  
  integrand_t = function(t)
  {
    y = dnorm(t)*pnorm((zk-t*r/sqrt(2)-c)/sqrt(1-r2),lower.tail=F)
    return(y)
  }
  
  denom = pnorm(gamma)
  err <- (1e-4 * baseline * denom) / (2 * (1-(1-denom)^n))
  internal_int = integrate(integrand_t,-Inf,gamma, abs.tol = err)$value
  numer = (1-(1-denom)^n) * internal_int
  term1 = numer/denom
  
  err <- (1e-4 * baseline) / (2 * (1-(1-denom)^n))
  internal_int = integrate(integrand_t,gamma,Inf)$value
  term2 = (1-denom)^(n-1) * internal_int
  
  risk = term1 + term2
  
  reduction = (baseline-risk)/baseline
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

# Not used, but oh well.
simulate_lowest_risk_two_traits = function(r2A,r2B,rho,KA,KB,ns,nfam=10000)
{
  # results = array(0,c(length(ns),4))
  results <- matrix(nrow = length(ns), ncol = 4)
  
  disease_countA = numeric(length(ns))
  disease_count_randomA = numeric(length(ns))
  disease_countB = numeric(length(ns))
  disease_count_randomB = numeric(length(ns))
  
  baselineA = numeric(length(ns))
  baselineB = numeric(length(ns))
  
  t_sickA = qnorm(1-KA)
  t_sickB = qnorm(1-KB)
  
  rA = sqrt(r2A)
  rB = sqrt(r2B)
  
  mu = c(0,0)
  sigma = 1/2 * matrix(c(r2A,rho*rA*rB,rho*rA*rB,r2B),nrow=2)
  
  # chol_sigma <- chol(sigma)
  cs <- matrix(0, nrow = nfam, ncol = 2)
  
  for (i in seq_along(ns))
  {
    # cat('\r',i)
    
    n = ns[i]
    rmvn(nfam, mu, sigma, A = cs)
    
    xs <- rmvn(nfam*n, mu, sigma)
    xsA = matrix(xs[,1],nrow=nfam)
    xsB = matrix(xs[,2],nrow=nfam)
    
    envsA = rnorm(nfam,0,sqrt(1-r2A))
    envsB = rnorm(nfam,0,sqrt(1-r2B))
    
    scoresA <- xsA + cs[,1]#csA
    scoresB <- xsB + cs[,2]#csB
    
    selected_ind <- max.col(-scoresA)
    
    liabA <- scoresA[cbind(1:nfam, selected_ind)] + envsA
    liabB <- scoresB[cbind(1:nfam, selected_ind)] + envsB
    disease_countA <- sum(liabA > t_sickA)
    disease_countB <- sum(liabB > t_sickB)
    
    liab_randomA <- scoresA[, 1] + envsA
    disease_count_randomA <- sum(liab_randomA > t_sickA)
    liab_randomB <- scoresB[, 1] + envsB
    disease_count_randomB <- sum(liab_randomB > t_sickB)
    
    rrrA = (disease_count_randomA - disease_countA)/disease_count_randomA
    arrA = (disease_count_randomA - disease_countA)/nfam
    rrrB = (disease_count_randomB - disease_countB)/disease_count_randomB
    arrB = (disease_count_randomB - disease_countB)/nfam
    
    results[i,1] = rrrA
    results[i,2] = arrA
    results[i,3] = rrrB
    results[i,4] = arrB
  }  
  return(results)
}

# I can probably refactor the next three function, since they are almost identical
# to the non random case, just with a slightly different inner integral.

#' Risk reduction using the lowest PRS strategy, when the number of live births is binomial
#' 
#' Calculate the relative risk reduction of the lowest PRS embryo based on the 
#' PRS accuracy ($r^2$), disease prevalence (K), number of embryos (n) and probability
#' of live birth (p).
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param n The number of live births. Any integer >= 1, though with 1 of course
#' the risk reduction is 0 (or close to it due to numerical reasons I didn't bother
#' handling separately).
#' @param p The probability of live birth, in range of 0-1.
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' risk_reduction_lowest_bin(0.05, 0.01, 10, 0.5)
#'
#' @export
risk_reduction_lowest_bin <- function(r2, K, n, p)
{
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, n >=1,
            p > 0, p <= 1)
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  integrand_lowest = function(t)
  {
    arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
    temp <- ((1-p*pnorm(arg))^n - (1-p)^n) / (1-(1-p)^n)
    y = dnorm(t)*temp
    return(y)
  }
  # risk/K, so if we want overall error < 1e-5, we need risk/K < 1e-5, so the integral error < K * 1e-5
  # Let's do 4 digits?
  risk = integrate(integrand_lowest,-Inf,Inf, abs.tol = K * 1e-4)$value
  reduction = (K-risk)/K
  # If abs risk: K-risk
  return(reduction)
}

#' Risk reduction using the lowest PRS strategy, when the number of live births is Poisson
#' 
#' Calculate the relative risk reduction of the lowest PRS embryo based on the 
#' PRS accuracy ($r^2$), disease prevalence (K) and expected number of live births (lambda).
#' In the calculator we use lambda = n * p, where p is the probability of live birth and
#' n the expected number of embryos.
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param lambda The expected number of live births. Any integer >= 1.
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' risk_reduction_lowest_pois(0.05, 0.01, 5)
#'
#' @export
risk_reduction_lowest_pois <- function(r2, K, lambda)
{
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, lambda > 0)
  
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  integrand_lowest = function(t)
  {
    arg = (zk-t*sqrt(1-r2/2)) / (r/sqrt(2))
    temp <- (exp(lambda*(-pnorm(arg))) - exp(-lambda)) / (1-exp(-lambda))
    y = dnorm(t)*temp
    return(y)
  }
  # risk/K, so if we want overall error < 1e-5, we need risk/K < 1e-5, so the integral error < K * 1e-5
  # Let's do 4 digits?
  risk = integrate(integrand_lowest,-Inf,Inf, abs.tol = K * 1e-4)$value
  reduction = (K-risk)/K
  # If abs risk: K-risk
  return(reduction)
}

#' Risk reduction using the lowest PRS strategy given the parents' PRS, when the number of live births is binomial
#' 
#' Calculate the relative risk reduction of the lowest PRS embryo, 
#' when the parents' PRS is known (qf and qm), based on the 
#' PRS accuracy ($r^2$), disease prevalence (K) and number of embryos (n) and probability
#' of live birth (p).
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param n The number of live births. Any integer >= 1, though with 1 of course
#' the risk reduction is 0 (or close to it due to numerical reasons I didn't bother
#' handling separately).
#' @param qf The father's PRS quantile, in range of 0-1. 
#' @param qm The mother's PRS quantile, in range of 0-1.
#' @param p The probability of live birth, in range of 0-1.
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' # The risk when one parent is in the top 5% of PRS, and the second in the median.
#' risk_reduction_lowest_conditional_bin(0.05, 0.01, 10, 0.95, 0.5, 0.5)
#'
#' @export
risk_reduction_lowest_conditional_bin = function(r2,K,n,qf,qm,p)
{
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, n >=1,
            qf > 0, qf < 1,
            qm > 0, qm < 1,
            p > 0, p <= 1)
  
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  c = (zqf+zqm)/2 * r
  
  baseline = pnorm((zk-c)/sqrt(1-r2/2),lower.tail=F)
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r2)) / (r/sqrt(2))
    temp <- ((1-p*pnorm(arg))^n - (1-p)^n) / (1-(1-p)^n)
    y = dnorm(t)*temp
    return(y)
  }
  # risk = integrate(integrand_lowest_cond,-Inf,Inf,rel.tol=1e-9)$value
  risk = integrate(integrand_lowest_cond,-Inf,Inf, abs.tol = baseline * 1e-4)$value
  reduction = (baseline-risk)/baseline
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

#' Risk reduction using the lowest PRS strategy given the parents' PRS, when the expected number of live births is poisson
#' 
#' Calculate the relative risk reduction of the lowest PRS embryo, 
#' when the parents' PRS is known (qf and qm), based on the 
#' PRS accuracy ($r^2$), disease prevalence (K) and expected number of live births (lambda).
#' In the calculator we use lambda = n * p, where p is the probability of live birth and
#' n the expected number of embryos.
#' 
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param lambda The expected number of live births. Any integer >= 1.
#' @param qf The father's PRS quantile, in range of 0-1. 
#' @param qm The mother's PRS quantile, in range of 0-1.
#' 
#' @return The relative risk reduction.
#' 
#' @examples
#' # The risk when one parent is in the top 5% of PRS, and the second in the median.
#' risk_reduction_lowest_conditional_pois(0.05, 0.01, 5, 0.95, 0.5)
#'
#' @export
risk_reduction_lowest_conditional_pois = function(r2,K,lambda,qf,qm)
{
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, lambda > 0,
            qf > 0, qf < 1,
            qm > 0, qm < 1)
  
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)

  c = (zqf+zqm)/2 * r
  
  baseline = pnorm((zk-c)/sqrt(1-r2/2),lower.tail=F)
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r2)) / (r/sqrt(2))
    temp <- (exp(lambda*(-pnorm(arg))) - exp(-lambda)) / (1-exp(-lambda))
    y = dnorm(t)*temp
    return(y)
  }
  # risk = integrate(integrand_lowest_cond,-Inf,Inf,rel.tol=1e-9)$value
  risk = integrate(integrand_lowest_cond,-Inf,Inf, abs.tol = baseline * 1e-4)$value
  reduction = (baseline-risk)/baseline
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

# Next are functions for when we condition on family history.

# Sample from truncated binomial/poisson, since we condition on n>0
sample_truncated_binomial <- function(num_samples, n, p) {
  pmf_vals <- dbinom(0:n, size = n, prob = p)
  prob_gt_zero <- 1 - pmf_vals[1]
  
  if (prob_gt_zero <= 1e-10) {
    warning("Probability of non-zero outcome is low. Returning 1.")
    return(rep(1, num_samples))
  }
  
  # Normalize the probabilities and calculate the CDF
  truncated_pmf <- pmf_vals[-1] / prob_gt_zero
  truncated_cdf <- cumsum(truncated_pmf)
  
  u <- runif(num_samples)
  
  # For each u, find the first index in the CDF that is >= u.
  # We add 1 because findInterval returns a 0-based index into the k=1...n space.
  samples <- findInterval(u, truncated_cdf) + 1
  return(samples)
}

sample_truncated_poisson <- function(num_samples, lambda) {
  prob_gt_zero <- 1 - exp(-lambda)
  
  if (prob_gt_zero <= 1e-10) {
    warning("Lambda is very small; probability of non-zero outcome is negligible. All samples will be 1.")
    return(rep(1, num_samples))
  }
  
  # Determine a safe upper limit for k to calculate the CDF.
  # qpois gives the quantile; we go far into the tail to ensure our CDF is ~1.
  # We add a buffer (e.g., +5) just to be extra safe.
  k_max <- qpois(1 - 1e-15, lambda) + 5 
  
  # Normalize the probabilities and calculate the CDF
  pmf_vals <- dpois(0:k_max, lambda)
  truncated_pmf <- pmf_vals[-1] / prob_gt_zero
  truncated_cdf <- cumsum(truncated_pmf)
  
  u <- runif(num_samples)
  
  # For each u, find the first index in the CDF that is >= u.
  samples <- findInterval(u, truncated_cdf) + 1
  return(samples)
}

# Sample from normal distribution with linear inequalities.
# i.e. the n samples x are so Gx > zk (or <, depending on directions)
# Could probably also use the function from tmvtnorm.
#' @import tmvtnorm
#' @import mvnfast
sample_func <- function(n, zk, G, Sigma, directions) {
  d <- nrow(Sigma)
  if (ncol(G) != d) {stop("The number of cols in G should be the same as the dimension of Sigma.")}
  if (length(directions) != length(zk)) {stop("The vector directions should be the same length as the vector zk.")}
  if (!all(directions %in% c(">", "<"))) {stop("directions should be one of \">\", or \"<\".")}
  
  # chol_sigma <- chol(Sigma)
  # Using chol here for G_chol_sigma sometimes breaks things for 
  # numerical reasons, though from my tests, not in the cases we care about.
  # G_chol_sigma <- tcrossprod(G %*% chol_sigma) # This is G %*% Sigma %*% t(G)
  # G_chol_sigma <- G %*% Sigma %*% t(G)
  tcross_Sigma_G <- tcrossprod(Sigma, G)
  G_chol_sigma <- G %*% tcross_Sigma_G
  temp <- tcross_Sigma_G %*% solve(G_chol_sigma)
  
  # Sample the truncated liabilities 'r'
  r <- tmvtnorm::rtmvnorm(n, sigma = G_chol_sigma,
                          lower = ifelse(directions == ">", zk, -Inf),
                          upper = ifelse(directions == "<", zk, Inf),
                          algorithm = "gibbs")
  
  # My function is a bit faster
  # r <- gibbs_sampling_truncated_normal(n, G_chol_sigma,
  #                                      ifelse(directions == ">", zk, -Inf),
  #                                      ifelse(directions == "<", zk, Inf))
  r <- t(r)
  
  # Sample an unconditional y ~ N(0, Sigma)
  # y <- rmvn(n, rep(0, d), sigma = chol_sigma, isChol = T)
  y <- rmvn(n, rep(0, d), sigma = Sigma, isChol = F)
  
  # Correct y to get a sample from the conditional distribution P(z | Gz=r)
  correction <- temp %*% (r - tcrossprod(G, y))
  y + t(correction)
}

#' Risk reduction when conditioning on the family disease status.
#' 
#' A generic function which estimates the risk in either the lowest PRS embryo
#' or high-risk exclusion strategy, and whether we use fixed number of live births,
#' binomial or poisson.
#' 
#' @param iter The number of monte-carlo samples to use in our estimate.
#' @param n Depending on the "random_strategy" variable, either the number of live births,
#' number of embryos or expected number of embryos.
#' @param r2 The PRS accuracy, in range of 0-1.
#' @param h2 The disease heritability, in range of 0-1 and should be greater than _r2_.
#' @param K The disease prevalence, in range of 0-1. Not including 0 and 1 since
#' the extremes are not really interesting.
#' @param sick_parents The number of parents which are known to be with the disease.
#' @param no_sick_parents The number of parents which are known not to have the disease.
#' @param sick_siblings The number of siblings which are known to be with the disease.
#' @param no_sick_siblings The number of siblings which are known not to have the disease.
#' @param qf The father's PRS quantile, in range of 0-1. NULL if unknown.
#' @param qm The mother's PRS quantile, in range of 0-1. NULL if unknown.
#' @param selection_strategy Either "lowest_prs" or "exclude_percentile".
#' If exclude_percentile, should also provide exclusion_q.
#' @param exclusion_q The threshold quantile, in range of 0-1. When possible, embryos with PRS greater than
#' then 1-q quantile are excluded from the selection, and a random embryo
#' is selected from the non-excluded ones. If there are no embryos after the
#' exclusion, one embryo is selected at random. NULL if the strategy is "lowest_prs".
#' @param random_strategy One of "Fixed", "Binomial" or "Poisson".
#' @importFrom matrixStats rowMins
#' @param p The probability of live birth (relevant only if random_strategy is either Binomial or Poisson), in range of 0-1.
risk_parents_offspring_generic <- function(iter, n, r2, h2, K,
                                           sick_parents = 0,
                                           no_sick_parents = 0,
                                           sick_siblings = 0,
                                           no_sick_siblings = 0,
                                           qf = NULL, qm = NULL,
                                           selection_strategy = "lowest_prs",
                                           exclusion_q = NULL,
                                           random_strategy = "Fixed",
                                           p=0) {
  
  # Some validations
  stopifnot(
    h2 >= r2, h2 <= 1, r2 >= 0, K > 0, K < 1, iter > 0, n > 0,
    all(c(sick_parents, no_sick_parents, sick_siblings, no_sick_siblings) >= 0),
    sick_parents + no_sick_parents <= 2
  )
  if (!is.null(qf) && is.null(qm) || is.null(qf) && !is.null(qm)) {
    stop("Please provide both qf and qm, or neither.")
  }
  if (selection_strategy == "exclude_percentile" && is.null(exclusion_q)) {
    stop("For 'exclude_percentile' strategy, you must provide 'exclusion_q'.")
  }
  if (!is.null(exclusion_q)) stopifnot(exclusion_q > 0, exclusion_q < 1)
  if (random_strategy != "Fixed") stopifnot(p > 0, p <= 1)
  
  zk <- qnorm(1 - K)
  n_known_parents <- sick_parents + no_sick_parents
  n_known_sibs <- sick_siblings + no_sick_siblings
  
  s_p_mean <- NULL
  w_p_mean <- NULL
  
  # Handle the r2 == h2 case explicitly.
  # The w component variance is 0, so it's a constant, not a random variable.
  is_w_zero <- abs(h2 - r2) < 1e-9
  
  # Sample the parental scores and non-score component based on the family
  # history
  
  # Not really pretty.
  if (is.null(qf)) {
    # Parental PRS is unknown
    if (n_known_parents > 0) {
      # Z = (s_f, [w_f], e_f, s_m, [w_m], e_m, sib_specific_1, ...)
      # The [w] components are included only if their variance is non-zero.
      
      parent_vars <- if (is_w_zero) c(r2, 1-h2) else c(r2, h2-r2, 1-h2)
      cols_per_parent <- length(parent_vars)
      s_idx <- 1 
      w_idx <- if(is_w_zero) NULL else 2
      
      n_vars <- 2 * cols_per_parent + n_known_sibs
      Sigma <- diag(c(rep(parent_vars, 2), rep(1-h2/2, n_known_sibs)))
      
      G_parents <- matrix(0, nrow=n_known_parents, ncol=n_vars)
      if (n_known_parents >= 1) G_parents[1, 1:cols_per_parent] <- 1
      if (n_known_parents == 2) G_parents[2, (cols_per_parent+1):(2*cols_per_parent)] <- 1
      
      G_sibs <- matrix(0, nrow=n_known_sibs, ncol=n_vars)
      if (n_known_sibs > 0) {
        sib_shared_indices <- c(s_idx, w_idx, s_idx + cols_per_parent, w_idx + cols_per_parent)
        G_sibs[, sib_shared_indices] <- 0.5
        G_sibs[, (2*cols_per_parent + 1):n_vars] <- diag(n_known_sibs)
      }
      G <- rbind(G_parents, G_sibs)
      
      samples <- sample_func(iter, rep(zk, n_known_parents + n_known_sibs), G, Sigma, 
                             c(rep(">", sick_parents), rep("<", no_sick_parents),
                               rep(">", sick_siblings), rep("<", no_sick_siblings)))
      
      s_p_mean <- (samples[, s_idx] + samples[, s_idx + cols_per_parent]) / 2
      w_p_mean <- if (is_w_zero) 0 else (samples[, w_idx] + samples[, w_idx + cols_per_parent]) / 2
      
    } else if (n_known_sibs > 0) {
      # Z = (s_shared, [w_shared], sib_specific_1, ...)
      shared_vars <- if (is_w_zero) c(r2/2) else c(r2/2, (h2-r2)/2)
      s_idx <- 1
      w_idx <- if(is_w_zero) NULL else 2
      
      n_vars <- length(shared_vars) + n_known_sibs
      Sigma <- diag(c(shared_vars, rep(1-h2/2, n_known_sibs)))
      
      G <- cbind(matrix(1, nrow=n_known_sibs, ncol=length(shared_vars)), diag(n_known_sibs))
      
      samples <- sample_func(iter, rep(zk, n_known_sibs), G, Sigma, 
                             c(rep(">", sick_siblings), rep("<", no_sick_siblings)))
      
      s_p_mean <- samples[, s_idx]
      w_p_mean <- if (is_w_zero) 0 else samples[, w_idx]
    }
  } else {
    # Parental PRS is known
    s_f_known <- sqrt(r2) * qf
    s_m_known <- sqrt(r2) * qm
    s_p_mean <- rep((s_f_known + s_m_known) / 2, iter) # Constant vector
    
    if (is_w_zero) {
      w_p_mean <- 0 # If w var is 0, its mean contribution is 0. Nothing to sample.
    } else {
      # Only sample w and e components if w variance > 0
      if (n_known_parents > 0) {
        # Z = (w_f, e_f, w_m, e_m, sib_specific_1, ...)
        n_vars <- 4 + n_known_sibs
        Sigma <- diag(c(h2-r2, 1-h2, h2-r2, 1-h2, rep(1-h2/2, n_known_sibs)))
        
        zk_adj_parents <- c(if(sick_parents>0) rep(zk-s_f_known, sick_parents) else NULL,
                            if(no_sick_parents>0) rep(zk-s_m_known, no_sick_parents) else NULL)
        zk_adj_sibs <- rep(zk - s_p_mean[1], n_known_sibs)
        
        G_parents <- matrix(0, nrow=n_known_parents, ncol=n_vars)
        if (n_known_parents >= 1) G_parents[1, 1:2] <- 1
        if (n_known_parents == 2) G_parents[2, 3:4] <- 1
        
        G_sibs <- matrix(0, nrow=n_known_sibs, ncol=n_vars)
        if (n_known_sibs > 0) {
          G_sibs[, c(1,3)] <- 0.5 # 0.5*w_f + 0.5*w_m
          G_sibs[, 5:n_vars] <- diag(n_known_sibs)
        }
        G <- rbind(G_parents, G_sibs)
        
        samples <- sample_func(iter, c(zk_adj_parents, zk_adj_sibs), G, Sigma, 
                               c(rep(">", sick_parents), rep("<", no_sick_parents),
                                 rep(">", sick_siblings), rep("<", no_sick_siblings)))
        w_p_mean <- (samples[, 1] + samples[, 3]) / 2
        
      } else if (n_known_sibs > 0) {
        # Z = (w_shared, sib_specific_1, ...)
        Sigma <- diag(c((h2-r2)/2, rep(1-h2/2, n_known_sibs)))
        zk_adj_sibs <- rep(zk - s_p_mean[1], n_known_sibs)
        G <- cbind(1, diag(n_known_sibs))
        
        samples <- sample_func(iter, zk_adj_sibs, G, Sigma,
                               c(rep(">", sick_siblings), rep("<", no_sick_siblings)))
        w_p_mean <- samples[, 1]
      }
    }
  }
  
  if (is.null(s_p_mean)) {
    # No family history provided at all. Shouldn't happen in reality.
    # Maybe change it to use the default risk calculation? since
    # then the selection risk is not K.
    return(c(baseline=K, selection=K, rel_red=0, abs_red=0, sd=0))
  }
  
  # Baseline risk (random embryo)
  liability_mean_baseline <- s_p_mean + w_p_mean
  sd_baseline <- sqrt(1 - h2/2)
  baseline_risks <- pnorm(zk, mean = liability_mean_baseline, sd = sd_baseline, 
                          lower.tail = FALSE)
  baseline <- mean(baseline_risks)
  
  # Selection risk
  # The NA trick in the matrix is somewhat ugly, but oh well.
  if (random_strategy == "Binomial") {
    ns <- sample_truncated_binomial(iter, n, p)
    embryo_s_mendelian <- sapply(ns, function(ns) c(rnorm(ns, sd=sqrt(r2/2)), 
                                                    rep(NA, n-ns))) |>
      t()
  }
  else if (random_strategy == "Poisson") {
    ns <- sample_truncated_poisson(iter, n * p)
    # ns <- sample_truncated_binomial(iter, n, p)
    max_ns <- max(ns)
    embryo_s_mendelian <- sapply(ns, function(ns) c(rnorm(ns, sd=sqrt(r2/2)), 
                                                    rep(NA, max_ns-ns))) |>
      t()
  }
  else {
    embryo_s_mendelian <- matrix(rnorm(iter * n, sd = sqrt(r2/2)), nrow = iter, 
                                 ncol = n)
  }
  
  embryo_scores <- s_p_mean + embryo_s_mendelian
  
  # Select one embryo based on the chosen strategy
  if (selection_strategy == "lowest_prs") {
    selected_score <- matrixStats::rowMins(embryo_scores, na.rm = T)
  } else if (selection_strategy == "exclude_percentile") {
    prs_threshold <- qnorm(exclusion_q, mean = 0, sd = sqrt(r2)) # Threshold on the absolute PRS scale
    
    # selected_score <- sapply(1:iter, function(i) {
    #   eligible_scores <- embryo_scores[i, embryo_scores[i,] < prs_threshold]
    #   if (length(eligible_scores) > 0) {
    #     # Sample randomly from the eligible embryos
    #     return(sample(eligible_scores, 1))
    #   } else {
    #     # If no embryo meets the criteria, sample one from the full set
    #     return(sample(embryo_scores[i,], 1))
    #   }
    # })
    # eligible_scores <- embryo_scores < prs_threshold
    # lengths <- rowSums(eligible_scores)
    # selected_score <- sapply(1:iter, function(i) {
    #   if(lengths[i] > 0) {
    #     return(sample(embryo_scores[i, which(eligible_scores[i, ]==1)], 1))
    #   }
    #   else {
    #     return(sample(embryo_scores[i, ], 1))
    #   }
    # })
    selected_score <- selected_score_rcpp(embryo_scores, iter, prs_threshold)
  } else {
    stop("Invalid selection_strategy provided.")
  }
  
  liability_mean_selection <- selected_score + w_p_mean
  sd_selection <- sqrt((h2 - r2)/2 + (1 - h2)) # Or just 1-h2/2-r2/2
  selection_risks <- pnorm(zk, mean = liability_mean_selection, 
                           sd = sd_selection, lower.tail = FALSE)
  selection <- mean(selection_risks)
  
  # print(cov(selection_risks, baseline_risks))
  return(c(baseline = baseline,
           selection = selection,
           relative_reduction = (baseline - selection) / baseline,
           absolute_reduction = baseline - selection,
           sd_of_estimate = sd(selection_risks) / sqrt(iter),
           sd_of_baseline = sd(baseline_risks) / sqrt(iter)))
}
