#' Unified Risk Prediction
#' 
#' Calculates disease risk and reduction using numerical integration.
#' Supports population-level or parental-conditional estimates across all strategies.
#' 
#' @param r2 PRS accuracy (0-1).
#' @param K Disease prevalence (0-1).
#' @param n Number of embryos (Fixed) or N parameter (Binomial/Poisson).
#' @param selection_strategy "lowest_prs" or "exclude_percentile".
#' @param exclusion_q Percentile threshold for exclusion (required if strategy is exclude).
#' @param qf Father's PRS quantile (0-1). Optional.
#' @param qm Mother's PRS quantile (0-1). Optional.
#' @param random_strategy "Fixed", "Binomial", or "Poisson".
#' @param p_lb Probability of live birth (required for Binomial/Poisson).
#' 
#' @return A list: baseline (random risk), selection (strategy risk), rr (relative reduction).
#' @examples
#' # Expected when selecting the minimum PRS embryo
#' risk_prediction_analytical(0.05, 0.01, 5)
#' # Excluding the top 5% PRS
#' risk_prediction_analytical(0.05, 0.01, 5, 
#' selection_strategy = "exclude", 
#' exclusion_q = 0.05)
#' # The risk when one parent is in the top 5% of PRS, and the second in the median.
#' risk_prediction_analytical(0.05, 0.01, 5, 
#' selection_strategy = "lowest", 
#' qf = 0.95, qm = 0.5)
#' # The risk when one parent is in the top 5% of PRS, and the second in the median.
#' # and excluding embryos below the 50% percentile.
#' risk_prediction_analytical(0.05, 0.01, 5, 
#' selection_strategy = "exclude",
#' exclusion_q = 0.5,
#' qf = 0.95, qm = 0.5)
#' Same as the previous, but in the binomial/poisson case.
#' risk_prediction_analytical(0.05, 0.01, 10, 
#' random_strategy = "Binomial", 
#' p_lb = 0.5)
#' risk_prediction_analytical(0.05, 0.01, 10, 
#' random_strategy = "Poisson", 
#' p_lb = 0.5)
#' risk_prediction_analytical(0.05, 0.01, 10, qf = 0.95, qm = 0.5,
#' random_strategy = "Binomial", 
#' p_lb = 0.5)
#' risk_prediction_analytical(0.05, 0.01, 10, qf = 0.95, qm = 0.5,
#' random_strategy = "Poisson", 
#' p_lb = 0.5)
#' @export
risk_prediction_analytical <- function(r2,
                                       K,
                                       n,
                                       selection_strategy = c("lowest_prs", "exclude_percentile"),
                                       exclusion_q = NULL,
                                       qf = NULL,
                                       qm = NULL,
                                       random_strategy = c("Fixed", "Binomial", "Poisson"),
                                       p_lb = 1) {
  if (n < 1) stop("'n' must be a positive integer.")
  
  if (r2 < 0 || r2 > 1) stop("'r2' must be between 0 and 1.")
  if (K <= 0 || K >= 1) stop("'K' must be strictly between 0 and 1.")
  strategy <- match.arg(selection_strategy)
  model <- match.arg(random_strategy)
  
  r <- sqrt(r2)
  zk <- qnorm(K, lower.tail = FALSE)
  
  if (!is.null(qf) && !is.null(qm)) {
    # Conditional on parents PRS case
    c_val <- (qnorm(qf) + qnorm(qm)) / 2 * r
    baseline <- pnorm((zk - c_val) / sqrt(1 - r2 / 2), lower.tail = FALSE)
    zk_adj <- zk - c_val
    
    # SD of the non-score liability component
    sigma_e <- sqrt(1 - r2) 
  } else {
    # Unconditional case
    baseline <- K
    zk_adj <- zk
    c_val <- 0
    
    # SD of the non-score liability component
    sigma_e <- sqrt(1 - r2 / 2) 
  }
  
  # SD of the Mendelian/Score component (t)
  sigma_v <- r / sqrt(2) 
  
  if (strategy == "lowest_prs") {
    integrand_lowest <- function(t) {
      arg <- (zk_adj - t * sigma_e) / sigma_v
      p_sick_given_t <- pnorm(arg, lower.tail = FALSE)
      
      # Probability weight based on embryo model
      weight <- switch(
        model,
        "Fixed" = p_sick_given_t^n,
        "Binomial" = ((1 - p_lb * (
          1 - p_sick_given_t
        ))^n - (1 - p_lb)^n) / (1 - (1 - p_lb)^n),
        "Poisson" = {
          lambda <- n * p_lb
          (exp(lambda * (-(
            1 - p_sick_given_t
          ))) - exp(-lambda)) / (1 - exp(-lambda))
        }
      )
      dnorm(t) * weight
    }
    
    risk <- integrate(integrand_lowest, -Inf, Inf, 
                      abs.tol = baseline * 1e-4)$value
  }
  else {
    # High-risk exclusion case
    if (is.null(qf) && is.null(qm)) {
      # Just calling the already defined function since
      # it is long enough.
      RRR <- risk_reduction_exclude(r2, K, exclusion_q, n)
      risk <- K - RRR * K
    }
    else {
      # Threshold in terms of the score deviation t
      zq <- qnorm(exclusion_q, lower.tail = FALSE)
      gamma <- (zq * r - c_val) / sigma_v
      
      integrand_exclude <- function(t) {
        dnorm(t) * pnorm((zk_adj - t * sigma_v) / sigma_e, lower.tail = FALSE)
      }
      
      denom <- pnorm(gamma)
      # If no embryos pass exclusion, we pick one at random (integrated over all t)
      
      int_low <- integrate(integrand_exclude, -Inf, gamma, abs.tol = 1e-6)$value
      term1 <- (1 - (1 - denom)^n) * (int_low / denom)
      
      int_high <- integrate(integrand_exclude, gamma, Inf, abs.tol = 1e-6)$value
      term2 <- (1 - denom)^(n - 1) * int_high
      
      risk <- term1 + term2
    }
  }
  
  return(list(
    baseline = baseline,
    selection = risk,
    rr = (baseline - risk) / baseline
  ))
}

# Older versions
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
risk_reduction_lowest_conditional = function(r2,K,n,qf,qm)
{
  stopifnot(r2 >=0, r2 <=1, K > 0, K < 1, n >=1,
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
#' # and excluding embryos below the 50% percentile.
#' risk_reduction_exclude_conditional(0.05, 0.01, 0.5, 5, 0.95, 0.5)
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
    s_f_known <- sqrt(r2) * qnorm(qf)
    s_m_known <- sqrt(r2) * qnorm(qm)
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

#' @import tmvtnorm
#' @import mvnfast
sample_func_optimized <- function(n, G, Sigma, lower, upper) {
  d <- nrow(Sigma)
  
  tcross_Sigma_G <- tcrossprod(Sigma, G) 
  G_tcross_Sigma_G <- G %*% tcross_Sigma_G

  weight_mat <- solve(G_tcross_Sigma_G, G %*% Sigma)
  
  r <- tmvtnorm::rtmvnorm(n, sigma = G_tcross_Sigma_G, 
                          lower = lower, upper = upper,
                          algorithm = "gibbs")
  
  # Sometimes chol fails in rmvn, so that's a way to "fix" it
  # Should probably also add the ability to use multicore
  # chol_Sigma <- chol(Sigma, pivot=T)
  # y <- mvnfast::rmvn(n, mu = rep(0, d), sigma = chol_Sigma, 
  #                    isChol = T)
  y <- MASS::mvrnorm(n, rep(0, d), Sigma)
  Gy <- tcrossprod(y, G)
  diff <- r - Gy
  
  y + (diff %*% weight_mat)
}

# Construct the equality and inequality matrices
construct_G_eq <- function(prs_data, p) {
  n <- length(unlist(prs_data))
  G_eq <- matrix(0, nrow=n, ncol=p)
  row_index <- 1
  vals <- numeric(n)
  
  add_row <- function(cols, weights, val) {
    G_eq[row_index, cols] <<- weights
    vals[row_index]       <<- val
    row_index             <<- row_index + 1
  }
  
  # Indexes of the needed scores, and their weights
  fixed_defs <- list(
    gp1a = list(c(1), 1),
    gp1b = list(c(4), 1),
    gp2a = list(c(7), 1),
    gp2b = list(c(10), 1),
    p1   = list(c(1, 4, 13), c(0.5, 0.5, 1)),
    p2   = list(c(7, 10, 16), c(0.5, 0.5, 1))
  )
  
  for (name in names(fixed_defs)) {
    val <- prs_data[[name]]
    if (!is.null(val) && !is.na(val)) {
      def <- fixed_defs[[name]]
      add_row(def[[1]], def[[2]], val)
    }
  }
  
  dyn_col <- 19
  
  sib_configs <- list(
    sib_p1   = list(fix=c(1, 4),         w=c(0.5, 0.5, 1)),
    sib_p2   = list(fix=c(7, 10),        w=c(0.5, 0.5, 1)),
    sib_self = list(fix=c(1,4,7,10,13,16), w=c(rep(0.25, 4), 0.5, 0.5, 1))
  )
  
  for (grp in names(sib_configs)) {
    prss <- prs_data[[grp]]
    if (!is.null(prss)) {
      cfg <- sib_configs[[grp]]
      for (val in prss) {
        if (!is.na(val)) {
          add_row(c(cfg$fix, dyn_col), cfg$w, val)
        }
        dyn_col <- dyn_col + 3 
      }
    }
  }
  
  # Trim unused rows
  if (row_index > 1) {
    G_eq <- G_eq[1:(row_index-1), , drop=FALSE]
    vals <- vals[1:(row_index-1)]
  } else {
    G_eq <- matrix(0, nrow=0, ncol=p)
    vals <- numeric(0)
  }
  
  return(list(G_eq=G_eq, vals=vals))
}

construct_G_in <- function(history, p, zk) {
  n <- length(unlist(history))
  G_in <- matrix(0, nrow=n, ncol=p)
  lower_in <- numeric(n)
  upper_in <- numeric(n)
  
  row_index <- 1
  
  add_row <- function(cols, weights, status) {
    G_in[row_index, cols] <<- weights
    lower_in[row_index]   <<- ifelse(status == 1, zk, -Inf)
    upper_in[row_index]   <<- ifelse(status == 1, Inf, zk)
    row_index             <<- row_index + 1
  }
  
  fixed_defs <- list(
    gp1a = list(c(1:3), rep(1, 3)),
    gp1b = list(c(4:6), rep(1, 3)),
    gp2a = list(c(7:9), rep(1, 3)),
    gp2b = list(c(10:12), rep(1, 3)),
    p1   = list(c(1,2, 4,5, 13:15), c(rep(0.5, 4), rep(1, 3))),
    p2   = list(c(7,8, 10,11, 16:18), c(rep(0.5, 4), rep(1, 3)))
  )
  
  for (name in names(fixed_defs)) {
    stat <- history[[name]]
    if (!is.null(stat) && !is.na(stat) && stat != 0) {
      def <- fixed_defs[[name]]
      add_row(def[[1]], def[[2]], history[[name]])
    }
  }
  
  dyn_col <- 19
  
  sib_configs <- list(
    sib_p1   = list(fix=c(1,2, 4,5),               w=c(rep(0.5, 4), rep(1, 3))),
    sib_p2   = list(fix=c(7,8, 10,11),             w=c(rep(0.5, 4), rep(1, 3))),
    sib_self = list(fix=c(1,2, 4,5, 7,8, 10,11, 13,14, 16,17), w=c(rep(0.25, 8), rep(0.5, 4), rep(1, 3)))
  )
  
  for (grp in names(sib_configs)) {
    stats <- history[[grp]]
    if (!is.null(stats)) {
      cfg <- sib_configs[[grp]]
      for (stat in stats) {
        if (!is.na(stat) & stat != 0) {
          cols <- c(cfg$fix, dyn_col:(dyn_col+2))
          add_row(cols, cfg$w, stat)
        }
        dyn_col <- dyn_col + 3
      }
    }
  }
  
  # Trim
  if (row_index > 1) {
    G_in <- G_in[1:(row_index-1), , drop=FALSE]
    lower_in <- lower_in[1:(row_index-1)]
    upper_in <- upper_in[1:(row_index-1)]
  } else {
    G_in <- matrix(0, nrow=0, ncol=p)
    lower_in <- numeric(0)
    upper_in <- numeric(0)
  }
  
  return(list(G_in=G_in, lower_in=lower_in, upper_in=upper_in))
}

#' Estimate relative risk reduction conditional on family history
#' 
#' Estimates the relative and absolute risk reduction of embryo selection 
#' strategies by modeling the liability threshold model conditional on 
#' specific family disease history and/or known Polygenic Risk Scores (PRS).
#' 
#' The function simulates family liabilities using a multivariate normal distribution,
#' applies constraints based on `history` (inequality constraints) and `prs_data`
#' (equality constraints), and estimates the disease risk of a selected embryo
#' compared to a random embryo.
#' 
#' @param iter Integer. The number of Monte-Carlo samples to use for the estimate.
#' @param n Integer. The number of available embryos. If `random_strategy`
#'   is "Fixed", this is the exact number. If "Binomial" or "Poisson", this acts 
#'   as the parameter N (trials) or basis for lambda, respectively.
#' @param r2 Numeric (0-1). The variance explained by the PRS (SNP-heritability).
#' @param h2 Numeric (0-1). The total narrow-sense heritability. Must be >= `r2`.
#' @param K Numeric (0-1). The disease prevalence in the population.
#' @param history A named list of clinical statuses. Allowed names: 
#'   "p1", "p2", "gp1a", "gp1b", "gp2a", "gp2b", "sib_p1", "sib_p2", "sib_self".
#'   Values should be integer vectors: 1 (Sick), -1 (Healthy), or 0/NA (Unknown).
#' @param prs_data A named list of known PRS Z-scores (standardized). 
#'   Allowed names matches `history`. Values should be numeric vectors.
#'   Note: You can use NA in cases where the PRS is unknown. For example, suppose you provide
#'   the following history for the siblings, c(1, -1, 1). If you don't know the 
#'   PRS of the second sibling, you can use c(PRS_1, NA, PRS_2) for the PRS.
#' @param selection_strategy Character. Either "lowest_prs" (select embryo with 
#'   lowest score) or "exclude_percentile" (exclude high risk, then random).
#' @param exclusion_q Numeric (0-1). Required if `selection_strategy` is 
#'   "exclude_percentile". Represents the quantile threshold (e.g., 0.9 for top 10% exclusion).
#' @param random_strategy Character. One of "Fixed", "Binomial", or "Poisson". 
#'   Determines how the number of viable embryos is simulated.
#' @param p_lb Numeric (0-1). The probability of an embryo being live-born/viable.
#'   Used only if `random_strategy` is "Binomial" or "Poisson".
#'   
#' @return A named vector containing:
#'   \item{baseline}{Risk of a random embryo.}
#'   \item{selection}{Risk of the selected embryo.}
#'   \item{rel_red}{Relative risk reduction.}
#'   \item{abs_red}{Absolute risk reduction.}
#'   \item{sd_of_estimate}{Standard error of the selection estimate.}
#'   \item{sd_of_baseline}{Standard error of the baseline estimate.}
#'   
#' @importFrom matrixStats rowMins
#' @export
risk_prediction_exact <- function(iter, n, r2, h2, K,
                                  history = list(),
                                  prs_data = list(),
                                  selection_strategy = c("lowest_prs", "exclude_percentile"),
                                  exclusion_q = NULL,
                                  random_strategy = c("Fixed", "Binomial", "Poisson"),
                                  p_lb = 0) {
  # Input validations
  selection_strategy <- match.arg(selection_strategy)
  random_strategy <- match.arg(random_strategy)
  if (iter < 1) stop("'iter' must be a positive integer.")
  if (n < 1) stop("'n' must be a positive integer.")

  if (r2 < 0 || r2 > 1) stop("'r2' must be between 0 and 1.")
  if (h2 < 0 || h2 > 1) stop("'h2' must be between 0 and 1.")
  if (r2 > h2) stop("'r2' cannot be greater than 'h2'.")
  if (K <= 0 || K >= 1) stop("'K' must be strictly between 0 and 1.")
  
  if (selection_strategy == "exclude_percentile") {
    if (is.null(exclusion_q)) {
      stop("For 'exclude_percentile', 'exclusion_q' must be provided.")
    }
    if (exclusion_q <= 0 || exclusion_q >= 1) {
      stop("'exclusion_q' must be between 0 and 1.")
    }
  }
  
  if (random_strategy != "Fixed") {
    if (p_lb < 0 || p_lb > 1) stop("'p_lb' must be between 0 and 1.")
  }
  
  valid_relatives <- c("p1", "p2", "gp1a", "gp1b", "gp2a", "gp2b",
                       "sib_p1", "sib_p2", "sib_self")
  
  if (!all(names(history) %in% valid_relatives)) {
    invalid <- setdiff(names(history), valid_relatives)
    stop("Unknown names in 'history': ", paste(invalid, collapse = ", "))
  }
  if (!all(names(prs_data) %in% valid_relatives)) {
    invalid <- setdiff(names(prs_data), valid_relatives)
    stop("Unknown names in 'prs_data': ", paste(invalid, collapse = ", "))
  }
  
  # Validate history and prs
  hist_vals <- unlist(history)
  if (!is.null(hist_vals) && !all(hist_vals %in% c(-1, 0, 1, NA))) {
    stop("Values in 'history' must be: 1 (Sick), -1 (Healthy), 0 (Unknown), or NA.")
  }
  
  prs_vals <- unlist(prs_data)
  if (!is.null(prs_vals) && !is.numeric(prs_vals)) {
    stop("Values in 'prs_data' must be numeric Z-scores.")
  }
  
  # Start of calculations.
  zk <- qnorm(1 - K)
  
  # Variances for the 4 grandparents and the rest
  var_founder <- c(r2, h2 - r2, 1 - h2)
  var_seg     <- c(r2/2, (h2 - r2)/2, 1 - h2)
  # And the variance matrix for the default which is
  # the grandparents+parents
  diag_sigma  <- c(rep(var_founder, 4), rep(var_seg, 2))
  
  # Need to pad the lists in case there are siblings/cousins
  # so the indexes would be consistent.
  sanitize_lists <- function(h, p, groups) {
    for(g in groups) {
      n <- max(length(h[[g]]), length(p[[g]]))
      if (n > 0) {
        diag_sigma <<- c(diag_sigma, rep(var_seg, n))
        if(length(h[[g]]) < n) h[[g]] <- c(h[[g]], rep(0, n - length(h[[g]])))
        if(length(p[[g]]) < n) p[[g]] <- c(p[[g]], rep(NA, n - length(p[[g]])))
      }
    }
    return(list(h=h, p=p))
  }
  
  clean <- sanitize_lists(history, prs_data, 
                          c("sib_p1", "sib_p2", "sib_self"))
  history <- clean$h
  prs_data <- clean$p
  
  # Construct the equality and inequality constraints matrices
  p <- length(diag_sigma)
  eq <- construct_G_eq(prs_data, p)
  G_eq <- eq$G_eq
  vals <- eq$vals
  
  inq <- construct_G_in(history, p, zk)
  G_in <- inq$G_in
  lower_in <- inq$lower_in
  upper_in <- inq$upper_in
  
  Sigma <- diag(diag_sigma)
  
  # Start with the equality constraints
  # which is basically X | G_eq X = vals
  # using standard results on multivariate normal
  mu_z <- rep(0, p)
  Sigma_curr <- Sigma
  
  if (nrow(G_eq) > 0) {
    Sigma_G_eq <- Sigma %*% t(G_eq)
    
    # New mean
    vals2 <- solve(G_eq %*% Sigma_G_eq, vals)
    mu_z <- as.numeric(Sigma_G_eq %*% vals2)
    
    # New Sigma
    # A <- solve(G_eq %*% Sigma_G_eq, G_eq %*% t(Sigma))
    A <- solve(G_eq %*% Sigma_G_eq, t(Sigma_G_eq))
    Sigma_curr <- Sigma - Sigma_G_eq %*% A
  }
  
  # Then apply the inequality constraints
  if (nrow(G_in) > 0) {
    # Shift thresholds based on the new mean from PRS
    shift <- as.vector(G_in %*% mu_z)
    
    lower_adj <- lower_in - shift
    upper_adj <- upper_in - shift

    # Samples given the inequality
    samples_resid <- sample_func_optimized(iter, G_in, 
                                           Sigma_curr, 
                                           lower_adj, 
                                           upper_adj)
    
    # Add the mean back in case there are equality constraints
    samples <- t(t(samples_resid) + mu_z)
    
  } else {
    # No history, just PRS constraints (or nothing)
    # Sample from the conditional normal
    # chol_sigma <- chol(Sigma_curr, pivot=T)
    # samples <- mvnfast::rmvn(iter, mu_z, 
    #                          sigma = chol_sigma + 1e-8,
    #                          isChol = T)
    samples <- MASS::mvrnorm(iter, mu_z, 
                             Sigma_curr)
  }
  
  
  # Next calculate the average parents' score and non-score
  # component, as it is used for the embryos.
  s_cols <- samples[, c(1,4,7,10, 13,16)]
  w_cols <- samples[, c(2,5,8,11, 14,17)]
  weights <- c(rep(0.25, 4), 0.5, 0.5)
  
  s_p_mean <- as.vector(s_cols %*% weights)
  w_p_mean <- as.vector(w_cols %*% weights)
  
  # Baseline
  liability_mean_baseline <- s_p_mean + w_p_mean
  sd_baseline <- sqrt(1 - h2/2) 
  baseline_risks <- pnorm(zk, mean = liability_mean_baseline, 
                          sd = sd_baseline, lower.tail = FALSE)
  baseline <- mean(baseline_risks)
  
  # Selection
  if (random_strategy == "Binomial") {
    ns <- sample_truncated_binomial(iter, n, p_lb)
    max_n <- max(ns)
    noise <- matrix(rnorm(iter * max_n, sd=sqrt(r2/2)), nrow=iter)
    embryo_s_mendelian <- matrix(NA, iter, max_n)
    for(i in 1:iter) if(ns[i]>0) embryo_s_mendelian[i, 1:ns[i]] <- noise[i, 1:ns[i]]
  } else if (random_strategy == "Poisson") {
    ns <- sample_truncated_poisson(iter, n * p_lb)
    max_n <- max(ns)
    noise <- matrix(rnorm(iter * max_n, sd=sqrt(r2/2)), nrow=iter)
    embryo_s_mendelian <- matrix(NA, iter, max_n)
    for(i in 1:iter) if(ns[i]>0) embryo_s_mendelian[i, 1:ns[i]] <- noise[i, 1:ns[i]]
  } else {
    embryo_s_mendelian <- matrix(rnorm(iter * n, sd = sqrt(r2/2)), nrow = iter, ncol = n)
  }
  
  embryo_scores <- s_p_mean + embryo_s_mendelian
  
  if (selection_strategy == "lowest_prs") {
    selected_score <- matrixStats::rowMins(embryo_scores, na.rm = T)
  } else if (selection_strategy == "exclude_percentile") {
    prs_threshold <- qnorm(exclusion_q, mean = 0, sd = sqrt(r2))
    selected_score <- numeric(iter)
    for(i in 1:iter) {
      row <- embryo_scores[i, ]
      row <- row[!is.na(row)]
      valid <- row[row < prs_threshold]
      if(length(valid) > 0) selected_score[i] <- sample(valid, 1)
      else selected_score[i] <- sample(row, 1)
    }
  }
  
  liability_mean_selection <- selected_score + w_p_mean
  sd_selection <- sqrt((h2 - r2)/2 + (1 - h2))
  selection_risks <- pnorm(zk, mean = liability_mean_selection, sd = sd_selection, lower.tail = FALSE)
  selection <- mean(selection_risks)
  
  return(c(baseline = baseline, selection = selection, 
           relative_reduction = (baseline - selection)/baseline,
           absolute_reduction = baseline - selection,
           sd_of_estimate = sd(selection_risks) / sqrt(iter),
           sd_of_baseline = sd(baseline_risks) / sqrt(iter)))
}
