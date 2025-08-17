# Add random number of embryos version for the rest

library(MASS)
library(mvnfast)
library(OwenQ)

# Some silly micro-optimization instead of using one more general log(pnorm(b) - pnorm(a)), 
# should also be more accurate though.
# log_dtruncnorm_left <- function(x, a) {
#   dnorm(x, log = T) - pnorm(a, lower.tail = F, log.p = T)
# }
# 
# log_dtruncnorm_right <- function(x, b) {
#   dnorm(x, log = T) - pnorm(b, log.p = T)
# }

# The next four functions are basically the same ones as before.

risk_reduction_lowest = function(r2,K,n)
{
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

risk_reduction_exclude = function(r2,K,q,n)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zq = qnorm(q, lower.tail=F)
  integrand_t = function(t,u)
  {
    y = dnorm(t)*pnorm((zk-r/sqrt(2)*(u+t))/sqrt(1-r2),lower.tail=F)
    return(y)
  }
  # Replacement for integrand_t, should be numerically better.
  inner_integral <- function(x, a, b) {
    if (is.infinite(x)) {
      OwenT(x, a/(x * sqrt(1+b^2))) + OwenT(a/sqrt(1+b^2), (x * sqrt(1+b^2)) / a) -
        # OwenT(x, b) - OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
        OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
    }
    else {
      OwenT(x, a/(x * sqrt(1+b^2))) + OwenT(a/sqrt(1+b^2), (x * sqrt(1+b^2)) / a) -
        OwenT(x, (a+b*x)/x) - OwenT(a/sqrt(1+b^2), (a*b+x*(1+b^2))/a) + pnorm(x) * pnorm(a / sqrt(1+b^2))
    }
  }
  
  integrand_u = function(us)
  {
    y = numeric(length(us))
    beta_vec = zq*sqrt(2)-us
    denom = pnorm(beta_vec)
    dnorm_u <- dnorm(us)
    denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
    
    a <- (zk - r/sqrt(2) * us) / sqrt(1-r2)
    b <- -(r/sqrt(2)) / sqrt(1-r2)
    
    for (i in seq_along(us))
    {
      u = us[i]
      beta <- beta_vec[i]
      # internal_int = integrate(integrand_t,-Inf,beta,u, abs.tol = 1e-10)$value
      
      internal_int = ifelse(r2 != 1, pnorm(beta) - inner_integral(beta, a[i], b) + inner_integral(-Inf, a[i], b),
                            ifelse(beta < sqrt(2) / r * zk - u, 0, pnorm(beta) - pnorm(sqrt(2) / r * zk - u)))
      numer = dnorm_u[i]*(1-(1-denom[i])^n) * internal_int
      term1 = numer/denom[i]
      
      # Probably doesn't matter, but also numerically more stable
      # term1 <- exp(dnorm(u, log = T) + VGAM::log1mexp(-n*pnorm(beta, lower.tail = F, log.p = T)) - 
      #                pnorm(beta, log.p = T)) * internal_int
      
      # internal_int = integrate(integrand_t,beta,Inf,u, abs.tol = 1e-10)$value
      internal_int <- ifelse(r2 != 1, 1 - pnorm(beta) - inner_integral(Inf, a[i], b) + inner_integral(beta, a[i], b),
                             ifelse(beta < sqrt(2) / r * zk - u, 1-pnorm(sqrt(2) / r * zk - u), 1-pnorm(beta)))
      term2 = dnorm_u[i]*(1-denom[i])^(n-1) * internal_int
      y[i] = term1 + term2
      # print(sprintf("%f", y[i]))
    }
    return(y)
  }
  risk = integrate(integrand_u,-Inf,Inf, abs.tol = K * 1e-4)$value
  reduction = (K-risk)/K
  # K-risk
  return(reduction)
}

risk_reduction_lowest_conditional = function(r2,K,n,qf,qm,relative=T,parental_avg_given=F)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  if (parental_avg_given)
  {
    # It is assumed that what is given is directly the parental average, so that the paternal and maternal quantiles are the same (both representing the quantile of the parental average)
    c = zqf * r/sqrt(2)
  } else {
    c = (zqf+zqm)/2 * r
  }
  baseline = pnorm((zk-c)/sqrt(1-r2/2),lower.tail=F)
  integrand_lowest_cond = function(t)
  {
    arg = (zk-c-t*sqrt(1-r2)) / (r/sqrt(2))
    y = dnorm(t)*pnorm(arg, lower.tail=F)^n
    return(y)
  }
  # risk = integrate(integrand_lowest_cond,-Inf,Inf,rel.tol=1e-9)$value
  risk = integrate(integrand_lowest_cond,-Inf,Inf, abs.tol = baseline * 1e-4)$value
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

risk_reduction_exclude_conditional = function(r2,K,q,n,qf,qm,relative=T)
{
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
  
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

# The next two are using monte carlo simulation

baseline_risk <- function(r2, h2, K, df, dm) {
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail=F)
  
  posterior = function(gm,gf)
  {
    y = 1
    y = y * dnorm(gm/h)/h
    y = y * dnorm(gf/h)/h
    arg = (zk-gm)/sqrt(1-h2)
    if (dm)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    arg = (zk-gf)/sqrt(1-h2)
    if (df)
    {
      y = y * pnorm(arg,lower.tail=F) / K
    } else {
      y = y * pnorm(arg) / (1-K)
    }
    return(y)
  }
  
  integrand_gm = function(gms,gf)
  {
    y = numeric(length(gms))
    for (i in seq_along(gms))
    {
      gm = gms[i]
      arg = (zk - (gm+gf)/2) / sqrt(1-h2/2)
      y[i] = pnorm(arg, lower.tail=F)
      post = posterior(gm,gf)
      y[i] = y[i] * post
    }
    return(y)
  }
  
  integrand_gf = function(gfs)
  {
    y = numeric(length(gfs))
    for (i in seq_along(gfs))
    {
      gf = gfs[i]
      y[i] = integrate(integrand_gm,-Inf,Inf,gf)$value
    }
    return(y)
  }
  
  risk_baseline = integrate(integrand_gf,-Inf,Inf)$value
  return(risk_baseline)
}

risk_reduction_lowest_family_history = function(r2, h2, K, n, df, dm, n_samples = 10000)
{
  r = sqrt(r2)
  h = sqrt(h2)
  zk = qnorm(K, lower.tail = F)
  
  integrand_lowest_given_parents = function(t,gm,gf)
  {
    arg = (zk - (gm+gf)/2 - t*r/sqrt(2)) / sqrt(1-h2/2-r2/2)
    y = log(n) + dnorm(t, log = T) + (n-1) * pnorm(t, lower.tail = F, log.p = T) + pnorm(arg, lower.tail = F, log.p = T)
    return(y)
  }
  
  posterior = function(gm,gf)
  {
    y = dnorm(gm, sd = h, log = T)
    y = y + dnorm(gf, sd = h, log = T)
    arg = (zk-gm)/sqrt(1-h2)
    if (dm)
    {
      y <- y + pnorm(arg, lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    arg = (zk-gf)/sqrt(1-h2)
    if (df)
    {
      y <- y + pnorm(arg, lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    return(y)
  }
  
  integrand_gm = function(t, gms,gfs)
  {
    y <- integrand_lowest_given_parents(t, gms, gfs)
    post <- posterior(gms,gfs)
    y <- y + post
    return(y)
  }  
  
  opt <- optim(c(1, 1, 1), function(theta) -integrand_gm(theta[1], theta[2], theta[3]), hessian = T)
  var_matrix <- solve(opt$hessian)
  means <- opt$par
  
  # Sample points from t distribution
  data <- rmvt(n_samples, sigma = var_matrix, df = 5, mu = means)
  y <- exp(integrand_gm(data[, 1], data[, 2], data[, 3]))

  risk_selection = mean(y / dmvt(data, var_matrix, df = 5, mu = means, log = F))
  sd <- sqrt(mean((y / dmvt(data, var_matrix, df = 5, mu = means, log = F) - risk_selection)^2) / n_samples)
  risk_baseline <- baseline_risk(r2, h2, K, df, dm)
  
  relative_reduction = (risk_baseline-risk_selection)/risk_baseline
  abs_reduction = risk_baseline-risk_selection
  return(c(risk_baseline,risk_selection,relative_reduction,abs_reduction, sd))
}

risk_reduction_exclude_family_history <- function(r2, h2, K, q, n, df, dm, n_samples = 10000)
{
  # We need h > r, so if not, subtract epsilon from r. Unless h is zero, and then we add epsilon to h
  if (h2 == 0) {h2 <- h2 + 0.0001}
  else if (h2 == r2) {r2 <- r2 - 0.0001}
  
  r <- sqrt(r2)
  h <- sqrt(h2)
  zk <- qnorm(K, lower.tail = F)
  zq <- qnorm(q, lower.tail = F)
  
  posterior <- function(gm,gf)
  {
    y <- dnorm(gm, sd = h, log = T)
    y <- y + dnorm(gf, sd = h, log = T)
    arg <- (zk-gm)/sqrt(1-h2)
    if (dm)
    {
      y <- y + pnorm(arg,lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    arg <- (zk-gf)/sqrt(1-h2)
    if (df)
    {
      y <- y + pnorm(arg,lower.tail = F, log.p = T) - log(K)
    } else {
      y <- y + pnorm(arg, log.p = T) - log(1-K)
    }
    return(y)
  }
  
  integrand_t <- function(t,gm,gf)
  {
    arg <- (zk-t*r/sqrt(2)-(gm+gf)/2)/sqrt(1-h2/2-r2/2)
    y <- dnorm(t, log = T) + pnorm(arg, lower.tail = F, log.p = T)
    return(y)
  }
  
  integrand_c1 <- function(cs, gm, gf, t)
  {
    gamma = zq*sqrt(2) - cs/(r/sqrt(2))
    denom <- pnorm(gamma)
    denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
    f <- integrand_t(t, gm, gf)
    
    numer <- log(1-(1-denom)^n)
    y <- numer-log(denom) + f
    y <- y + dnorm(cs, mean=r2/h2 * (gm+gf)/2, sd=r/h * sqrt((h2-r2)/2), log = T)
    
    return(y)
  }
  
  integrand_c2 <- function(cs, gm, gf, t)
  {
    gamma = zq*sqrt(2) - cs/(r/sqrt(2))
    denom <- pnorm(gamma)
    denom <- ifelse(denom==0, 1e-300, denom) # Avoid dividing by zero
    f <- integrand_t(t, gm, gf)
    
    numer <- log(1-(1-denom)^n)
    y <- (n-1) * log(1-denom) + f
    y <- y + dnorm(cs, mean=r2/h2 * (gm+gf)/2, sd=r/h * sqrt((h2-r2)/2), log = T)
    
    return(y)
  }
  
  integrand_gm1 = function(gms,gfs, c, t)
  {
    y <- integrand_c1(c, gms, gfs, t)
    post <- posterior(gms, gfs)
    y <- y + post
    return(y)
  }
  
  integrand_gm2 = function(gms,gfs, c, t)
  {
    y <- integrand_c2(c, gms, gfs, t)
    post <- posterior(gms, gfs)
    y <- y + post
    return(y)
  }
  # df didn't really matter, but I chose higher one since otherwise we get a bit more NA (due to pnorm(gamma) being 1, and then
  # the truncated distribution might sample infinities)
  degrees_of_freedom <- 5
  
  opt1 <- optim(c(1, 1, zq*r, 1), function(theta) -integrand_gm1(theta[1], theta[2], theta[3], theta[4]), hessian = T)
  var_matrix1 <- solve(opt1$hessian)
  means1 <- opt1$par
  
  data1 <- rmvt(n_samples, sigma = var_matrix1[1:3, 1:3], df = degrees_of_freedom, mu = means1[1:3])
  data1 <- cbind(data1, qnorm(runif(n_samples) * pnorm(zq*sqrt(2) - data1[, 3]/(r/sqrt(2)), means1[4], sqrt(var_matrix1[4, 4])),
                              means1[4], sqrt(var_matrix1[4, 4])))
  y1 <- exp(integrand_gm1(data1[, 1], data1[, 2], data1[, 3], data1[, 4])) / 
    exp(dmvt(data1[, 1:3], var_matrix1[1:3, 1:3], df = degrees_of_freedom, mu = means1[1:3], log = T) + 
      dnorm(data1[, 4], means1[4], sqrt(var_matrix1[4, 4]), log = T) - 
      pnorm(zq*sqrt(2) - data1[, 3]/(r/sqrt(2)), means1[4], sqrt(var_matrix1[4, 4]), log.p = T))
  
  opt2 <- optim(c(1, 1, zq*r, 1), function(theta) -integrand_gm2(theta[1], theta[2], theta[3], theta[4]), hessian = T)
  var_matrix2 <- solve(opt2$hessian)
  means2 <- opt2$par
  
  data2 <- rmvt(n_samples, sigma = var_matrix2[1:3, 1:3], df = degrees_of_freedom, mu = means2[1:3])
  pnorm_temp <- pnorm(zq*sqrt(2) - data2[, 3]/(r/sqrt(2)), means2[4], sqrt(var_matrix2[4, 4]))
  data2 <- cbind(data2, qnorm(pnorm_temp+runif(n_samples)*(1-pnorm_temp), 
                              opt2$par[4], sqrt(var_matrix2[4, 4])))
  y2 <- exp(integrand_gm2(data2[, 1], data2[, 2], data2[, 3], data2[, 4])) /
    exp(dmvt(data2[, 1:3], var_matrix2[1:3, 1:3], df = degrees_of_freedom, mu = means2[1:3], log = T) +
       dnorm(data2[, 4], means2[4], sqrt(var_matrix2[4, 4]), log = T) - 
       pnorm(zq*sqrt(2) - data2[, 3]/(r/sqrt(2)), means2[4], sqrt(var_matrix2[4, 4]), lower.tail = F, log.p = T))

  risk_selection = mean(y1, na.rm = T) + mean(y2, na.rm = T)
  sd <- sqrt(((mean((y1 - mean(y1, na.rm=T))^2, na.rm = T) + (mean((y2 - mean(y2, na.rm = T))^2, na.rm = T))) / n_samples))
  risk_baseline <- baseline_risk(r2, h2, K, df, dm)
  
  relative_reduction = (risk_baseline-risk_selection)/risk_baseline
  abs_reduction = risk_baseline-risk_selection
  return(c(risk_baseline,risk_selection,relative_reduction,abs_reduction, sd))
}

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

sample_func <- function(n, zk, G, Sigma, directions) {
  # Need to extend it so different variable sum should be zk.
  d <- nrow(Sigma)
  if (ncol(G) != d) {stop("The number of cols in G should be the same as the dimension of Sigma.")}
  if (length(directions) != length(zk)) {stop("The vector directions should be the same length as the vector zk.")}
  if (!all(directions %in% c(">", "<"))) {stop("directions should be one of \">\", or \"<\".")}
  
  chol_sigma <- chol(Sigma)
  G_chol_sigma <- tcrossprod(G %*% chol_sigma)
  
  # temp <- Sigma %*% t(G) %*% MASS::ginv(G %*% Sigma %*% t(G))
  temp <- Sigma %*% t(G) %*% MASS::ginv(G_chol_sigma)
  # print(temp)
  # r <- truncnorm::rtruncnorm(n, zk, sd = sqrt(sum(Sigma)))
  
  # Maybe need 
  # temp_G_sigma <- rowSums(G %*% Sigma)
  
  # ARGH, the problem is that when the sums are not
  # independent, the truncated normals here are not either.
  # r <- t(sapply(1:length(zk), function(i) {
  #   # Note: won't work with non block structure.
  #   if (directions[i] == ">") {
  #     # truncnorm::rtruncnorm(n, a=zk[i], sd = sqrt(sum(Sigma)))
  #     truncnorm::rtruncnorm(n, a=zk[i], sd = sqrt(temp_G_sigma[i]))
  #   }
  #   else {
  #     truncnorm::rtruncnorm(n, b=zk[i], sd = sqrt(temp_G_sigma[i]))
  #   }
  # }))
  # print(dim(r))
  # Ok, works. But really slow.
  # r <- tmvtnorm::rtmvnorm(n, sigma = G %*% Sigma %*% t(G), 
  #                         lower = ifelse(directions == ">", zk, -Inf),
  #                         upper = ifelse(directions == "<", zk, Inf),
  #                         algorithm = "gibbs")
  r <- tmvtnorm::rtmvnorm(n, sigma = G_chol_sigma, 
                          lower = ifelse(directions == ">", zk, -Inf),
                          upper = ifelse(directions == "<", zk, Inf),
                          algorithm = "gibbs")
  
  # print(summary(r))
  r <- t(r)
  # print(dim(r))
  
  # y <- mvnfast::rmvn(n, rep(0, d), sigma = Sigma)
  y <- mvnfast::rmvn(n, rep(0, d), sigma = chol_sigma, isChol = T)
  
  # y + matrix((r - as.numeric(tcrossprod(G, y))), nrow=n, ncol=d) / d
  y + matrix(temp %*% (r-tcrossprod(G, y)), nrow=n, ncol=d, byrow=T)
  # t(sapply(1:n, function(i) y[i, ] + Sigma %*% t(G) %*% MASS::ginv(G %*% Sigma %*% t(G)) %*% (r[i]-G %*% y[i, ])))
}

risk_parents_offspring_exclude <- function(iter, n, r2, h2, K,
                                           q,
                                           sick_parents,
                                           no_sick_parents,
                                           sick_siblings,
                                           no_sick_siblings,
                                           qf=NULL,
                                           qm=NULL) {
  # Isn't exact at the moment. Argh
  # I think that's fixed ^
  # OK so here is the thing - when comparing the estimation from
  # this function to risk_calc_parent (and here we condition)
  # on one sick parent - we get the same result. 
  # but when comapring it to the older risk_reduction_lowest_family_history
  # with one sick parent, we get a slightly different result.
  
  # For two parents, we get the same thing
  # and also for when sibling when comparing to risk_calc
  if(any(c(sick_parents,
           no_sick_parents,
           sick_siblings,
           no_sick_siblings) < 0)) {
    stop("")
  }
  if (sick_parents + no_sick_parents > 2) {
    stop("You can only condition on 2 parents.")
  }
  
  zk <- qnorm(1-K)
  Sigma <- NULL
  samples <- NULL
  G <- NULL
  
  # For parents, we have
  # s_f, w_f, e_f
  # s_m, w_m, e_m
  # With no covariance and variance r2, h2-r2, 1-h2
  
  # For offspring we have
  # s_c = (s_f + s_m)/2
  # w_c = (w_f + w_m)/2
  # s_i, w_i, e_i
  # with covariance due to s_c, w_c
  # and variance r2/2, (h2-r2)/2, 1-h2/2?
  # Let's start by supporting only 0 parents or 2 parents. In this case:
  if (is.null(qf)) {
  if (sick_parents + no_sick_parents == 2) {
    Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
    G <- rbind(c(1, 1, 1, rep(0, 3 + sick_siblings+no_sick_siblings)),
               c(0, 0, 0, 1, 1, 1, rep(0, sick_siblings+no_sick_siblings)),
               cbind(matrix(c(0.5, 0.5, 0, 0.5, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
    
    samples <- sample_func(iter, rep(zk, 2+sick_siblings+no_sick_siblings), G, 
                           Sigma, c(rep(">", sick_parents),
                                    rep("<", no_sick_parents),
                                    rep(">", sick_siblings),
                                    rep("<", no_sick_siblings)))
  }
  else if (sick_parents + no_sick_parents == 1) {
    Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
    G <- rbind(c(1, 1, 1, rep(0, 3 + sick_siblings+no_sick_siblings)),
               cbind(matrix(c(0.5, 0.5, 0, 0.5, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
    
    samples <- sample_func(iter, rep(zk, 1+sick_siblings+no_sick_siblings), G, 
                           Sigma, c(rep(">", sick_parents),
                                    rep("<", no_sick_parents),
                                    rep(">", sick_siblings),
                                    rep("<", no_sick_siblings)))
  }
  else if (sick_parents + no_sick_parents == 0) {
    # So generate only siblings.
    Sigma <- diag(c(r2/2, (h2-r2)/2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
    G <- cbind(matrix(1, nrow=sick_siblings+no_sick_siblings, ncol=2),
               diag(sick_siblings+no_sick_siblings))
    
    samples <- sample_func(iter, rep(zk, sick_siblings+no_sick_siblings), G, 
                           Sigma, c(rep(">", sick_siblings),
                                    rep("<", no_sick_siblings)))
  }
  }
  else {
    zqf <- sqrt(r2) * qf
    zqm <- sqrt(r2) * qm
    if (sick_parents + no_sick_parents == 2) {
      Sigma <- diag(c(h2-r2, 1-h2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
      G <- rbind(c(1, 1, rep(0, 2+sick_siblings+no_sick_siblings)),
                 c(0, 0, 1, 1, rep(0, sick_siblings+no_sick_siblings)),
                 cbind(matrix(c(0.5, 0, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
      
      samples <- sample_func(iter, c(zk-zqf, zk-zqm,
                                     rep(zk-(zqf+zqm)/2, sick_siblings+no_sick_siblings)), G, 
                             Sigma, c(rep(">", sick_parents),
                                      rep("<", no_sick_parents),
                                      rep(">", sick_siblings),
                                      rep("<", no_sick_siblings)))
    }
    else if (sick_parents + no_sick_parents == 1) {
      # Didn't change that
      Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
      G <- rbind(c(1, 1, 1, rep(0, 3 + sick_siblings+no_sick_siblings)),
                 cbind(matrix(c(0.5, 0.5, 0, 0.5, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
      
      samples <- sample_func(iter, rep(zk, 1+sick_siblings+no_sick_siblings), G, 
                             Sigma, c(rep(">", sick_parents),
                                      rep("<", no_sick_parents),
                                      rep(">", sick_siblings),
                                      rep("<", no_sick_siblings)))
    }
    else if (sick_parents + no_sick_parents == 0) {
      # So generate only siblings.
      Sigma <- diag(c(r2/2, (h2-r2)/2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
      G <- cbind(matrix(1, nrow=sick_siblings+no_sick_siblings, ncol=2),
                 diag(sick_siblings+no_sick_siblings))
      
      samples <- sample_func(iter, rep(zk-(zqf+zqm)/2, sick_siblings+no_sick_siblings), G, 
                             Sigma, c(rep(">", sick_siblings),
                                      rep("<", no_sick_siblings)))
    }
  }
  
  
  # Two parents
  # Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2))
  # samples <- sample_func(iter, zk, Sigma)
  # Sigma <- diag(c(r2, h2-r2, 1-h2))
  # samples <- cbind(sample_func2(iter, zk, Sigma),
  # sample_func2(iter, zk, Sigma))
  # samples <- sample_func(iter, c(zk, zk), rbind(c(1, 1, 1, 0, 0, 0),
  # c(0, 0, 0, 1, 1, 1)),
  # Sigma, c(">", ">"))
  
  # Two siblings
  # Sigma <- diag(c(r2/2, (h2-r2)/2, 1-h2/2, 1-h2/2))
  # G <- rbind(c(1, 1, 1, 0),
  #            c(1, 1, 0, 1))
  
  # sample_r_mean <- rep(dnorm(zk) / pnorm(zk, lower.tail=F), 2)
  # samples_means <- as.numeric(Sigma %*% t(G) %*% MASS::ginv(G %*% Sigma %*% t(G)) %*% sample_r_mean)
  
  # samples <- sample_func(iter, c(zk, zk), G, 
  #                        Sigma, c(">", ">"))
  
  scores <- matrix(rnorm(nrow(samples)*n, sd = sqrt(r2/2)), 
                   nrow=nrow(samples),
                   ncol=n) + samples[, 1]
  mean_1 <- samples[, 1]+samples[, 2]
  
  if(sick_parents + no_sick_parents %in% c(1, 2)) {
    scores <- matrix(rnorm(nrow(samples)*n, sd = sqrt(r2/2)), 
                     nrow=nrow(samples), 
                     ncol=n) + (samples[, 1]+samples[, 4])/2
    mean_1 <- (samples[, 1]+samples[, 4])/2 + 
      (samples[, 2]+samples[, 5])/2
  }
  
  sds <- sqrt(r2/2+(h2-r2)/2+1-h2)
  
  baseline <- mean(pnorm(zk, mean_1, sds, lower.tail = F))
  
  # Select randomly out of those with prs < z_q
  q <- qnorm(1-q)
  min_score <- apply(scores, 1, function(x) {
    if (all(x > q)) {
      return(sample(x, 1))
    }
    if (length(which(x < q)) == 1) {
      return(x[which(x < q)])
    }
    return(sample(x[which(x < q)], 1))
  })
  
  mean_1 <- min_score + samples[, 2]
  
  if(sick_parents + no_sick_parents %in% c(1, 2)) {
    mean_1 <- min_score + (samples[, 2]+samples[, 5])/2
  }
  
  # P(e_min > z_k | mu = mean_1, sd = sds)
  
  sds <- sqrt((h2-r2)/2 + 1-h2)
  selection <- mean(pnorm(zk, mean_1, sds, lower.tail = F))
  
  # return(list(baseline=baseline, selection=selection))
  return(c(baseline, selection, (baseline-selection)/selection,
           baseline-selection, sd(pnorm(zk, mean_1, sds, lower.tail = F))))
}

risk_parents_offspring <- function(iter, n, r2, h2, K,
                                   sick_parents,
                                   no_sick_parents,
                                   sick_siblings,
                                   no_sick_siblings,
                                   qf=NULL, qm=NULL) {
  # I think 1 sick parent+1 non-sick parent is incorrect since I still
  # need to condition on the two parents, just the second one is with zk <...
  
  # Isn't exact at the moment. Argh
  # I think that's fixed ^
  # OK so here is the thing - when comparing the estimation from
  # this function to risk_calc_parent (and here we condition)
  # on one sick parent - we get the same result. 
  # but when comparing it to the older risk_reduction_lowest_family_history
  # with one sick parent, we get a slightly different result.
  
  # For two parents, we get the same thing
  # and also for when sibling when comparing to risk_calc
  if(any(c(sick_parents,
           no_sick_parents,
           sick_siblings,
           no_sick_siblings) < 0)) {
    stop("")
  }
  if (sick_parents + no_sick_parents > 2) {
    stop("You can only condition on 2 parents.")
  }
  
  zk <- qnorm(1-K)
  Sigma <- NULL
  samples <- NULL
  G <- NULL
  
  # For parents, we have
  # s_f, w_f, e_f
  # s_m, w_m, e_m
  # With no covariance and variance r2, h2-r2, 1-h2
  
  # For offspring we have
  # s_c = (s_f + s_m)/2
  # w_c = (w_f + w_m)/2
  # s_i, w_i, e_i
  # with covariance due to s_c, w_c
  # and variance r2/2, (h2-r2)/2, 1-h2/2?
  # Let's start by supporting only 0 parents or 2 parents. In this case:
  if (is.null(qf)) {
  if (sick_parents + no_sick_parents == 2) {
    Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
    G <- rbind(c(1, 1, 1, rep(0, 3 + sick_siblings+no_sick_siblings)),
               c(0, 0, 0, 1, 1, 1, rep(0, sick_siblings+no_sick_siblings)),
               cbind(matrix(c(0.5, 0.5, 0, 0.5, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
    
    samples <- sample_func(iter, rep(zk, 2+sick_siblings+no_sick_siblings), G, 
                           Sigma, c(rep(">", sick_parents),
                                    rep("<", no_sick_parents),
                                    rep(">", sick_siblings),
                                    rep("<", no_sick_siblings)))
  }
  else if (sick_parents + no_sick_parents == 1) {
    Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
    G <- rbind(c(1, 1, 1, rep(0, 3 + sick_siblings+no_sick_siblings)),
               cbind(matrix(c(0.5, 0.5, 0, 0.5, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
    
    samples <- sample_func(iter, rep(zk, 1+sick_siblings+no_sick_siblings), G, 
                           Sigma, c(rep(">", sick_parents),
                                    rep("<", no_sick_parents),
                                    rep(">", sick_siblings),
                                    rep("<", no_sick_siblings)))
  }
  else if (sick_parents + no_sick_parents == 0) {
    # So generate only siblings.
    Sigma <- diag(c(r2/2, (h2-r2)/2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
    G <- cbind(matrix(1, nrow=sick_siblings+no_sick_siblings, ncol=2),
               diag(sick_siblings+no_sick_siblings))
    
    samples <- sample_func(iter, rep(zk, sick_siblings+no_sick_siblings), G, 
                           Sigma, c(rep(">", sick_siblings),
                                    rep("<", no_sick_siblings)))
  }
  }
  else {
    zqf <- sqrt(r2) * qf
    zqm <- sqrt(r2) * qm
    if (sick_parents + no_sick_parents == 2) {
      Sigma <- diag(c(h2-r2, 1-h2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
      
      G <- rbind(c(1, 1, rep(0, 2+sick_siblings+no_sick_siblings)),
                 c(0, 0, 1, 1, rep(0, sick_siblings+no_sick_siblings)),
                 cbind(matrix(c(0.5, 0, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
      
      samples <- sample_func(iter, c(zk-zqf, zk-zqm,
                                     rep(zk-(zqf+zqm)/2, sick_siblings+no_sick_siblings)), G, 
                             Sigma, c(rep(">", sick_parents),
                                      rep("<", no_sick_parents),
                                      rep(">", sick_siblings),
                                      rep("<", no_sick_siblings)))
    }
    else if (sick_parents + no_sick_parents == 1) {
      # Didn't change that
      Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
      G <- rbind(c(1, 1, 1, rep(0, 3 + sick_siblings+no_sick_siblings)),
                 cbind(matrix(c(0.5, 0.5, 0, 0.5, 0.5, 0), ncol=6, nrow=sick_siblings+no_sick_siblings, byrow=T), diag(sick_siblings+no_sick_siblings)))
      
      samples <- sample_func(iter, rep(zk, 1+sick_siblings+no_sick_siblings), G, 
                             Sigma, c(rep(">", sick_parents),
                                      rep("<", no_sick_parents),
                                      rep(">", sick_siblings),
                                      rep("<", no_sick_siblings)))
    }
    else if (sick_parents + no_sick_parents == 0) {
      # So generate only siblings.
      Sigma <- diag(c(r2/2, (h2-r2)/2, rep(1-h2/2, sick_siblings+no_sick_siblings)))
      G <- cbind(matrix(1, nrow=sick_siblings+no_sick_siblings, ncol=2),
                 diag(sick_siblings+no_sick_siblings))
      
      samples <- sample_func(iter, rep(zk-(zqf+zqm)/2, sick_siblings+no_sick_siblings), G, 
                             Sigma, c(rep(">", sick_siblings),
                                      rep("<", no_sick_siblings)))
    }
  }
  
  # Two parents
  # Sigma <- diag(c(r2, h2-r2, 1-h2, r2, h2-r2, 1-h2))
  # samples <- sample_func(iter, zk, Sigma)
  # Sigma <- diag(c(r2, h2-r2, 1-h2))
  # samples <- cbind(sample_func2(iter, zk, Sigma),
  # sample_func2(iter, zk, Sigma))
  # samples <- sample_func(iter, c(zk, zk), rbind(c(1, 1, 1, 0, 0, 0),
  # c(0, 0, 0, 1, 1, 1)),
  # Sigma, c(">", ">"))
  
  # Two siblings
  # Sigma <- diag(c(r2/2, (h2-r2)/2, 1-h2/2, 1-h2/2))
  # G <- rbind(c(1, 1, 1, 0),
  #            c(1, 1, 0, 1))
  
  # sample_r_mean <- rep(dnorm(zk) / pnorm(zk, lower.tail=F), 2)
  # samples_means <- as.numeric(Sigma %*% t(G) %*% MASS::ginv(G %*% Sigma %*% t(G)) %*% sample_r_mean)
  
  # samples <- sample_func(iter, c(zk, zk), G, 
  #                        Sigma, c(">", ">"))
  
  scores <- matrix(rnorm(nrow(samples)*n, sd = sqrt(r2/2)), 
                   nrow=nrow(samples),
                   ncol=n) + samples[, 1]
  mean_1 <- samples[, 1]+samples[, 2]
  
  if(sick_parents + no_sick_parents %in% c(1, 2)) {
    scores <- matrix(rnorm(nrow(samples)*n, sd = sqrt(r2/2)), 
                     nrow=nrow(samples), 
                     ncol=n) + (samples[, 1]+samples[, 4])/2
    mean_1 <- (samples[, 1]+samples[, 4])/2 + 
      (samples[, 2]+samples[, 5])/2
  }
  
  sds <- sqrt(r2/2+(h2-r2)/2+1-h2)
  
  baseline <- mean(pnorm(zk, mean_1, sds, lower.tail = F))
  
  min_score <- matrixStats::rowMins(scores)
  
  mean_1 <- min_score + samples[, 2]
  
  if(sick_parents + no_sick_parents %in% c(1, 2)) {
    mean_1 <- min_score + (samples[, 2]+samples[, 5])/2
  }
  
  # P(e_min > z_k | mu = mean_1, sd = sds)
  
  sds <- sqrt((h2-r2)/2 + 1-h2)
  selection <- mean(pnorm(zk, mean_1, sds, lower.tail = F))
  
  # return(list(baseline=baseline, selection=selection))
  return(c(baseline, selection, (baseline-selection)/selection,
           baseline-selection, sd(pnorm(zk, mean_1, sds, lower.tail = F))))
}

# binomial_random <- function(selection_risk, p, n) {
#   ((1-p*selection_risk) - (1-p)^n) / (1-(1-p)^n)
# }

risk_reduction_lowest_bin <- function(r2, K, n, p)
{
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

risk_reduction_lowest_pois <- function(r2, K, lambda)
{
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

risk_reduction_lowest_conditional_bin = function(r2,K,n,qf,qm,relative=T,parental_avg_given=F, p)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  if (parental_avg_given)
  {
    # It is assumed that what is given is directly the parental average, so that the paternal and maternal quantiles are the same (both representing the quantile of the parental average)
    c = zqf * r/sqrt(2)
  } else {
    c = (zqf+zqm)/2 * r
  }
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
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(list(rr=reduction, baseline=baseline, risk=risk))
}

risk_reduction_lowest_conditional_pois = function(r2,K,lambda,qf,qm,relative=T,parental_avg_given=F)
{
  r = sqrt(r2)
  zk = qnorm(K, lower.tail=F)
  zqf = qnorm(qf)
  zqm = qnorm(qm)
  if (parental_avg_given)
  {
    # It is assumed that what is given is directly the parental average, so that the paternal and maternal quantiles are the same (both representing the quantile of the parental average)
    c = zqf * r/sqrt(2)
  } else {
    c = (zqf+zqm)/2 * r
  }
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
  if (relative) {
    reduction = (baseline-risk)/baseline
  } else {
    reduction = baseline-risk
  }
  return(list(rr=reduction, baseline=baseline, risk=risk))
}





sample_truncated_binomial <- function(num_samples, n, p) {
  # Calculate all relevant PMFs at once
  pmf_vals <- dbinom(0:n, size = n, prob = p)
  
  prob_gt_zero <- 1 - pmf_vals[1] # pmf_vals[1] corresponds to k=0
  
  if (prob_gt_zero <= 1e-10) {
    warning("Probability of non-zero outcome is low. Returning 1.")
    return(rep(1, num_samples))
  }
  
  # Get the probabilities for k=1...n and normalize them
  truncated_pmf <- pmf_vals[-1] / prob_gt_zero
  
  # Calculate the CDF of the truncated distribution
  truncated_cdf <- cumsum(truncated_pmf)
  
  # Generate all random numbers at once
  u <- runif(num_samples)
  
  # For each u, find the first index in the CDF that is >= u.
  # The `findInterval` function is perfect and highly optimized for this.
  # We add 1 because findInterval returns a 0-based index into the k=1...n space.
  samples <- findInterval(u, truncated_cdf) + 1
  
  return(samples)
}

sample_truncated_poisson <- function(num_samples, lambda) {
  
  # 1. Calculate the probability of the outcome being > 0.
  prob_gt_zero <- 1 - exp(-lambda)
  
  # Handle the edge case where lambda is so small that P(X>0) is nearly zero.
  if (prob_gt_zero <= 1e-10) {
    warning("Lambda is very small; probability of non-zero outcome is negligible. All samples will be 1.")
    return(rep(1, num_samples))
  }
  
  # 2. Determine a safe upper limit for k to calculate the CDF.
  # qpois gives the quantile; we go far into the tail to ensure our CDF is ~1.
  # We add a buffer (e.g., +5) just to be extra safe.
  k_max <- qpois(1 - 1e-15, lambda) + 5 
  
  # 3. Calculate the PMF for k=1...k_max and normalize it.
  # We calculate from k=0 to reuse dpois, then discard the k=0 term.
  pmf_vals <- dpois(0:k_max, lambda)
  truncated_pmf <- pmf_vals[-1] / prob_gt_zero # pmf_vals[-1] is P(X=1), P(X=2),...
  
  # 4. Calculate the CDF of the truncated distribution.
  truncated_cdf <- cumsum(truncated_pmf)
  
  # 5. Generate all random numbers at once.
  u <- runif(num_samples)
  
  # 6. For each u, find the first index in the CDF that is >= u.
  # `findInterval` is extremely fast for this task.
  # The result is the sample value k (since our CDF is for k=1, 2, ...).
  samples <- findInterval(u, truncated_cdf) + 1
  
  return(samples)
}

sample_func <- function(n, zk, G, Sigma, directions) {
  d <- nrow(Sigma)
  if (ncol(G) != d) {stop("The number of cols in G should be the same as the dimension of Sigma.")}
  if (length(directions) != length(zk)) {stop("The vector directions should be the same length as the vector zk.")}
  if (!all(directions %in% c(">", "<"))) {stop("directions should be one of \">\", or \"<\".")}
  
  chol_sigma <- chol(Sigma)
  # G_chol_sigma <- tcrossprod(G %*% chol_sigma) # This is G %*% Sigma %*% t(G)
  G_chol_sigma <- G %*% Sigma %*% t(G)
  
  # The conditional expectation operator: Sigma * t(G) * inv(G*Sigma*t(G))
  temp <- Sigma %*% t(G) %*% solve(G_chol_sigma)
  
  # 1. Sample the truncated liabilities 'r'
  r <- tmvtnorm::rtmvnorm(n, sigma = G_chol_sigma, 
                          lower = ifelse(directions == ">", zk, -Inf),
                          upper = ifelse(directions == "<", zk, Inf),
                          algorithm = "gibbs")
  r <- t(r)
  
  # 2. Sample an unconditional z (named y here) ~ N(0, Sigma)
  y <- mvnfast::rmvn(n, rep(0, d), sigma = chol_sigma, isChol = T)
  
  # 3. Correct y to get a sample from the conditional distribution P(z | Gz=r)
  # Formula: E[z|Gz=r] = temp * r. Correction: z = y - E[y|Gy] + E[z|Gz=r]
  correction <- temp %*% (r - tcrossprod(G, y))
  
  # Return the corrected sample
  y + t(correction)
}

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
  
  # --- 1. Input Validation and Basic Setup ---
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
  
  zk <- qnorm(1 - K)
  n_known_parents <- sick_parents + no_sick_parents
  n_known_sibs <- sick_siblings + no_sick_siblings
  
  s_p_mean <- NULL
  w_p_mean <- NULL
  
  # FIX: Handle the r2 == h2 case explicitly.
  # The w component variance is 0, so it's a constant, not a random variable.
  is_w_zero <- abs(h2 - r2) < 1e-9
  
  # --- 2. Sample Parental Genetic Components Conditional on Family History ---
  
  if (is.null(qf)) {
    # --- CASE A: Parental PRS is UNKNOWN ---
    if (n_known_parents > 0) {
      # Z = (s_f, [w_f], e_f, s_m, [w_m], e_m, sib_specific_1, ...)
      # The [w] components are included only if their variance is non-zero.
      
      parent_vars <- if (is_w_zero) c(r2, 1-h2) else c(r2, h2-r2, 1-h2)
      cols_per_parent <- length(parent_vars)
      s_idx <- 1; w_idx <- if(is_w_zero) NULL else 2;
      
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
      s_idx <- 1; w_idx <- if(is_w_zero) NULL else 2;
      
      n_vars <- length(shared_vars) + n_known_sibs
      Sigma <- diag(c(shared_vars, rep(1-h2/2, n_known_sibs)))
      
      G <- cbind(matrix(1, nrow=n_known_sibs, ncol=length(shared_vars)), diag(n_known_sibs))
      
      samples <- sample_func(iter, rep(zk, n_known_sibs), G, Sigma, 
                             c(rep(">", sick_siblings), rep("<", no_sick_siblings)))
      
      s_p_mean <- samples[, s_idx]
      w_p_mean <- if (is_w_zero) 0 else samples[, w_idx]
    }
    
  } else {
    # --- CASE B: Parental PRS is KNOWN ---
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
  
  # --- 3. Calculate Risk for Baseline and Selected Embryo ---
  
  if (is.null(s_p_mean)) {
    # No family history provided at all.
    return(c(baseline=K, selection=K, rel_red=0, abs_red=0, sd=0))
  }
  
  # A. Baseline risk (random embryo)
  liability_mean_baseline <- s_p_mean + w_p_mean
  sd_baseline <- sqrt(1 - h2/2)
  baseline_risks <- pnorm(zk, mean = liability_mean_baseline, sd = sd_baseline, lower.tail = FALSE)
  baseline <- mean(baseline_risks)
  
  # B. Selection risk
  if (random_strategy == "Binomial") {
    ns <- sample_truncated_binomial(iter, n, p)
    embryo_s_mendelian <- sapply(ns, function(ns) c(rnorm(ns, sd=sqrt(r2/2)), 
                                                    rep(NA, n-ns))) |>
      t()
  }
  else if (random_strategy == "Poisson") {
    # ns <- sample_truncated_poisson(iter, n * p)
    ns <- sample_truncated_binomial(iter, n, p)
    max_ns <- max(ns)
    embryo_s_mendelian <- sapply(ns, function(ns) c(rnorm(ns, sd=sqrt(r2/2)), 
                                                    rep(NA, max_ns-ns))) |>
      t()
  }
  else {
    embryo_s_mendelian <- matrix(rnorm(iter * n, sd = sqrt(r2/2)), nrow = iter, ncol = n)
  }
  
  embryo_scores <- s_p_mean + embryo_s_mendelian
  
  # Select one embryo based on the chosen strategy
  if (selection_strategy == "lowest_prs") {
    selected_score <- matrixStats::rowMins(embryo_scores, na.rm = T)
  } else if (selection_strategy == "exclude_percentile") {
    prs_threshold <- qnorm(exclusion_q, mean = 0, sd = sqrt(r2)) # Threshold on the absolute PRS scale
    
    selected_score <- sapply(1:iter, function(i) {
      eligible_scores <- embryo_scores[i, embryo_scores[i,] < prs_threshold]
      if (length(eligible_scores) > 0) {
        # Your logic: sample randomly from the eligible embryos
        return(sample(eligible_scores, 1))
        # Alternative powerful strategy: pick the best of the eligible embryos
        # return(min(eligible_scores)) 
      } else {
        # Fallback: if no embryo meets the criteria, sample one from the full set
        # This follows the logic in your 'risk_parents_offspring_exclude' function
        return(sample(embryo_scores[i,], 1))
      }
    })
  } else {
    stop("Invalid selection_strategy provided.")
  }
  
  liability_mean_selection <- selected_score + w_p_mean
  sd_selection <- sqrt((h2 - r2)/2 + (1 - h2))
  selection_risks <- pnorm(zk, mean = liability_mean_selection, sd = sd_selection, lower.tail = FALSE)
  selection <- mean(selection_risks)
  
  # print(cov(selection_risks, baseline_risks))
  # --- 4. Return Results ---
  return(c(baseline = baseline,
           selection = selection,
           relative_reduction = (baseline - selection) / baseline,
           absolute_reduction = baseline - selection,
           sd_of_estimate = sd(selection_risks) / sqrt(iter),
           sd_of_baseline = sd(baseline_risks) / sqrt(iter)))
}

# risk_parents_offspring_generic(1e4, 5, 0.05, 0.4, 0.001, 1, 0, 0, 0,
#                                random_strategy = "Fixed")
# 
# risk_reduction_lowest_family_history(0.1, 0.4, 0.05, 5, T, F, n_samples = 1e4)
# 
# risk_parents_offspring_generic(1e4, 5, 0.1, 0.4, 0.05, 1, 0, 0, 0, 
#                                random_strategy = "Poisson", p=0.5)
