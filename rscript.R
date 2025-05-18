library(tidyverse)
library(cmdstanr)

# > str(rain)
# 'data.frame':	23203 obs. of  2 variables:
#  $ t: POSIXct, format: "1960-01-01" "1960-01-02" "1960-01-03" "1960-01-04" ...
#  $ y: num  0 1.37 0.127 0 0 ...



# EGP ----

model_egp <- cmdstan_model("egp.stan", cpp_options = list(stan_threads = TRUE))

season <- yday(rain$t) / (365 + leap_year(rain$t)) * 2*pi

thresh <- 0.5 * 1e-3 / sd_rain # 0.xx mm threshold

exceed <- rain$y - thresh
exceed[exceed < 0] <- 0

shards <- 8

opt_egp <- model_egp$optimize(
  data = list(K = shards,
              N = length(rain$y),
              y = rain$y,
              s = season,
              u = thresh,
              egp_id = 3,
              season_form = 1),
  threads = shards)



# SWG ----

model_swg <- cmdstan_model("swg.stan", cpp_options = list(stan_threads = TRUE))

len <- length(rain$t)
shards <- 2

opt_swg <- model_swg$optimize(
  data = list(K = shards,
              N = length(rain$y[1:len]),
              y = rain$y[1:len],
              s = season[1:len],
              u = thresh,
              egp_estim = opt_egp$mle(),
              vc_dim = nrow(cop_id),
              cop_id = cop_id,
              rotation = rotation),
  init = list(list(lambda_a = 0.05,
                   alpha = 0.05,
                   lambda_b = 7,
                   beta = 0.05,
                   gamma = 3,
                   delta1 = 2,
                   delta2 = 0.3,
                   omega1 = 7,
                   omega2 = 0.005,
                   omega3 = 0.002)),
  threads = shards,
  refresh = 20
)



# Simulation ----

sim_swg <- function(n, par, egp_par, mempreci = 1e-4, diag_intens = FALSE, obs_marks = NULL) {
  
  library(copula)
  library(rvinecopulib)
  library(progressr)
  library(DescTools)
  
  tic()
  
  par <- c(egp_par, par)
  
  # Generalized Pareto formulas.
  gp_pdf <- function(m, sigma, xi) {
    # GP pdf.
    if (xi != 0)
      dens <- 1/sigma * (1 + xi*m/sigma)^(-1/xi-1)
    else
      dens <- 1/sigma * exp(-m/sigma)
    return(dens)
  }
  gp_cdf <- function(m, sigma, xi) {
    # GP cdf.
    if (xi != 0)
      prob <- 1 - (1 + xi*m/sigma)^(-1/xi)
    else
      prob <- 1 - exp(-m/sigma)
    return(prob)
  }
  gp_q <- function(p, sigma, xi) {
    # GP quantile function.
    if (xi != 0)
      quant <- sigma * ((1 - p)^(-xi) - 1)/xi
    else
      quant <- sigma * log(1 - p)
    return(quant)
  }
  
  t <- 1:n
  memkern_alpha <- exp(par["alpha"] - par["alpha"] * rev(t))
  memkern_gamma <- exp(par["gamma"] - par["gamma"] * rev(t))
  
  memkern_alpha <- memkern_alpha[memkern_alpha > mempreci]
  memkern_gamma <- memkern_gamma[memkern_gamma > mempreci]
  memlen_alpha <- length(memkern_alpha)
  memlen_gamma <- length(memkern_gamma)
  
  impact <- function(ti, m, memlen, delta1, delta2) {
    # Impact function of the marked Hawkes process.
    out <- rep(NA, min(memlen, ti - 1))
    for (i in 1:length(out))
      out[i] = delta1 * (m[ti - 1 - length(out) + i])^delta2;
    return(out)
  }
  
  logis1 <- function(x, par1, par2) {
    # https://en.wikipedia.org/wiki/Logistic_function
    1 / (1 + exp(-par2*(x - par1)))
  }
  
  logis2_cdf <- function(x, omega1, omega2, omega3) {
    # cdf of the generalized logistic type II
    # https://en.wikipedia.org/wiki/Generalized_logistic_distribution
    1 - exp(-omega3 * (x - omega1)/omega2) / (1 + exp(-(x - omega1)/omega2))^omega3
  }
  
  lambda <- function(x, ti) {
    lambda_alphastar <- par["lambda_a"] * (1 - sum(tail(memkern_alpha, min(ti - 1, memlen_alpha)) * tail(x[1:(ti - 1)] > 0, min(ti - 1, memlen_alpha))))
    
    imp <- impact(ti, x, memlen_gamma, par["delta1"], par["delta2"])
    
    lambda_betastar <- par["lambda_b"] * exp(par["beta"] * (par["lambda_a"] - lambda_alphastar)) + sum(tail(memkern_gamma, min(ti - 1, memlen_gamma)) * imp * tail(x[1:(ti - 1)] > 0, min(ti - 1, memlen_gamma)))

    return(list(
      p = logis2_cdf(lambda_betastar, par["omega1"], par["omega2"], par["omega3"]),
      lambda_alphastar = lambda_alphastar,
      lambda_betastar = lambda_betastar)
    )
  }
  
  egp3_pdf <- function(y, sigma, xi, k1, k2, lb) {
    # EGP-beta pdf
    ub <- (k1*k2-1) / (k2-2)
    term1 <- 1.0 / (pbeta(0.5, k1, k2) - pbeta(lb, k1, k2))
    term2 <- dbeta((0.5 - lb) * gp_cdf(y, sigma, xi) + lb, k1, k2)
    dens <- term1 * term2 * gp_pdf(y, sigma, xi)
    dens <- dens * (0.5 - lb)
    return(dens)
  }
  egp3_cdf <- function(y, sigma, xi, k1, k2, lb) {
    # EGP-beta cdf
    ub <- (k1*k2-1) / (k2-2)
    p <- gp_cdf(y, sigma, xi)
    prob <- (pbeta((ub-lb)*p + lb, k1, k2) - pbeta(lb, k1, k2)) / (pbeta(ub, k1, k2) - pbeta(lb, k1, k2))
    return(prob)
  }
  egp3_q <- function(p, sigma, xi, k1, k2, lb) {
    # EGP-beta quantile function
    ub <- (k1*k2-1) / (k2-2)
    quant <- gp_q(
      (qbeta(pbeta(lb, k1, k2) + p*(pbeta(ub, k1, k2) - pbeta(lb, k1, k2)), k1, k2) - lb) / (ub - lb),
      sigma, xi
    )
    return(quant)
  }
  
  egp4_pdf <- function(y, sigma, xi, k1, k2, lb, k3) {
    # EGP-genbeta pdf
    gbeta1k_pdf <- function(u, a, b, c) c * beta(a, b)^(-1) * u^(a*c-1) * (1-u^c)^(b-1)
    gbeta1k_cdf <- function(u, a, b, c) pbeta(u^c, a, b)
    ub <-  ((k1*k3 - 1) / (k1*k3 + k2*k3 - k3 - 1))^(1/k3)
    term1 <- gbeta1k_pdf(gp_cdf(y, sigma, xi)*(ub-lb) + lb, k1, k2, k3)*(ub-lb) / (gbeta1k_cdf(ub, k1, k2, k3) - gbeta1k_cdf(lb, k1, k2, k3))
    dens <- term1 * gp_pdf(y, sigma, xi)
    return(dens)
  }
  egp4_cdf <- function(y, sigma, xi, k1, k2, lb, k3) {
    # EGP-genbeta cdf
    gbeta1k_cdf <- function(u, a, b, c) pbeta(u^c, a, b)
    ub <- ((k1*k3 - 1) / (k1*k3 + k2*k3 - k3 - 1))^(1/k3)
    prob <- (gbeta1k_cdf(gp_cdf(y, sigma, xi)*(ub-lb) + lb, k1, k2, k3) - gbeta1k_cdf(lb, k1, k2, k3)) / (gbeta1k_cdf(ub, k1, k2, k3) - gbeta1k_cdf(lb, k1, k2, k3))
    return(prob)
  }
  egp4_q <- function(p, sigma, xi, k1, k2, lb, k3) {
    # EGP-genbeta quantile function
    gbeta1k_cdf <- function(u, a, b, c) pbeta(u^c, a, b)
    gbeta1k_q <- function(p, a, b, c) qbeta(p, a, b)^(1/c)
    ub <- ((k1*k3 - 1) / (k1*k3 + k2*k3 - k3 - 1))^(1/k3)
    quant <- gp_q(
      (gbeta1k_q(gbeta1k_cdf(lb, k1, k2, k3) + p*(gbeta1k_cdf(ub, k1, k2, k3) - gbeta1k_cdf(lb, k1, k2, k3)), k1, k2, k3) - lb) / (ub - lb),
      sigma, xi
    )
    return(quant)
  }
  
  theta1 <- par[which(substr(names(par), 1, 7) == "theta1[")]
  theta2 <- par[which(substr(names(par), 1, 7) == "theta2[")]
  
  d_paircop <- length(theta1)
  
  vcdim <- par["vcdim"]
  
  cop_id <- matrix(par[which(substr(names(par), 1, 5) == "copid")], vcdim, vcdim)
  rotation <- matrix(par[which(substr(names(par), 1, 5) == "rotat")], vcdim, vcdim)
  vc_mat <- matrix(par[which(substr(names(par), 1, 5) == "vcmat")], vcdim, vcdim)
  cop_pos <- matrix(par[which(substr(names(par), 1, 6) == "coppos")], vcdim, vcdim)
  
  # Create the vine copula object for each possible exceedance or none pattern (0 meaning none).
  exc_pattern <- list()
  for (i in 1:vcdim)
    exc_pattern[[i]] <- 0:1
  exc_pattern <- expand.grid(exc_pattern)
  exc_pattern <- exc_pattern[exc_pattern[, 1] == 1, ] # First dim of the canonical VC always has an exceedance.
  
  bicop_code <- c("0" = "indep",
                  "1" = "gaussian",
                  "3" = "clayton",
                  "4" = "gumbel",
                  "5" = "frank",
                  "6" = "joe",
                  "91" = "bb1",
                  "96" = "bb6",
                  "97" = "bb7",
                  "98" = "bb8")
  
  # List of canonical vine copulas for each missingness pattern.
  vc_ls <- list()
  for (h in 1:nrow(exc_pattern)) {
    bicops <- list()
    if (any(exc_pattern[h, ] == 0)) {
      cop_pos_sub <- cop_pos[-which(exc_pattern[h, ] == 0), -which(exc_pattern[h, ] == 0)]
      cop_id_sub <- cop_id[-which(exc_pattern[h, ] == 0), -which(exc_pattern[h, ] == 0)]
    } else {
      cop_pos_sub <- cop_pos
      cop_id_sub <- cop_id
    }
    vcdim_sub <- sum(exc_pattern[h, ])
    
    if (vcdim_sub == 1) {
      vc_ls[[h]] <- NA
    }
    else {
      for (i in 1:(vcdim_sub - 1)) { # iterates vcmat rows and vine copula trees
        temp <- list()
        for (j in vcdim_sub:(i + 1)) { # iterates vcmat cols and ith tree pair-copulas
          bicoppars <- c(theta1[cop_pos_sub[i, j]], theta2[cop_pos_sub[i, j]])
          bicoppars <- bicoppars[bicoppars != 0]
          if (rotation[i, j] == 0)
            rot <- 0
          else if (rotation[i, j] == 1)
            rot <- 90
          temp <- list.append(temp, bicop_dist(
            family = bicop_code[as.character(cop_id_sub[i, j])],
            rotation = rot,
            parameters = bicoppars))
        }
        bicops[[i]] <- temp
      }
      vc_ls[[h]] <- vinecop_dist(pair_copulas = bicops, structure = cvine_structure(vcdim_sub:1))
    }
  }
  
  # Sequence on [0, 1] where the pdf of the vine copula is evaluated to find the maxima.
  # More points close to 0 and 1 because the density could spike there.
  len_seqd <- 20
  seqd <- c(exp(seq(-5, 0, length.out = len_seqd/2)) / 2,
            1 - exp(seq(0, -5, length.out = len_seqd/2)[-1]) / 2)
  
  cegp_rand <- function(p_prev, sigma, xi, k_vec) {
    
    if (all(p_prev == 0)) # No previous exceedance covered by the vine copula.
      draw <- runif(1)
    else {
      excpatt <- c(1, p_prev > 0)
      idx_patt <- 1 # find the exceedance pattern index (stupid way to do that but works ...)
      for (p in 1:nrow(exc_pattern)) {
        if (!all(exc_pattern[p, ] == excpatt))
          idx_patt <- idx_patt + 1
        else
          break
      }
      # Rejection sampling in the vine copula pdf.
      p_prev_exc <- p_prev[p_prev > 0]
      vc <- vc_ls[[idx_patt]]
      
      dens <- rep(NA, length(seqd))
      for (d in 1:length(dens))
        dens[d] <- dvinecop(u = c(1, p_prev_exc), vinecop = vc)
      up <- max(dens)
      d_propos <- -1
      while (runif(1, min = 0, max = up) > d_propos) {
        draw <- runif(1)
        d_propos <- dvinecop(u = c(draw, p_prev_exc), vine = vc)
      }
    }
    
    if (par["egpid"] == 3) {
      out <- egp3_q(p = draw, sigma, xi, k_vec[1], k_vec[2], k_vec[3])
    } else if (par["egpid"] == 4) {
      out <- egp4_q(p = draw, sigma, xi, k_vec[1], k_vec[2], k_vec[3], k_vec[4])
    }
    return(out)
  }
  
  tvec <- seq.POSIXt(as.POSIXct("0000-01-01"), by = "hour", length.out = n)
  season <- yday(tvec) / (365 + leap_year(tvec)) * 2*pi
  
  if (par["seasonform"] == 1)
    sigma_s <- par["sigma_0"] + par["season_a"] * (1/par["season_t"] * tanh(par["season_t"]*sin(season - par["season_d"]) / (1 - par["season_t"]*cos(season - par["season_d"]))))
  else
    sigma_s <- rep(par["sigma_0"], n)
  
  if (par["seasonform"] == 2)
    xi_s <- par["xi_0"] + par["season_a"] * (1/par["season_t"] * tanh(par["season_t"]*sin(season - par["season_d"]) / (1 - par["season_t"]*cos(season - par["season_d"]))))
  else
    xi_s <- rep(par["xi_0"], n)
  
  sim <- rep(NA, n)
  sim[1:(d_paircop - 1)] <- 0
  ptvec <- rep(NA, n)
  ptvec[1:(d_paircop - 1)] <- 0
  lambda1 <- rep(NA, n)
  lambda1[1:(d_paircop - 1)] <- 0
  lambda2 <- rep(NA, n)
  lambda2[1:(d_paircop - 1)] <- 0
  p_egp <- rep(NA, n)
  p_egp[1:(d_paircop - 1)] <- 0
  
  print("computing the simulations")
  pb <- txtProgressBar(min = 1, max = n, initial = 1)
  for (i in vcdim:n) {
    
    setTxtProgressBar(pb, i)
    
    temp <- lambda(sim, i)
    
    p <- temp[[1]]
    
    lambda1[i] <- temp[[2]]
    lambda2[i] <- temp[[3]]
    
    exceed <- runif(1) < p
    
    ptvec[i] <- p
    
    if (exceed) {
      
      sim[i] <- cegp_rand(p_prev = p_egp[(i - 1):(i - vcdim + 1)],
                          sigma = sigma_s[i], xi = xi_s[i],
                          k_vec = c(par["kappa1"], par["kappa2"], par["kappa3"], par["kappa4"])
      )
      if (par["egpid"] == 3) {
        p_egp[i] <- egp3_cdf(sim[i], sigma_s[i], xi_s[i],
                             par["kappa1"], par["kappa2"], par["kappa3"])
      } else if (par["egpid"] == 4) {
        p_egp[i] <- egp4_cdf(sim[i], sigma_s[i], xi_s[i],
                             par["kappa1"], par["kappa2"], par["kappa3"], par["kappa4"])
      }
    } else {
      sim[i] <- 0
      p_egp[i] <- 0
    }
  }
  close(pb)
  
  # Computing observations intensity. No need to compute for different simu or copula, so done here...
  if (diag_intens) {
    n_obs <- length(obs_marks)
    ptvec_obs <- rep(NA, n_obs)
    lambda1_obs <- rep(NA, n_obs)
    lambda2_obs <- rep(NA, n_obs)
    
    ptvec_obs[1] <- 0
    lambda1_obs[1] <- 0
    lambda2_obs[1] <- 0
    
    print("computing the obs intensity")
    pb <- txtProgressBar(min = 1, max = n_obs, initial = 1) 
    
    for (i in 2:n_obs) {
      setTxtProgressBar(pb, i)
      temp <- lambda(obs_marks, i)
      ptvec_obs[i] <- temp[[1]]
      lambda1_obs[i] <- temp[[2]]
      lambda2_obs[i] <- temp[[3]]
    }
    close(pb)
  }
  toc()
  
  if (diag_intens) {
    return(list(tvec = tvec,
                mark = sim,
                prob = ptvec,
                lambda1 = lambda1,
                lambda2 = lambda2,
                pegp = p_egp,
                prob_o = ptvec_obs,
                lambda1_o = lambda1_obs,
                lambda2_o = lambda2_obs))
  } else {
    return(list(tvec = tvec,
                mark = sim,
                prob = ptvec,
                lambda1 = lambda1,
                lambda2 = lambda2,
                pegp = p_egp))
  }
}

simu <- sim_swg(
  n = length(exceed),
  par = opt_swg$mle(),
  egp_par = opt_egp$mle(),
  diag_intens = FALSE
)