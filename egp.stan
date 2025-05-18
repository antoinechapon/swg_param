// extended generalized Pareto
functions {
  // _helpers ----
  int count_elem(array[] int test, int elem) {
    // Count number times elem appears in test set.
    // https://discourse.mc-stan.org/t/dealing-with-data-subsetting-in-stan/9842/8
    int count;
    count = 0;
    for (i in 1:num_elements(test))
      if (test[i] == elem)
        count = count + 1;
    return(count);
  }
  array[] int which_elem(array[] int test, int elem) {
    // Find elements in test which are equal to elem.
    // https://discourse.mc-stan.org/t/dealing-with-data-subsetting-in-stan/9842/8
    array[count_elem(test, elem)] int res;
    int ci;
    ci = 1;
    for (i in 1:num_elements(test))
      if (test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
  
  // _generalized Pareto ----
  real gp_pdf(real m, real sigma, real xi) {
    // GP pdf.
    real dens;
    if (xi != 0)
      dens = 1/sigma * (1 + xi*m/sigma)^(-1/xi-1);
    else
      dens = 1/sigma * exp(-m/sigma);
    return dens;
  }
  real gp_cdf(real m, real sigma, real xi) {
    // GP cdf.
    real prob;
    if (xi != 0)
      prob = 1 - (1 + xi*m/sigma)^(-1/xi);
    else
      prob = 1 - exp(-m/sigma);
    return prob;
  }
  real gp_q(real p, real sigma, real xi) {
    // GP quantile function.
    real quant;
    if (xi != 0)
      quant = sigma * ((1 - p)^(-xi) - 1)/xi;
    else
      quant = sigma * log(1 - p);
    return quant;
  }

  // _e0: no EGP ----
  real egp0_pdf(real d_gp, real p_gp, real kappa) {
    // No extending function pdf.
    return d_gp;
  }
  
  // _e1: EGP-power ----
  real egp1_pdf(real d_gp, real p_gp, real kappa) {
    // EGP-power pdf.
    real dens;
    dens = kappa * p_gp^(kappa-1) * d_gp;
    return dens;
  }

  // _e2: EGP-normal ----
  real egp2_pdf(real d_gp, real p_gp, real kappa) {
    // EGP-normal pdf (Gamet 2022 DOI: 10.1002/env.2744).
    real dens;
    if (kappa == 0)
      dens = 1;
    else {
      real term1;
      real term2;
      term1 = 2*sqrt(kappa) / (2*normal_cdf(sqrt(kappa) | 0, 1) - 1);
      term2 = exp(normal_lpdf(sqrt(kappa)*(p_gp - 1) | 0, 1));
      dens =  term1 * term2 * d_gp;
    }
    return dens;
  }
  
  // _e3: EGP-beta ----
  real egp3_pdf(real d_gp, real p_gp, real k1, real k2, real lb) {
    real dens;
    real ub;
    ub = (k1*k2-1) / (k2-2); // point in [0, 1] where derivative of the density is 0
    real term1;
    term1 = exp(beta_lpdf(p_gp*(ub-lb) + lb | k1, k2)) * (ub-lb) / (beta_cdf(ub | k1, k2) - beta_cdf(lb | k1, k2));
    dens = term1 * d_gp;
    return dens;
  }
  
  // _e4: EGP-genbeta ----
  real gbeta1k_pdf(real u, real a, real b, real c) {
    // Generalized beta of the first kind pdf. Alexander 2012 doi:10.1016/j.csda.2011.11.015
    real dens;
    dens = c * beta(a, b)^(-1) * u^(a*c-1) * (1-u^c)^(b-1);
    return dens;
  }
  real gbeta1k_cdf(real u, real a, real b, real c) {
    // Generalized beta of the first kind cdf.
    real prob;
    prob = beta_cdf(u^c | a, b);
    return prob;
  }
  real egp4_pdf(real d_gp, real p_gp, real k1, real k2, real lb, real k3) {
    // EGP-genbeta pdf.
    real dens;
    real ub; // upper truncation, point in (0, 1) where the derivative is 0
    ub = ((k1*k3 - 1) / (k1*k3 + k2*k3 - k3 - 1))^(1/k3);
    real term1;
    term1 = gbeta1k_pdf(p_gp*(ub-lb) + lb, k1, k2, k3) * (ub-lb) / (gbeta1k_cdf(ub | k1, k2, k3) - gbeta1k_cdf(lb | k1, k2, k3));
    dens = term1 * d_gp;
    return dens;
  }
  
  real egp_pdf(
    real m, real s,
    int egp_id, int season_form,
    real sigma_0, real xi_0,
    real season_a, real season_t, real season_d,
    real k1, real k2, real lb, real k3) {
    // EGP pdf wrapper.
    
    real sigma;
    if (season_form == 1) { // https://math.stackexchange.com/questions/2430564/equation-of-a-tilted-sine
      if (season_t == 0)
        sigma = sigma_0 + season_a * sin(s - season_d);
      else
        sigma = sigma_0 + season_a * 1/season_t * tanh(season_t*sin(s - season_d) / (1 - season_t*cos(s - season_d)));
    }
    else
      sigma = sigma_0;
    
    real xi;
    if (season_form == 2) {
      if (season_t == 0)
        xi = xi_0 + season_a * sin(s - season_d);
      else
        xi = xi_0 + season_a * 1/season_t * tanh(season_t*sin(s - season_d) / (1 - season_t*cos(s - season_d)));
    }
    else
      xi = xi_0;
      
    real gp_dens;
    real gp_prob;
    gp_dens = gp_pdf(m, sigma, xi);
    gp_prob = gp_cdf(m | sigma, xi);
    
    real dens;
    if (egp_id == 0) {// no EGP
      dens = egp0_pdf(gp_dens, gp_prob, k1);
    }
    else if (egp_id == 1) { // EGP-power
      dens = egp1_pdf(gp_dens, gp_prob, k1);
    }
    else if (egp_id == 2) { // EGP-normal
      dens = egp2_pdf(gp_dens, gp_prob, k1);
    }
    else if (egp_id == 3) { // EGP-beta
      dens = egp3_pdf(gp_dens, gp_prob, k1, k2, lb);
    }
    else if (egp_id == 4) { // EGP-genbeta
      dens = egp4_pdf(gp_dens, gp_prob, k1, k2, lb, k3);
    }
    else if (egp_id == 5) { // EGP-beta (Gamet and Jalbert article)
      dens = egp3_pdf(gp_dens, gp_prob, k1, k1, 1.0/32.0); // kappa1==kappa2
    }
    return dens;
  }
  
  vector egp_mr(vector param, vector dummy, data array[] real x_r, data array[] int x_i) {
    // Mapping function for parallelism.
    int n;
    n = dims(x_r)[1] %/% 2; // 2nd half of the x_r vector is the season covariate
    vector[n] egp_dens;
    
    for (i in 1:n) {
      if (x_r[i] > 0) {
        egp_dens[i] = egp_pdf(
          x_r[i], // ith mark
          x_r[i + n], // season covariate of the ith mark
          x_i[1], x_i[2],
          param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9]
        );
      }
      else {
        egp_dens[i] = 1; // removed in log
      }
    }
    return [sum(log(egp_dens))]';
  }
  
}

data {
  int K; // number of shards for parallelization
  int N;
  vector[N] y; // time series
  vector[N] s; // season on [0, 2*pi]
  
  real u; // threshold
  int egp_id; // identifier of the distribution extending the GP
  int season_form; //0: no season effect, 1: sigma, 2: xi
}

transformed data {
  vector[N] m; // mark
  m = y - u;
  array[N] int exc; // counting exceedances to compute the BIC
  for (i in 1:N) {
    if (m[i] < 0) {
      m[i] = 0;
    }
    exc[i] = m[i] > 0;
  }
  
  int<lower=0> J = N %/% K; // length of each shard subset
  array[K, J * 2] real x_r; // times 2 to add the season
  array[K, 2] int x_i;
  {
    int pos = 1;
    for (k in 1:K) {
      int end = pos + J - 1;
      
      matrix[J, 2] temp;
      temp[, 1] = m[pos:end];
      temp[, 2] = s[pos:end];
      
      x_r[k, ] = to_array_1d(temp);
      
      pos += J;
      x_i[k, 1] = egp_id;
      x_i[k, 2] = season_form;
    }
  }
  
  real eps;
  eps = 1e-9;
}

parameters {
  real sigma_0; // GP scale
  real xi_0; // GP shape
  
  real<lower=0> season_a; // amplitude
  real<lower=-1, upper=1> season_t; // skewness
  real<lower=0> season_d; // phase shift
  
  real<lower=0+eps, upper=1> kappa1;
  real<lower=0+eps, upper=1> kappa2;
  real<lower=0+eps, upper=0.1> kappa3; // lb for truncated beta distribution
  real<lower=0+eps, upper=1> kappa4;
  
  array[K] vector[0] dummy;
}

model {
  if (season_form == 1) {
    sigma_0 ~ normal(0, 1);
    season_a ~ gamma(2, 1);
    season_t ~ normal(0, 1);
    season_d ~ gamma(2, 1);
  }
  else
    sigma_0 ~ gamma(2, 1);
  
  if (season_form == 2) {
    season_a ~ gamma(2, 1);
    season_t ~ normal(0, 1);
    season_d ~ gamma(2, 1);
  }
  xi_0 ~ normal(0, 0.2);
  
  kappa1 ~ gamma(5, 15); // for param in [0, 1]
  kappa2 ~ gamma(5, 15);
  kappa3 ~ gamma(2, 40); // for lb
  kappa4 ~ gamma(5, 15);
  
  target += sum(
    map_rect(
      egp_mr,
      [sigma_0, xi_0,
      season_a, season_t, season_d,
      kappa1, kappa2, kappa3, kappa4]',
      dummy, x_r, x_i)
  );
}


generated quantities {
  int egpid;
  egpid = egp_id;
  int seasonform;
  seasonform = season_form;
  
 real loglik;
 loglik = sum(
    map_rect(
      egp_mr,
      [sigma_0, xi_0,
      season_a, season_t, season_d,
      kappa1, kappa2, kappa3, kappa4]',
      dummy, x_r, x_i)
  );

 int npar_egp;
 if (egp_id == 0) {
   npar_egp = 0;
 } else if (egp_id == 1) {
   npar_egp = 1;
 } else if (egp_id == 2) {
   npar_egp = 1;
 } else if (egp_id == 3) {
   npar_egp = 3;
 } else if (egp_id == 4) {
   npar_egp = 4;
 } else if (egp_id == 5) {
   npar_egp = 1;
 }
 
 if (season_form != 0)
  npar_egp = npar_egp + 3;
 
 real AIC;
 AIC = 2*npar_egp - 2*loglik;
 real BIC;
 BIC = (2 + npar_egp) * log(count_elem(exc, 1)) - 2*loglik;
}
