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
  real egp0_cdf(real p_gp, real kappa) {
    // No extending function cdf.
    return p_gp;
  }
  
  // _e1: EGP-power ----
  real egp1_pdf(real d_gp, real p_gp, real kappa) {
    // EGP-power pdf.
    real dens;
    dens = kappa * p_gp^(kappa-1) * d_gp;
    return dens;
  }
  real egp1_cdf(real p_gp, real kappa) {
    // EGP-power cdf.
    real prob;
    prob = p_gp^kappa;
    return prob;
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
  real egp2_cdf(real p_gp, real kappa) {
    // EGP-normal cdf.
    real prob;
    prob = 2/(2*normal_cdf(sqrt(kappa) | 0, 1)-1) * (normal_cdf(sqrt(kappa)*(p_gp - 1) | 0, 1) - (1-normal_cdf(sqrt(kappa) | 0, 1)) );
    return prob;
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
  real egp3_cdf(real p_gp, real k1, real k2, real lb) {
    real prob;
    real ub;
    ub = (k1*k2-1) / (k2-2); // point in [0, 1] where derivative of the density is 0
    prob = (beta_cdf(p_gp*(ub-lb) + lb | k1, k2) - beta_cdf(lb | k1, k2)) / (beta_cdf(ub | k1, k2) - beta_cdf(lb | k1, k2));
    return prob;
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
  real egp4_cdf(real p_gp, real k1, real k2, real lb, real k3) {
    // EGP-genbeta cdf.
    real prob;
    real ub; // upper truncation, point in (0, 1) where the derivative is 0
    ub = ((k1*k3 - 1) / (k1*k3 + k2*k3 - k3 - 1))^(1/k3);
    prob = (gbeta1k_cdf(p_gp*(ub-lb) + lb | k1, k2, k3) - gbeta1k_cdf(lb | k1, k2, k3)) / (gbeta1k_cdf(ub | k1, k2, k3) - gbeta1k_cdf(lb | k1, k2, k3));
    return prob;
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
    
    real egp_dens;
    if (egp_id == 0) // no EGP
      egp_dens = egp0_pdf(gp_dens, gp_prob, k1);
    else if (egp_id == 1) // EGP-power
      egp_dens = egp1_pdf(gp_dens, gp_prob, k1);
    else if (egp_id == 2) // EGP-normal
      egp_dens = egp2_pdf(gp_dens, gp_prob, k1);
    else if (egp_id == 3) // EGP-beta
      egp_dens = egp3_pdf(gp_dens, gp_prob, k1, k2, lb);
    else if (egp_id == 4) // EGP-genbeta
      egp_dens = egp4_pdf(gp_dens, gp_prob, k1, k2, lb, k3);
    else if (egp_id == 5) // EGP-beta (Gamet and Jalbert article)
      egp_dens = egp3_pdf(gp_dens, gp_prob, k1, k1, 1.0/32.0); // kappa1==kappa2
      
    return egp_dens;
  }
  
  real egp_cdf(
    real m, real s,
    int egp_id, int season_form,
    real sigma_0, real xi_0,
    real season_a, real season_t, real season_d,
    real k1, real k2, real lb, real k3) {
    // EGP cdf wrapper.
    
    real sigma;
    if (season_form == 1) { // https://math.stackexchange.com/questions/2430564/equation-of-a-tilted-sine
      if (season_t == 0)
        sigma = sigma_0 + season_a * sin(s - season_d);
      else
        sigma = sigma_0 + season_a * (1/season_t * tanh(season_t*sin(s - season_d) / (1 - season_t*cos(s - season_d))));
    }
    else
      sigma = sigma_0;
    
    real xi;
    if (season_form == 2) {
      if (season_t == 0)
        xi = xi_0 + season_a * sin(s - season_d);
      else
        xi = xi_0 + season_a * (1/season_t * tanh(season_t*sin(s - season_d) / (1 - season_t*cos(s - season_d))));
    }
    else
      xi = xi_0;

    real gp_prob;
    gp_prob = gp_cdf(m | sigma, xi);
    
    real egp_prob;
    if (egp_id == 0)// no EGP
      egp_prob = egp0_cdf(gp_prob | k1);
    else if (egp_id == 1) // EGP-power
      egp_prob = egp1_cdf(gp_prob | k1);
    else if (egp_id == 2) // EGP-normal
      egp_prob = egp2_cdf(gp_prob | k1);
    else if (egp_id == 3) // EGP-beta
      egp_prob = egp3_cdf(gp_prob | k1, k2, lb);
    else if (egp_id == 4) // EGP-genbeta
      egp_prob = egp4_cdf(gp_prob | k1, k2, lb, k3);
    else if (egp_id == 5) // EGP-beta (Gamet and Jalbert article)
      egp_prob = egp3_cdf(gp_prob | k1, k1, 1.0/32.0); // kappa1==kappa2
      
    return egp_prob;
  }
  
  
  // Schepsmeier 2014 DOI: 10.1007/s00362-013-0498-x for the pair-copula pdf.
  
  // _c0: independence copula ----
  real cop0_pdf(real u1, real u2) {
    // Independence copula pdf.
    return 1;
  }
  real cop0_hfun(real u1, real u2) {
    // Independence copula h-function.
    return u2;
  }
  
  // _c1: Gaussian copula ----
  real cop1_pdf(real u1, real u2, real theta) {
    // Gaussian copula pdf.
    real dens;
    real x1;
    real x2;
    x1 = inv_Phi(u1);
    x2 = inv_Phi(u2);
    dens = 1/sqrt(1-theta^2) * exp(-((theta^2*(x1^2+x2^2)-(2*theta*x1*x2)) / (2*(1-theta^2))));
    return dens;
  }
  real cop1_hfun(real u1, real u2, real theta) {
    real prob;
    prob = normal_cdf((inv_Phi(u2) - theta*inv_Phi(u1)) / sqrt(1 - theta^2) | 0, 1);
    return prob;
  }
  
  // _c3: Clayton copula ----
  real cop3_pdf(real u1, real u2, real theta) {
    // Clayton copula pdf.
    real dens;
    dens = (1 + theta) * (u1*u2)^(-1 - theta) * (u1^(-theta) + u2^(-theta) - 1)^(-1/theta - 2);
    return dens;
  }
  real cop3_hfun(real u1, real u2, real theta) {
    // Clayton copula h-function.
    real prob;
    prob = (1 + u1^theta * (u2^-theta - 1))^(-1 - 1/theta);
    return prob;
  }
  
  // _c4: Gumbel copula ----
  real cop4_pdf(real u1, real u2, real theta) {
    // Gumbel copula pdf.
    real dens;
    real t1;
    real t2;
    t1 = (-log(u1))^theta;
    t2 = (-log(u2))^theta;
    real cdf;
    cdf = exp(-(t1 + t2)^(1/theta));
    dens = cdf * 1/(u1*u2) * (t1+t2)^(-2+2/theta) * (log(u1)*log(u2))^(theta-1) * (1 + (theta-1)*(t1+t2)^(-1/theta));
    return dens;
  }
  real cop4_hfun(real u1, real u2, real theta) {
    // Gumbel copula h-function.
    real prob;
    real x;
    real y;
    x = -log(u1);
    y = -log(u2);
    prob = u1^-1 * exp(-(x^theta + y^theta)^(1/theta)) * (1 + (y/x)^theta)^(1/theta-1);
    return prob;
  }
  
  // _c5: Frank copula ----
  real cop5_pdf(real u1, real u2, real theta) {
    // Franck copula pdf.
    real dens;
    dens = (-theta*exp(-theta*(u1+u2))*(exp(-theta)-1)) / (exp(-theta) - exp(-theta*u1) - exp(-theta*u2) + exp(-theta*(u1+u2)))^2;
    return dens;
  }
  real cop5_hfun(real u1, real u2, real theta) {
    // Franck copula h-function.
    real prob;
    prob = exp(-theta*u1) * ((1-exp(-theta)) * (1-exp(-theta*u2))^-1 - (1-exp(-theta*u1)))^-1;
    return prob;
  }
  
  // _c6: Joe copula ----
  real cop6_pdf(real u1, real u2, real theta) {
    // Joe copula pdf.
    real dens;
    real t1;
    real t2;
    t1 = (1 - u1)^theta;
    t2 = (1 - u2)^theta;
    dens = (t1+t2-t1*t2)^(1/theta-2) * (theta-1+t1+t2-t1*t2) * (1-u1)^(theta-1) * (1-u2)^(theta-1);
    return dens;
  }
  real cop6_hfun(real u1, real u2, real theta) {
    // Joe copula h-function.
    real prob;
    prob = (1 + (1-u2)^theta * (1-u1)^-theta - (1-u2)^theta)^(-1 + 1/theta) * (1 - (1-u2)^theta);
    return prob;
  }
  
  // _c7: Galambos copula ----
  real cop7_pdf(real u1, real u2, real theta) {
    // Galambos copula pdf.
    real dens;
    real x;
    real y;
    x = -log(u1);
    y = -log(u2);
    real cdf;
    cdf = u1*u2*exp( ( (-log(u1))^-theta + (-log(u2))^-theta )^(-1/theta) );
    dens = cdf/(u1*u2) * (1 - (x^-theta + y^-theta)^(-1-1/theta) * (x^(-theta-1) + y^(-theta-1)) + (x^-theta + y^-theta)^(-2-1/theta) * (x*y)^(-theta-1) * (1+theta+(x^-theta+y^-theta)^(-1/theta)) );
    return dens;
  }
  real cop7_hfun(real u1, real u2, real theta) {
    // Galambos copula h-function.
    real prob;
    real x;
    real y;
    x = -log(u1);
    y = -log(u2);
    prob = u2 * exp((x^-theta + y^-theta)^(-1/theta)) * (1 - (1 + (x/y)^theta)^(-1-1/theta));
    return prob;
  }
  
  // _c8: Hüsler-Reiss copula ----
  real cop8_pdf(real u1, real u2, real theta) {
    // Hüsler-Reiss copula pdf.
    real dens;
    real x;
    real y;
    x = -log(u1);
    y = -log(u2);
    real cdf;
    cdf = exp(-x*normal_cdf(theta^-1 + 0.5*theta*log(x/y) | 0, 1) -y*normal_cdf(theta^-1 + 0.5*theta*log(y/x) | 0, 1));
    dens = cdf/(u1*u2) * (normal_cdf(theta^-1 + 0.5*theta*log(x/y) | 0, 1) * normal_cdf(theta^-1 + 0.5*theta*log(y/x) | 0, 1) + 0.5*theta*y^-1 * exp(normal_lpdf(theta^-1 + 0.5*theta*log(x/y)| 0, 1)));
    return dens;
  }
  real cop8_hfun(real u1, real u2, real theta) {
    // Hüsler-Reiss copula h-function.
    real prob;
    real x;
    real y;
    x = -log(u1);
    y = -log(u2);
    real cdf;
    cdf = exp(-x*normal_cdf(theta^-1 + 0.5*theta*log(x/y) | 0, 1) -y*normal_cdf(theta^-1 + 0.5*theta*log(y/x) | 0, 1));
    prob = cdf * u1^-1 * normal_cdf(theta^-1 + 0.5*theta*log(x/y) | 0, 1);
    return prob;
  }
  
  // _c91: BB1 copula ----
  real cop91_pdf(real u1, real u2, real t1, real t2) {
    // BB1 copula pdf.
    real dens;
    real x;
    real y;
    x = (u1^-t1-1)^t2;
    y = (u2^-t1-1)^t2;
    dens = (1+(x+y)^(1/t2))^(-1/t1-2) * (x+y)^(1/t2-2) * (t1*(t2-1) + (t1*t2+1) * (x+y)^(1/t2)) * (x*y)^(1-1/t2) * (u1*u2)^(-t1-1);
    return dens;
  }
  real cop91_hfun(real u1, real u2, real t1, real t2) {
  // BB1 copula h-function.
    real prob;
    real x;
    real y;
    x = (u1^-t1-1)^t2;
    y = (u2^-t1-1)^t2;
    prob = (1 + (x+y)^(1/t2))^(-1/t1-1) * (x+y)^(1/t2-1) * x^(1-1/t2) * u1^(-t1-1);
    return prob;
  }
  
  // _c96: BB6 copula ----
  real cop96_pdf(real u1, real u2, real t1, real t2) {
    // BB6 copula pdf.
    real dens;
    real x;
    real y;
    real w;
    x = -log(1-(1-u1)^t1);
    y = -log(1-(1-u2)^t1);
    w = exp(-(x^t2 + y^t2)^(1/t2));
    dens = (1-w)^(1/t1-2)*w * (x^t2 + y^t2)^(1/t2-2) * ((t1-w)*(x^t2 + y^t2)^(1/t2) + t1*(t2-1)*(1-w)) * (x*y)^(t2-1) * (1-(1-u1)^t1)^-1 * (1-(1-u2)^t1)^-1 * ((1-u1)*(1-u2))^(t1-1);
    return dens;
  }
  real cop96_hfun(real u1, real u2, real t1, real t2) {
  // BB6 copula h-function.
    real prob;
    real x;
    real y;
    real w;
    x = -log(1-(1-u1)^t1);
    y = -log(1-(1-u2)^t1);
    w = exp(-(x^t2 + y^t2)^(1/t2));
    prob = (1-w)^(1/t1-1) * w * (x^t2 + y^t2)^(1/t2-1) * x^(t2-1) * exp(x) * (1-exp(-x))^(1-1/t1);
    return prob;
  }
  
  // _c97: BB7 copula ----
  real cop97_pdf(real u1, real u2, real t1, real t2) {
    // BB7 copula pdf.
    real dens;
    real x;
    real y;
    x = (1 - (1 - u1)^t1)^(-t2) - 1;
    y = (1 - (1 - u2)^t1)^(-t2) - 1;
    dens = (1 - (x + y + 1)^(-1/t2))^(1/t1-2) * (x + y + 1)^(-1/t2-2) * ((x + 1)*(y + 1))^(1 + 1/t2) * (t1*(t2 + 1) - (t1*t2 + 1)*(x + y + 1)^(-1/t2)) * ((1 - u1)*(1 - u2))^(t1 - 1);
    return dens;
  }
  real cop97_hfun(real u1, real u2, real t1, real t2) {
  // BB7 copula h-function.
    real prob;
    real x;
    real y;
    x = (1 - (1 - u1)^t1)^(-t2) - 1;
    y = (1 - (1 - u2)^t1)^(-t2) - 1;
    prob = (1 - (x+y+1)^(-1/t2))^(1/t1-1) * (x+y+1)^(-1/t2-1) * (x+1)^(1+1/t2) * (1-u1)^(t1-1);
    return prob;
  }
  
  // _c98: BB8 copula ----
  real cop98_pdf(real u1, real u2, real t1, real t2) {
    // BB8 copula pdf.
    real dens;
    real eta;
    real x;
    real y;
    eta = 1 - (1 - t2)^t1;
    x = 1 - (1 - t2 * u1)^t1;
    y = 1 - (1 - t2 * u2)^t1;
    dens = 1/eta * t2*(1 - 1/eta * x*y)^(1/t1-2) * (t1 - 1/eta * x*y) * (1 - t2*u1)^(t1-1) * (1 - t2*u2)^(t1-1);
    return dens;
  }
  real cop98_hfun(real u1, real u2, real t1, real t2) {
  // BB8 copula h-function.
    real prob;
    real eta;
    real x;
    real y;
    eta = 1 - (1 - t2)^t1;
    x = 1 - (1 - t2 * u1)^t1;
    y = 1 - (1 - t2 * u2)^t1;
    prob = (eta^(-1) * y * (1 - eta^(-1)*x*y)^(1/t1 - 1)) / ((1 - x)^(1/t1 - 1));
    return prob;
  }
  
  vector cop_par_map(int cop_id, real theta1_g, real theta2_g) {
    // Maps the generic pair-copula paramters from [0, 1] to a family specific interval.
    real theta1;
    real theta2;
    if (cop_id == 0) { // independence
      theta1 = 0.0;
      theta2 = 0.0;
    }
    else if (cop_id == 1) { // Gaussian
      theta1 = theta1_g * 1.9 - 1 + 1e-9;
      theta2 = 0.0;
    }
    else if (cop_id == 3) { // Clayton
      theta1 = theta1_g * 5 + 1e-9;
      theta2 = 0.0;
    }
    else if (cop_id == 4) { // Gumbel
      theta1 = theta1_g * 5 + 1 + 1e-9;
      theta2 = 0.0;
    }
    else if (cop_id == 5) { // Frank
      theta1 = theta1_g * 20 - 10 + 1e-9;
      theta2 = 0.0;
    }
    else if (cop_id == 6) { // Joe
      theta1 = theta1_g * 4 + 1 + 1e-9;
      theta2 = 0.0;
    }
    else if (cop_id == 7) { // Galambos
      theta1 = theta1_g * 5 + 1e-9;
      theta2 = 0.0;
    }
    else if (cop_id == 8) { // Hüsler-Reiss
      theta1 = theta1_g * 5 + 1e-9;
      theta2 = 0.0;
    }
    else if (cop_id == 91) { // BB1
      theta1 = theta1_g * 5 + 1e-9;
      theta2 = theta2_g * 5 + 1 + 1e-9;
    }
    else if (cop_id == 96) { // BB6
      theta1 = theta1_g / 2 + 1 + 1e-9;
      theta2 = theta2_g / 2 + 1 + 1e-9;
    }
    else if (cop_id == 97) { // BB7
      theta1 = theta1_g * 1 + 1 + 1e-9;
      theta2 = theta2_g * 0.5 + 1e-9;
    }
    else if (cop_id == 98) { // BB8
      theta1 = theta1_g * 7 + 1 + 1e-9;
      theta2 = theta2_g * 2 + 1e-9;
    }
    return [theta1, theta2]';
  }
  
  real cop_pdf(real u1, real u2, int cop_id, int rotation, real theta1_g, real theta2_g) {
    // Copula pdf wrapper.
    real u1_r;
    real u2_r;
    if (rotation == 0) {
      u1_r = u1;
      u2_r = u2;
    }
    else if (rotation == 1) {
      u1_r = 1 - u1;
      u2_r = 1 - u2;
    }
    vector[2] theta;
    theta = cop_par_map(cop_id, theta1_g, theta2_g);
    real dens;
    if (cop_id == 0) // independence
      dens = cop0_pdf(u1_r, u2_r);
    else if (cop_id == 1) // Gaussian
      dens = cop1_pdf(u1_r, u2_r, theta[1]);
    else if (cop_id == 3) // Clayton
      dens = cop3_pdf(u1_r, u2_r, theta[1]);
    else if (cop_id == 4) // Gumbel
      dens = cop4_pdf(u1_r, u2_r, theta[1]);
    else if (cop_id == 5) // Frank
      dens = cop5_pdf(u1_r, u2_r, theta[1]);
    else if (cop_id == 6) // Joe
      dens = cop6_pdf(u1_r, u2_r, theta[1]);
    else if (cop_id == 7) // Galambos
      dens = cop7_pdf(u1_r, u2_r, theta[1]);
    else if (cop_id == 8) // Hüsler-Reiss
      dens = cop8_pdf(u1_r, u2_r, theta[1]);
    else if (cop_id == 91) // BB1
      dens = cop91_pdf(u1_r, u2_r, theta[1], theta[2]);
    else if (cop_id == 96) // BB6
      dens = cop96_pdf(u1_r, u2_r, theta[1], theta[2]);
    else if (cop_id == 97) // BB7
      dens = cop97_pdf(u1_r, u2_r, theta[1], theta[2]);
    else if (cop_id == 98) // BB8
      dens = cop98_pdf(u1_r, u2_r, theta[1], theta[2]);
    return dens;
  }
  real cop_hfun(real u1, real u2, int cop_id, int rotation, real theta1_g, real theta2_g) {
    // Copula h-function wrapper. Always probability of u2 given u1.
    real u1_r;
    real u2_r;
    if (rotation == 0) {
      u1_r = u1;
      u2_r = u2;
    }
    else if (rotation == 1) {
      u1_r = 1 - u1;
      u2_r = 1 - u2;
    }
    if (u1_r == 0.0)
      return u2_r;
    if (u2_r == 0.0)
      return 0.0;
    vector[2] theta;
    theta = cop_par_map(cop_id, theta1_g, theta2_g);
    real prob;
    if (cop_id == 0) // independence
      prob = cop0_hfun(u1_r, u2_r);
    else if (cop_id == 1) // Gaussian
      prob = cop1_hfun(u1_r, u2_r, theta[1]);
    else if (cop_id == 3) // Clayton
      prob = cop3_hfun(u1_r, u2_r, theta[1]);
    else if (cop_id == 4) // Gumbel
      prob = cop4_hfun(u1_r, u2_r, theta[1]);
    else if (cop_id == 5) // Frank
      prob = cop5_hfun(u1_r, u2_r, theta[1]);
    else if (cop_id == 6) // Joe
      prob = cop6_hfun(u1_r, u2_r, theta[1]);
    else if (cop_id == 7) // Galambos
      prob = cop7_hfun(u1_r, u2_r, theta[1]);
    else if (cop_id == 8) // Hüsler-Reiss
      prob = cop8_hfun(u1_r, u2_r, theta[1]);
    else if (cop_id == 91) // BB1
      prob = cop91_hfun(u1_r, u2_r, theta[1], theta[2]);
    else if (cop_id == 96) // BB6
      prob = cop96_hfun(u1_r, u2_r, theta[1], theta[2]);
    else if (cop_id == 97) // BB7
      prob = cop97_hfun(u1_r, u2_r, theta[1], theta[2]);
    else if (cop_id == 98) // BB8
      prob = cop98_hfun(u1_r, u2_r, theta[1], theta[2]);
    return prob;
  }
  
  
  // _canonical vine copula ----
  real cvc_pdf(
    array[,,] real effx_t, array[] int bool_miss,
    int vc_dim, array[,] int cop_id, array[,] int rotation, array[,] int vc_mat, array[,] int cop_pos,
    vector theta1_g, vector theta2_g) {
    // Canonical vine copula pdf, with truncation according to missing dimensions.
    matrix[vc_dim, vc_dim] dens_c;
    for (i in 1:vc_dim) { // Iterating vine copula trees, rows of the vc matrix.
      for (j in 1:vc_dim) {
        if (j > i && bool_miss[i] == 0 && bool_miss[j] == 0) {
          // Upper right corner and not missing dimension.
          dens_c[i, j] = cop_pdf(
            effx_t[i, j, 1], effx_t[i, j, 2], cop_id[i, j], rotation[i, j],
            theta1_g[cop_pos[i, j]],
            theta2_g[cop_pos[i, j]]
            );
        }
        else {
          dens_c[i, j] = 1; // 1s removed in matrix product
        }
      }
    }
    real dens_prod;
    dens_prod = prod(dens_c);
    return dens_prod;
  }
  
  real cegp_pdf(
    array[,,] real effx_t, real d0,
    int vc_dim,
    array[,] int cop_id, array[,] int rotation, array[,] int vc_mat, array[,] int cop_pos,
    vector theta1_g, vector theta2_g
    ) {
    // Assembles EGP and copula densities to get the CEGP density.
    real cegp_dens;
    real vc_dens;
    array[vc_dim] int bool_miss;
    for (i in 1:vc_dim)
      bool_miss[i] = effx_t[1, i, 2] == 0;
    
    vc_dens = cvc_pdf(effx_t, bool_miss, vc_dim, cop_id, rotation, vc_mat, cop_pos, theta1_g, theta2_g);
    
    cegp_dens = vc_dens * d0;
    return cegp_dens;
  }
  
  vector impact(int ti, vector m, int memlen, real delta1, real delta2) {
    // Impact function of the marked Hawkes process.
    vector[min(memlen, ti - 1)] out;
    for (i in 1:num_elements(out)) {
      out[i] = delta1 * (m[ti-1 - num_elements(out)+i])^delta2;
    }
    return out;
  }
  
  real logis2_cdf(real x, real omega1, real omega2, real omega3) {
    // cdf of the generalized logistic type II
    // https://en.wikipedia.org/wiki/Generalized_logistic_distribution
    real prob;
    prob = 1 - exp(-omega3 * (x - omega1)/omega2) / (1 + exp(-(x - omega1)/omega2))^omega3;
    return prob;
  }
  
  real hawkes_lpdf(
    vector m, vector exc, vector p_egp, vector d_egp,
    
    real lambda_a, real alpha,
    real lambda_b, real beta, real gamma,
    
    real delta1, real delta2,
    real omega1, real omega2, real omega3,
    
    int vc_dim,
    array[,] int cop_id, array[,] int rotation,
    array[,] int vc_mat, array[,] int cop_pos,
    vector theta1_g, vector theta2_g
    ) {
    
    real mempreci;
    mempreci = 1e-4;
    // Chosen to be close enough to 0 (compared to 1 at t-1).
    
    int n;
    n = num_elements(m);
    
    vector[n] t;
    for (i in 1:n) {
      t[i] = i;
    }
    
    vector[n] memt_alpha = exp(alpha - alpha * reverse(t)); // 'infinite' memory kernel
    vector[n] memt_gamma = exp(gamma - gamma * reverse(t));
    
    array[n] int lowmem_alpha;
    array[n] int  lowmem_gamma;
    for (i in 1:n) {
      lowmem_alpha[i] = memt_alpha[i] > mempreci;
      lowmem_gamma[i] = memt_gamma[i] > mempreci;
    }
    
    int n_mem_alpha;
    n_mem_alpha = count_elem(lowmem_alpha, 1);
    vector[n_mem_alpha] memkern_alpha; // censored memory kernel alpha
    memkern_alpha = memt_alpha[which_elem(lowmem_alpha, 1)];
    
    int n_mem_gamma;
    n_mem_gamma = count_elem(lowmem_gamma, 1);
    vector[n_mem_gamma] memkern_gamma; // censored memory kernel gamma
    memkern_gamma = memt_gamma[which_elem(lowmem_gamma, 1)];
    
    vector[n] lambda_betastar;
    for (i in 1:(vc_dim - 1))
      lambda_betastar[i] = 0;
    
    for (i in vc_dim:n) {
      int n_alpha;
      n_alpha = min(i - 1, n_mem_alpha);
      vector[n_alpha] inhib;
      inhib = tail(memkern_alpha, n_alpha) .* tail(exc[1:(i - 1)], n_alpha);
      real lambda_alphastar;
      lambda_alphastar = lambda_a * (1 - sum(inhib));
      vector[min(i - 1, n_mem_gamma)] imp;
      imp = impact(i, m, n_mem_gamma, delta1, delta2);
      int n_mem_trunc;
      n_mem_trunc = min(i - 1, n_mem_gamma);
      vector[n_mem_trunc] imp_mem;
      imp_mem = tail(memkern_gamma, n_mem_trunc) .* imp .* tail(exc[1:(i - 1)], n_mem_trunc);
      lambda_betastar[i] = lambda_b * exp(beta * (lambda_a - lambda_alphastar)) + sum(imp_mem);
    }
  
    // Effective probabilities for the vine copula.
    array[vc_dim, vc_dim, n, 2] real effx_vc;
    // Indexing works with the canonical vine copula matrix in upper right corner.
    // Second place of the 4th dimension of effx_vc is always on the diagonal of vc_mat.
    for (i in 1:vc_dim) {
      for (j in 1:vc_dim) {
        if (vc_mat[i, j] > 0 && j > i) {
          if (i == 1) { // first tree
            for (k in 1:n) {
              vector[vc_dim] x_vc;
              for (l in 1:vc_dim) {
                if (k < l)
                  x_vc[l] = 0;
                else
                  x_vc[l] = p_egp[k - l + 1];
              }
              effx_vc[i, j, k, 1] = x_vc[vc_mat[i, j]];
              effx_vc[i, j, k, 2] = x_vc[vc_mat[j, j]];
            }
          }
          else {
            for (k in 1:n) {
              effx_vc[i, j, k, 1] = cop_hfun(
                effx_vc[i-1, j-1, k, 1], effx_vc[i-1, j-1, k, 2],
                cop_id[i-1, j-1], rotation[i-1, j-1],
                theta1_g[cop_pos[i-1, j-1]], theta2_g[cop_pos[i-1, j-1]]
                );
              effx_vc[i, j, k, 2] = cop_hfun(
                effx_vc[i-1, j, k, 1], effx_vc[i-1, j, k, 2],
                cop_id[i-1, j], rotation[i-1, j],
                theta1_g[cop_pos[i-1, j]], theta2_g[cop_pos[i-1, j]]
                );
            }
          }
        }
      }
    }
    
    vector[n] ll1; // mark term of the log-likelihood
    vector[n] ll2; // point process term of the log-likelihood
    for (i in 1:(vc_dim-1)) {
      ll1[i] = 0;
      ll2[i] = 0;
    }
    for (i in vc_dim:n) { // Starts at t2 because of the copula-EGP.
      real pt; // exceedance probability at time t
      pt = logis2_cdf(lambda_betastar[i] | omega1, omega2, omega3);
      
      if (exc[i] == 1) {
        ll1[i] = log(
          cegp_pdf(effx_vc[ , , i, ], d_egp[i], vc_dim, cop_id, rotation, vc_mat,
          cop_pos, theta1_g, theta2_g)
          );
      }
      else { // y[i]<=u, so it does not contribute to the mark likelihood.
        ll1[i] = 0;
      }
      real term;
      term = (1 - exc[i]) * log(1 - pt);
      if (is_nan(term)) { // avoiding NaNs
        term = 0; // (1 - exc[i]) is always 0 when log(1 - pt[i]) is -inf
      }
      ll2[i] = exc[i] * log(pt) + term;
    }
    return sum(ll1) + sum(ll2);
  }
  
  
  vector hawkes_mr(vector param, vector dummy, data array[] real x_r, data array[] int x_i) {
    // Mapping function.
    vector[1] dens;
    int vcdim;
    vcdim = x_i[1];
    int J;
    J = dims(x_r)[1] %/% 4;
    
    int dpaircop; // numbers of pair-copulas in the vine
    dpaircop = vcdim * (vcdim - 1) / 2;
    
    array[vcdim, vcdim] int copid;
    array[vcdim, vcdim] int rotat;
    array[vcdim, vcdim] int vcmat;
    array[vcdim, vcdim] int coppos;
    for (i in 1:vcdim) {
      for (j in 1:vcdim) {
        copid[j, i] = x_i[(2+0*vcdim*vcdim):(1+1*vcdim*vcdim)][i+vcdim*(j-1)];
        rotat[j, i] = x_i[(2+1*vcdim*vcdim):(1+2*vcdim*vcdim)][i+vcdim*(j-1)];
        vcmat[j, i] = x_i[(2+2*vcdim*vcdim):(1+3*vcdim*vcdim)][i+vcdim*(j-1)];
        coppos[j, i] = x_i[(2+3*vcdim*vcdim):(1+4*vcdim*vcdim)][i+vcdim*(j-1)];
      }
    }
    
    dens[1] = hawkes_lpdf(
      to_vector(x_r[1:J]) | to_vector(x_r[(J+1):(2*J)]), to_vector(x_r[(2*J+1):(3*J)]), to_vector(x_r[(3*J+1):(4*J)]),
      
      param[1], param[2], param[3], param[4], param[5],
      param[6], param[7], param[8], param[9], param[10],
      
      vcdim, copid, rotat, vcmat, coppos,
      
      param[(11 + 0*dpaircop):(10 + 1*dpaircop)],
      param[(11 + 1*dpaircop):(10 + 2*dpaircop)]
    );
    
    return dens;
  }
  
}


data {
  int K; // number of shards for map_rect parallelization
  int N;
  vector[N] y; // time series
  vector[N] s; // season on [0, 2*pi]
  real u; // threshold
  
  vector[15] egp_estim; // estimates of the 9 param of the EGP plus the 2 model indicators
  
  int vc_dim;
  array[vc_dim, vc_dim] int cop_id; // identifier of the copula
  array[vc_dim, vc_dim] int rotation; // 0: no copula rotation, 1: 180° rotation
}

transformed data {
  
  int egp_id; // identifier of the distribution extending the GP
  int season_form; //0: no season effect, 1: sigma, 2: xi
  real sigma_0; // GP scale
  real xi_0; // GP shape
  real season_a; // amplitude
  real season_t; // skewness
  real season_d; // phase shift
  real kappa1; // parameters of the EGP probability integral transform
  real kappa2;
  real kappa3;
  real kappa4;
  
  egp_id = to_int(egp_estim[10]);
  season_form = to_int(egp_estim[11]);
  sigma_0 = egp_estim[1]; 
  xi_0 = egp_estim[2];
  season_a = egp_estim[3];
  season_t = egp_estim[4];
  season_d = egp_estim[5];
  kappa1 = egp_estim[6];
  kappa2 = egp_estim[7];
  kappa3 = egp_estim[8];
  kappa4 = egp_estim[9];
  
  vector[N] m; // mark
  vector[N] exc; // bool exceedance (as real to avoid casting)
  vector[N] p_egp; // probability of the EGP
  vector[N] d_egp; // density of the EGP
  m = y - u;
  for (i in 1:N) {
    if (m[i] < 0) {
      m[i] = 0;
      p_egp[i] = 0;
      d_egp[i] = 0;
    }
    else {
      p_egp[i] = egp_cdf(
        m[i] | s[i], egp_id, season_form,
        sigma_0, xi_0,
        season_a, season_t, season_d,
        kappa1, kappa2, kappa3, kappa4
      );
      d_egp[i] = egp_pdf(
        m[i], s[i], egp_id, season_form,
        sigma_0, xi_0,
        season_a, season_t, season_d,
        kappa1, kappa2, kappa3, kappa4
      );
    }
    exc[i] = m[i] > 0;
  }
  
  array[vc_dim, vc_dim] int vc_mat; // canonical vine copula matrix (upper right corner)
  for (i in 1:vc_dim) {
    for (j in 1:vc_dim) {
      if (j >= i)
        vc_mat[i, j] = i;
      else
        vc_mat[i, j] = 0;
    }
  }
  
  array[vc_dim, vc_dim] int cop_pos; // position in the parameter vector
  int d_paircop; // numbers of pair-copulas in the vine
  d_paircop = 0;
  for (i in 1:vc_dim) {
    for (j in 1:vc_dim) {
      if (j > i) {
        d_paircop = d_paircop + 1;
        cop_pos[i, j] = d_paircop;
      }
      else
        cop_pos[i, j] = 0;
    }
  }
  
  array[d_paircop] int npar_cop;
  for (i in 1:vc_dim) {
    for (j in 1:vc_dim) {
      if (cop_pos[i, j] > 0) {
        if (cop_id[i, j] > 10) // COULD NEED MODIF IF MORE COP FAMILY IMPLEMENTED
          npar_cop[cop_pos[i, j]] = 2;
        else
          npar_cop[cop_pos[i, j]] = 1;
      }
    }
  }
  
  
  int<lower=0> J = N %/% K;
  array[K, J * 4] real x_r;
  array[K, 1+4*vc_dim*vc_dim] int x_i;
  {
    int pos = 1;
    for (k in 1:K) {
      int end = pos + J - 1;
      
      matrix[J, 4] temp;
      temp[, 1] = m[pos:end];
      temp[, 2] = exc[pos:end];
      temp[, 3] = p_egp[pos:end];
      temp[, 4] = d_egp[pos:end];
      
      x_r[k, ] = to_array_1d(temp);
      
      x_i[k, 1] = vc_dim;
      x_i[k, 2:(1+vc_dim*vc_dim)] = to_array_1d(cop_id);
      x_i[k, (2+1*vc_dim*vc_dim):(1+2*vc_dim*vc_dim)] = to_array_1d(rotation);
      x_i[k, (2+2*vc_dim*vc_dim):(1+3*vc_dim*vc_dim)] = to_array_1d(vc_mat);
      x_i[k, (2+3*vc_dim*vc_dim):(1+4*vc_dim*vc_dim)] = to_array_1d(cop_pos);
      
      pos += J;
    }
  }
  
  real eps;
  eps = 1e-9;
}

parameters {
  real<lower=0+eps, upper=1-eps> lambda_a; // long clusters base intensity
  
  real<lower=0+eps, upper=0.1> alpha; // long clusters memory kernel
  
  real<lower=0+eps> lambda_b; // short clusters base intensity
  real<lower=0+eps> beta;
  real<lower=0.6+eps> gamma; // short clusters memory kernel
  
  real<lower=0+eps> delta1; // impact function
  real<lower=0+eps> delta2;

  real<lower=0+eps> omega1; // non-linear intensity function
  // real omega1; // for 2 param link
  real<lower=0+eps> omega2;
  real<lower=0+eps> omega3;
  
  vector<lower=0, upper=1>[d_paircop] theta1_g; // generic copula parameter, in row order of
  vector<lower=0, upper=1>[d_paircop] theta2_g; // the top-right corner of the vine copula matrix
  
  array[K] vector[0] dummy;
}

model {
  theta1_g[d_paircop] ~ uniform(0, 1);
  for (i in 1:d_paircop) {
    if (npar_cop[i] == 2)
      theta2_g[i] ~ uniform(0, 1);
  }
  
  lambda_a ~ gamma(2, 1);
  alpha ~ gamma(2, 60);
  
  lambda_b ~ gamma(2, 1);
  beta ~ gamma(2, 1);
  gamma ~ gamma(8, 2);
  
  delta1 ~ gamma(2, 1);
  delta2 ~ gamma(2, 1);
  
  omega1 ~ gamma(2, 1);
  omega2 ~ gamma(2, 1);
  omega3 ~ gamma(2, 1);
  
  vector[10 + 2*d_paircop] par; // vector of parameters to feed to map_rect()
  
  par[1] = lambda_a;
  par[2] = alpha;
  par[3] = lambda_b;
  par[4] = beta;
  par[5] = gamma;

  par[6] = delta1;
  par[7] = delta2;
  
  par[8] = omega1;
  par[9] = omega2;
  par[10] = omega3;
  
  par[(11 + 0*d_paircop):(10 + 1*d_paircop)] = theta1_g;
  par[(11 + 1*d_paircop):(10 + 2*d_paircop)] = theta2_g;
  
  target += sum(map_rect(hawkes_mr, par, dummy, x_r, x_i));
}

generated quantities {
  int vcdim;
  array[vc_dim, vc_dim] int copid;
  array[vc_dim, vc_dim] int rotat;
  array[vc_dim, vc_dim] int vcmat;
  array[vc_dim, vc_dim] int coppos;
  vcdim = vc_dim;
  copid = cop_id;
  rotat = rotation;
  vcmat = vc_mat;
  coppos = cop_pos;
  
  vector[d_paircop] theta1;
  vector[d_paircop] theta2;
  for (i in 1:vc_dim) {
    for (j in 1:vc_dim) {
      if (cop_pos[i, j] > 0) {
        vector[2] temp;
        temp = cop_par_map(cop_id[i, j], theta1_g[cop_pos[i, j]], theta2_g[cop_pos[i, j]]);
        theta1[cop_pos[i, j]] = temp[1];
        theta2[cop_pos[i, j]] = temp[2];
      }
    }
  }
}
