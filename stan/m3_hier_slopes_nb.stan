data {
    int< lower=0> N;
    int< lower=0> J;
    array[N] int< lower=0, upper=J> hospital;
    array[N] int< lower=0> y;
    vector[N] severity;
}

parameters {
    real mu_alpha;
    real mu_beta;
    vector< lower=0>[2] tau;
    real< lower=0> phi;
    cholesky_factor_corr[2] L_Omega;
    matrix[2,J] z;
}

transformed parameters {
   matrix[2, J] u = diag_matrix(tau) * L_Omega * z;
   vector[J] alpha = mu_alpha + u[1,]';
   vector[J] beta_j = mu_beta + u[2,]';
}

model {
    mu_alpha ~ normal(0, 2.5);
    mu_beta ~ normal(0, 1);
    tau ~ normal(0, 1);
    phi ~ exponential(1);
    L_Omega ~ lkj_corr_cholesky(2);
    to_vector(z) ~ normal(0, 1);

    y ~ neg_binomial_2_log(alpha[hospital] + beta_j[hospital] .* severity, phi);
}

generated quantities {
  corr_matrix[2] Omega = multiply_lower_tri_self_transpose(L_Omega);
  real rho             = Omega[1, 2];
  array[N] int y_rep;
  vector[N] log_lik;
  for (n in 1:N) {
    real eta   = alpha[hospital[n]] + beta_j[hospital[n]] * severity[n];
    y_rep[n]   = neg_binomial_2_log_rng(eta, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta, phi);
  }
}