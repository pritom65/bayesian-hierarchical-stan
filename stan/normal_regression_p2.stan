data {
    int<lower=0> N;
    vector[N] age;
    vector[N] bmi;
    vector[N] sbp;
}

parameters {
    real prms_alpha;
    real prms_beta_age;
    real prms_beta_bmi;
    real< lower=0> prms_sigma;
}

model {
    prms_alpha ~ normal(120,10);
    prms_beta_age ~ normal(0,5);
    prms_beta_bmi ~ normal(0,5);
    prms_sigma ~ normal(0,10);

    sbp ~ normal(prms_alpha + prms_beta_age*age + prms_beta_bmi*bmi,prms_sigma);
}

generated quantities {
  vector[N] mu;
  vector[N] y_rep;
  vector[N] log_lik;
  mu = prms_alpha + prms_beta_age*age + prms_beta_bmi*bmi;
  
  for (n in 1:N) {
    y_rep[n]   = normal_rng(mu[n], prms_sigma);
    log_lik[n] = normal_lpdf(sbp[n] | mu[n], prms_sigma);
  }
}