data {
    int<lower=0> N;
    vector[N] age;
    vector[N] sbp;
}

parameters {
    real prms_alpha;
    real prms_beta;
    real< lower=0> prms_sigma;
}

model {
    prms_alpha ~ normal(100,20);
    prms_beta ~ normal(0,5);
    prms_sigma ~ normal(0,10);

    sbp ~ normal(prms_alpha + prms_beta*age,prms_sigma);
}

generated quantities {
  vector[N] mu;
  vector[N] y_rep;
  vector[N] log_lik;
  mu = prms_alpha + prms_beta * age;
  
  for (n in 1:N) {
    y_rep[n]   = normal_rng(mu[n], prms_sigma);
    log_lik[n] = normal_lpdf(sbp[n] | mu[n], prms_sigma);
  }
}