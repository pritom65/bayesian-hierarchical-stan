data {
    int< lower=0> N;
    vector[N] severity;
    array[N] int< lower=0> y;
}

parameters {
    real prms_alpha;
    real prms_beta;
    real < lower=0> phi;
}

model {
    prms_alpha ~ normal(0, 2.5);
    prms_beta ~ normal(0, 1);
    phi ~ exponential(1);
    y ~ neg_binomial_2_log(prms_alpha + prms_beta * severity, phi);
}

generated quantities {
   array[N] int< lower=0> y_rep;
   vector[N] log_lik;
   for(n in 1:N){
    real eta = prms_alpha + prms_beta * severity[n];
    y_rep[n] = neg_binomial_2_log_rng(eta, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(y[n]|eta, phi);
   }
}