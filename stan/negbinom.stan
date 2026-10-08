data {
    int N;
    vector[N] severity;
    array[N] int< lower=0> y;
}

parameters {
    real prms_alpha;
    real prms_beta;
    real< lower=0> phi;
}

model {
    prms_alpha ~ normal(0, 2.5);
    prms_beta  ~ normal(0, 2.5);
    phi ~ exponential(1);

    y ~ neg_binomial_2_log(prms_alpha + prms_beta * severity, phi);
}

generated quantities {
   vector[N] mu;
   array[N] int y_rep;
   vector[N] log_lik;

   for(n in 1:N){
    mu[n] = exp(prms_alpha + prms_beta * severity[n]);
    y_rep[n] = neg_binomial_2_rng(mu[n], phi);
    log_lik[n] = neg_binomial_2_lpmf(y[n]|mu[n], phi);
   }
}