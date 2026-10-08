data {
    int N;
    vector[N] severity;
    array[N] int<lower=0> y;
}

parameters {
    real prms_alpha;
    real prms_beta;
}

model {
    prms_alpha ~ normal(0, 2.5);
    prms_beta  ~ normal(0, 2.5);
    y ~ poisson_log(prms_alpha + prms_beta * severity);
}

generated quantities {
   vector[N] mu;
   array[N] int y_rep;
   vector[N] log_lik;

   for(n in 1:N){
    mu[n] = exp(prms_alpha + prms_beta * severity[n]);
    y_rep[n] = poisson_rng(mu[n]);
    log_lik[n] = poisson_lpmf(y[n]| mu[n]);
   }
}