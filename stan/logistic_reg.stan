data{
    int<lower=0> N;
    vector[N] apache;
    array[N] int<lower=0,upper=1> y;
    real sigma_alpha;
    real sigma_beta;
}

parameters {
    real prms_alpha;
    real prms_beta;
}

model {
    prms_alpha ~ normal(0,sigma_alpha);
    prms_beta ~  normal(0,sigma_beta);

    y ~ bernoulli_logit(prms_alpha + prms_beta * apache);
}

generated quantities {
   vector[N] p;
   vector[N] log_lik;
   vector[N] y_rep;

   for(n in 1:N){
    p[n] = inv_logit(prms_alpha + prms_beta * apache[n]);
    log_lik[n] = bernoulli_logit_lpmf(y[n]|prms_alpha + prms_beta * apache[n]);
    y_rep[n] = bernoulli_rng(p[n]);
   }
}