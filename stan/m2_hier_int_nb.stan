data {
    int< lower=0> N;
    int< lower=0> J;
    array[N] int< lower=0, upper=J> hospital;
    array[N] int< lower=0> y;
    vector[N] severity;
}

parameters {
    real mu_alpha;
    real< lower=0> tau;
    real prms_beta;
    real< lower=0> phi;
    vector[J] alpha_tilde;
}

transformed parameters {
   vector[J] alpha = mu_alpha + tau * alpha_tilde;
}

model {
    mu_alpha ~ normal(0, 2.5);
    tau ~ normal(0, 1);
    prms_beta ~ normal(0, 1);
    phi ~ exponential(1);
    alpha_tilde ~ normal(0, 1);

    y ~ neg_binomial_2_log(alpha[hospital] + prms_beta * severity, phi);
}

generated quantities {
   array[N] int y_rep;
   vector[N] log_lik;
   for(n in 1:N){
    real eta = alpha[hospital[n]] + prms_beta*severity[n];
    y_rep[n] = neg_binomial_2_log_rng(eta, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(y[n]|eta, phi);
   }
}