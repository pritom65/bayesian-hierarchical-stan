data {
    int< lower=0> N;
    int< lower=0> J;
    array[N] int< lower=0, upper=J> hospital;
    vector[N] severity;
    array[N] int< lower=0, upper=1> y;
}

parameters {
    real mu_alpha;
    real< lower=0> tau;
    real prms_beta;
    vector[J] alpha_tilde;
}

transformed parameters {
   vector[J] alpha = mu_alpha + tau * alpha_tilde;
}

model {
    mu_alpha ~ normal(-1.5,1);
    tau ~ normal(0, 0.5);
    prms_beta ~ normal(0, 2.5);
    alpha_tilde ~ normal(0, 1);

    y ~ bernoulli_logit(alpha[hospital] + prms_beta * severity);
}

generated quantities {
   vector[J] p_hospital;
   vector[N] log_lik;
   array[N] int y_rep;

   for(j in 1:J){
    p_hospital[j] = inv_logit(alpha[j]);
   }

   for(n in 1:N){
    real eta = alpha[hospital[n]] + prms_beta * severity[n];
    log_lik[n] = bernoulli_logit_lpmf(y[n]|eta);
    y_rep[n] = bernoulli_rng( inv_logit(eta));
   }
}