data {
    int< lower=0> N;
    int< lower=0> J;
    array[N] int< lower=0, upper=J> site;
    vector[N] y;
}

parameters {
    real mu;
    real< lower=0> tau;
    real< lower=0> sigma;
    vector[J] theta_tilde;
}

transformed parameters {
   vector[J] theta = mu + tau * theta_tilde;
}

model {
    mu ~ normal(0, 10);
    tau ~ normal(0, 10);
    sigma ~ normal(0, 10);
    theta_tilde ~ normal(0, 10);

    y ~ normal(theta[site], sigma);
}

generated quantities {
   vector[N] log_lik;
   vector[N] y_rep;
   for(n in 1:N){
    log_lik[n] = normal_lpdf(y[n]| theta[site[n]], sigma);
    y_rep[n] = normal_rng(theta[site[n]], sigma);
   }
}