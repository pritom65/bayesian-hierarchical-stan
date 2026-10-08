data {
    int< lower=0> N;
    int< lower=0> J;
    array[N] int< lower=0, upper=J> site;
    vector[N] y;
}

parameters {
    vector[J] theta;
    real< lower=0> sigma;
}

model {
    theta ~ normal(0, 100);
    sigma ~ normal(0, 100);
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