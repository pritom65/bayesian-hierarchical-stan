data {
    int< lower=0> N;
    int< lower=0> J;
    array[N] int< lower=0, upper=J> site;
    vector[N] y;
}

parameters {
    real mu;
    real< lower=0> sigma;
}

model {
    mu ~ normal(0, 10);
    sigma ~ normal(0, 10);
    y ~ normal(mu, sigma);
}

generated quantities {
    vector[N] log_lik;
    vector[N] y_rep;

    for(n in 1:N){
        log_lik[n] = normal_lpdf(y[n]|mu, sigma);
        y_rep[n] = normal_rng(mu, sigma);
    }
}