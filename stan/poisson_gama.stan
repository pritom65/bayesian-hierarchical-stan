data {
    int<lower=0> y;
    real<lower=0> t;
    real<lower=0> prms_alpha;
    real<lower=0> prms_beta;
}

parameters {
    real<lower=0> lambda;
}

model {
    lambda ~ gamma(prms_alpha,prms_beta);
    y ~ poisson(lambda*t);
}

generated quantities {
   int y_rep = poisson_rng(lambda*t);
   real rate_per_1000 = lambda*1000;
   int above_threshold = rate_per_1000 > 4;
}
