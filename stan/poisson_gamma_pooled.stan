data {
    int<lower=0> J;
    array[J] int<lower=0> y;
    vector<lower=0>[J] t;
    real<lower=0> prms_alpha;
    real<lower=0> prms_beta;
}

parameters {
    real lambda;
}

model {
    lambda ~ gamma(prms_alpha, prms_beta);
    y ~ poisson(lambda*t);
}


generated quantities {
   array[J] int y_rep;
   real rate_per_1000 = lambda * 1000;
   for(j in 1:J){
    y_rep[j] = poisson_rng(lambda*t[j]);
   }
}