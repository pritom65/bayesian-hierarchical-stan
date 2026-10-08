data {
    int<lower=0> n;
    array[n] int<lower=0, upper=1> y;
    real< lower=0 > prms_alpha;
    real< lower=0 > prms_beta;
}

parameters{
    real< lower=0, upper=1> theta;
}

model {
    theta ~ beta(prms_alpha,prms_beta);
    y ~bernoulli(theta);
}


generated quantities {
   int y_rep = binomial_rng(n, theta);
   int above_threshold = theta > .3;
}

