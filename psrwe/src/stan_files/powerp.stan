//
//  TYPICAL POWER PRIOR FOR CONTINUOUS DATA
//  Power parameter A is fixed
//

data {
  //existing data
  int<lower = 1>  N0;
  real            YBAR0;
  real<lower = 0> SD0;

  //current data
  int<lower = 1>  TN1;
  real            Y1[TN1];

  //target borrowing
  real<lower = 0> A;
}

transformed data {
  real<lower = 0> a0;
  real<lower = 0> sn0;

  a0  = 1 < A/N0 ? 1 : A/N0;
  sn0 = SD0/sqrt(N0*1.0);
}

parameters {
  real          theta;
  real<lower=0> tau1;
}

model {
  //prior
  theta ~ normal(0, 1000);
  tau1  ~ cauchy(0, 2.5);

  //likelihood
  if (N0 > 0) {
    target +=  normal_lpdf(YBAR0 | theta, sn0) * a0;
  }

  Y1 ~ normal(theta, tau1);
}
