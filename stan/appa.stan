data {
  int<lower=0> N;                  // number of observations
  int<lower=1> R;                  // number of brain regions
  int<lower=1> J;                  // number of slides
  vector[N] y;                     // log-transformed signal
  int<lower=0,upper=1> is_wt[N];  // genotype indicator (0 = KO, 1 = WT)
  int<lower=1,upper=J> slide[N];  // slide index
  int<lower=1,upper=R> region[N]; // region index
}

parameters {
  vector[R] mu_region;            // region baseline
  vector[J] alpha_slide;          // slide random effects
  vector[R] delta_region;         // region-specific genotype effects
  real<lower=0> sigma;            // residual SD
}

model {
  mu_region ~ normal(0, 2);
  alpha_slide ~ normal(0, 1);
  delta_region ~ normal(0, 1);
  sigma ~ exponential(1);

  for (n in 1:N) {
    y[n] ~ normal(
      mu_region[region[n]] +
      alpha_slide[slide[n]] +
      delta_region[region[n]] * is_wt[n],
      sigma
    );
  }
}
