data {
  int N;         // Total number of functions
  int M;         // Number of observations along the domain

  int Q;    // Dimension of spline bases
  int P;    // Number of covariates

  matrix[N, M] Y;     // Original data
  matrix[N, P] X;     // Matrix of covariates
  matrix[M, Q] B;     // Orthogonalized basis
  matrix[Q, Q] P_alpha;     // Penalty matrix for splines
}

parameters {
  real<lower=0> sigma2; // Error in observation

  // Fixed-effect components
  matrix[Q,P] w_mu;          // B-spline weights for overall mean
  vector<lower=0>[P] H_mu;   // Smoothness parameters for FE
}

transformed parameters{
  // Fixed effects
  matrix[N, M] FE = X * (B * w_mu)';
}

model {
  // Smoothing weight priors
  H_mu ~ gamma(0.001, 0.001);

  // Variance component priors
  sigma2 ~ inv_gamma(0.001, 0.001);

  // Smoothing additions to the target density
  for(p in 1:P){
    target += Q/2.0 * log(H_mu[p]) - H_mu[p] / 2.0 * w_mu[,p]' * P_alpha * w_mu[,p];
  }

  // Likelihood
  {
    to_vector(Y) ~ normal(to_vector(FE), sqrt(sigma2));
  }
}
