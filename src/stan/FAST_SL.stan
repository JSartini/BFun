data {
  int N;             // Total number of functions
  int I;             // Number of groups
  array[N] int ID;   // Index-based ID per function
  int M;             // Number of observations along the domain

  int Q;    // Dimension of spline bases
  int P;    // Number of covariates
  int K1;   // Number of Eigenfunctions at group level

  matrix[N, M] Y;     // Original data
  matrix[N, P] X;     // Matrix of covariates
  matrix[M, Q] B;     // Orthogonalized basis
  matrix[Q, Q] P_alpha;     // Penalty matrix for splines
}

transformed data {
  real EV1 = eigenvalues_sym(P_alpha)[1]; // Leading eigenvalue
}

parameters {
  real<lower=0> sigma2; // Error in observation

  // Fixed-effect components
  matrix[Q,P] w_mu;          // B-spline weights for overall mean
  vector<lower=0>[P] H_mu;   // Smoothness parameters for FE

  // Subject-specific components
  vector<lower=0>[K1] lambda_1;  // Variances of each eigen-function
  vector<lower=0>[K1] H_1;        // Smoothing multipliers
  matrix[Q, K1] X_1;              // Full-rank representation decomposed into group EF
  matrix[I, K1] Scores1;          // Raw group scores
}

transformed parameters{
  // Fixed effects
  matrix[N, M] FE = X * (B * w_mu)';

  // Polar decomposition
  matrix[Q, K1] Psi_1;
  {
    vector[K1] trans_eval_V1;   // Inverse sqrt of eigenvalues of X_1'*X_1
    matrix[K1,K1] evec_V1;      // Eigenvectors of X_1'*X_1

    trans_eval_V1 = 1/sqrt(eigenvalues_sym(X_1'*X_1));
    evec_V1 = eigenvectors_sym(X_1'*X_1);
    Psi_1 = X_1*evec_V1*diag_matrix(trans_eval_V1)*evec_V1';
  }
}

model {
  // Smoothing weight priors
  H_mu ~ gamma(0.001, 0.001);
  H_1 ~ gamma(0.01, EV1/2 + 0.01);

  // Variance component priors
  lambda_1 ~ inv_gamma(0.001, 0.001);
  sigma2 ~ inv_gamma(0.001, 0.001);

  // Smoothing additions to the target density
  for(p in 1:P){
    target += Q/2.0 * log(H_mu[p]) - H_mu[p] / 2.0 * w_mu[,p]' * P_alpha * w_mu[,p];
  }
  for(i in 1:K1){
    target += Q/2.0 * log(H_1[i]) - H_1[i] / 2.0 * Psi_1[,i]' * P_alpha * Psi_1[,i];
  }

  // Uniform priors through matrix normals
  to_vector(X_1) ~ normal(0, 1);

  // Score priors
  for(k in 1:K1){
    to_vector(Scores1[,k]) ~ normal(0, sqrt(lambda_1[k]));
  }

  // Likelihood
  {
    matrix[I, M] UDV_1 = Scores1 * (B * Psi_1)';
    to_vector(Y) ~ normal(to_vector(FE + UDV_1[ID]), sqrt(sigma2));
  }
}
