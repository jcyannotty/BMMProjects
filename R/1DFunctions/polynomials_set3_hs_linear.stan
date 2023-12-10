data{
  int<lower=1> N; // number of observations
  int<lower=1> d; // number of input variables
  int<lower=0> d_discrete ; // number of discrete dummy inputs
  int<lower=2> K; // number of models

  matrix[N,d] X;
  
  // predictors
  // including continuous and discrete in dummy variables , no constant
  matrix[N,K] lpd_point; // the input pointwise predictive density
  real<lower=0> tau_mu;
  real<lower=0> tau_discrete; // global regularization for discrete x
  real<lower=0> tau_con; // overall regularization for continuous x
}

transformed data{
  matrix[N,K] exp_lpd_point = exp(lpd_point);
}

parameters{
  vector[K-1] mu;
  real mu_0;
  vector<lower=0>[K-1] sigma;
  vector<lower=0>[K-1] sigma_con;
  vector[d-d_discrete] beta_con[K-1];
}

transformed parameters {
  simplex[K] w[N];
  matrix[N,K] f;

  for(k in 1:(K-1))
    f[,k] = X*beta_con[k]*sigma_con[k] + (mu_0*tau_mu + mu[k]*tau_mu);
  f[,K] = rep_vector(0,N);
  for(n in 1:N)
    w[n] = softmax(to_vector(f[n,1:K]));
}


model{
  for(k in 1:(K-1)){
    beta_con[k] ~ std_normal();
  }
  mu ~ std_normal();
  mu_0 ~ std_normal();
  sigma_con ~ normal(0,tau_con);
  for(i in 1:N)
    target += log(exp_lpd_point[i,]*w[i]); // log likelihood
}

// optional block : needed if an extra layer of LOO ( eq .28) is called

generated quantities{
  vector[N] log_lik;
  for(i in 1:N)
    log_lik[i] = log(exp_lpd_point[i,]*w[i]);
}
