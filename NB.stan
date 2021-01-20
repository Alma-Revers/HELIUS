data {
  int<lower=0> N_subj;            //Total number of samples/subjects
  int<lower=0> N_otu;             //Total number of otus
  int<lower=0> N_family;          //Total number of families
  int<lower=0> N_variables;       //Total number of variables + intercept
  int<lower=0> y[N_subj,N_otu];   //Otu count matrix
  matrix[N_subj, N_variables] x;  //Design matrix with the food item score
  vector[N_subj] offset;          //Log of the total counts per subject
  int index_family[N_otu];        //Per OTU the index of its phylogenetic family
}
parameters {
  vector<lower=0>[N_variables] sigma_beta;
  vector[N_variables] mu_beta;
  cholesky_factor_corr[N_variables] L_Omega;
 
  matrix<lower=0>[N_variables, N_family] sigma_beta_f;
  matrix[N_variables, N_family] mu_beta_f;
 
  matrix[N_variables, N_otu] beta_otu;
  
  vector<lower=0>[N_otu] phi;
  real<lower=0> sigma_phi;
  real a1;
  real alpha_0;
}

transformed parameters {
  matrix[N_subj, N_otu] mu;
  vector[N_otu] mu_phi;
  matrix[N_variables, N_variables] L_Sigma[N_family];
 
  for ( f in 1:N_family){
    L_Sigma[f] = diag_pre_multiply(sigma_beta_f[,f], L_Omega);
  }
 
  for ( j in 1:N_otu){
    mu[,j] = (x * beta_otu[,j]) + offset;
    mu_phi[j] = a1/mean(exp(mu[,j])) + alpha_0;
  }
}

model {
  sigma_beta ~ exponential(1);
  mu_beta ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(2);
 
  for ( f in 1:N_family){
    sigma_beta_f[,f] ~ exponential(1);
    mu_beta_f[,f] ~ normal(mu_beta, sigma_beta);
  }
  
  a1 ~ std_normal();
  alpha_0 ~ std_normal();
  sigma_phi ~ exponential(1);
 
  for(j in 1:N_otu){
    beta_otu[,j] ~ multi_normal_cholesky(mu_beta_f[,(index_family[j])], L_Sigma[(index_family[j])]);
    phi[j] ~ lognormal(mu_phi[j], sigma_phi);
    y[,j] ~ neg_binomial_2_log(mu[,j], phi[j]);
  }
} 
