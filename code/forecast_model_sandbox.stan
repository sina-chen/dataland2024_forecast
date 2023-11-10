data {
  
  int<lower=1> N;                           // number of poll results (one for each poll and party)
  int<lower=1> K;                           // number of unique parties
  int<lower=1> R;                           // number of unique elections
  int<lower=1> KR;                          // number of unique party-election combinations

  int<lower=1> P;                           // number of polls
  
  vector<lower=0, upper=1>[N] poll;         // poll estimate for each poll and party
  vector<lower=0, upper=1>[KR] prev_vote;   // previous election results for each party (starting value for random walk)

  vector<lower=0>[P] n;                     // sample size
  int<lower=1> max_D;                       // max. distance between poll and election for random walk length

  int<lower=1, upper=K> k_id[N];            // party identifier
  int<lower=1, upper=R> r_id[P];            // election id for each poll
  int<lower=1, upper = max_D> t_id[P];      // days to election id for each poll
  int<lower=1, upper=P> p_id[N];            // poll id for each poll party support
  vector<lower=0, upper=R>[R] st_id;        // Synapse Territories id
  
  vector[K-1] zero_r[R];                    // zero mean matrix for mvn distribution of party-election bias
  vector[K] zero_rK[R];                     // zero mean matrix for mvn distribution

}

transformed data {
  
  // rescale previous election result to log scale for random walk
  vector[K] log_prev_vote[R];
  for(r in 1:R) {
    if(st_id[r] == 0) {
      log_prev_vote[r, 1:(K-1)] = log(prev_vote[((r-1)*K+1):((r*K)-1)]);
      log_prev_vote[r, K] = 0;
    } else {
      log_prev_vote[r] = log(prev_vote[((r-1)*K+1):(r*K)]);
    }
  }

}

parameters {
  
  // Hyper parameters
  
  vector[K-1] alpha1_raw[R];                     // party-election bias for K-1 parties due to redundant parametrisation
  cholesky_factor_corr[K-1] lkj_alpha1_corr;     // cholesky factor of corr. matrix 
  real<lower=0> sd_tau_alpha1;                   // sd for non-centerred common party sd
  vector<lower=0>[K-1] tau_alpha1_sc;            // scale for non-centerred common party sd
  
  vector[K] theta_raw[R, (max_D+1)];             // unscaled party support for every day and every election
  cholesky_factor_corr[K] lkj_theta_corr[R];     // cholesky factor of corr. matrix for every election
  real<lower=0> sd_tau_theta[R];                 // sd for non-centerred common party sd for every election
  vector<lower=0>[K] tau_theta_sc[R];            // scale for non-centerred common party sd for every election
  
  vector[K] phi[R];                              // party-election excess variance
  cholesky_factor_corr[K] lkj_phi_corr;          // cholesky factor of corr. matrix 
  real<lower=0> sd_tau_phi;                      // sd for non-centerred common party sd
  vector<lower=0>[K] tau_phi_sc;                 // scale for non-centerred common party sd
  
}

transformed parameters {
  
  // parameters for mean equation
  vector[K] alpha1[R];
  vector[K] theta[R, (max_D+1)];         

  // party specific sd for mvn distribution
  vector<lower=0>[K-1] tau_alpha1; 
  vector<lower=0>[K] tau_theta[R]; 
  vector<lower=0>[K] tau_phi;   

  
  // non-centered parametrisation for party specific sd
  tau_alpha1 = tau_alpha1_sc * sd_tau_alpha1;
  tau_phi = tau_phi_sc * sd_tau_phi;
  for(r in 1:R) {
    if(st_id[r] == 0){ // elections outside Synapse Territories with K-1 parties
      tau_theta[r, 1:(K-1)] = tau_theta_sc[r, 1:(K-1)] * sd_tau_theta[r];
      tau_theta[r, K] = 1;
    } else { // elections inside Synapse Territories with K parties
      tau_theta[r, :] = tau_theta_sc[r, :] * sd_tau_theta[r];
    }
  }
  
  // softmax transformation of party support to ensure sum to 1 for every day
  for(r in 1:R) {
    for(d in 1:max_D+1) {
      if(st_id[r] == 0){
        theta[r, d, 1:(K-1)] = softmax(theta_raw[r, d, 1:(K-1)]);
        theta[r, d, K] = 0;
      } else {
      theta[r, d, :] = softmax(theta_raw[r, d, :]);
      }
    }
  }
  
  
  // redundant parametrisation: add 0 for k-th party such that softmax transformation is identifiable
  for(r in 1:R) {alpha1[r] = append_row(alpha1_raw[r],0);}  

}

model {

  vector[N] log_mu;       // unscaled mean on log scale
  vector[N] mu;           // scaled mean on original scale
  vector[N] log_sigma2;   // unscaled variance on log scale
  vector[N] sigma2;       // variance

  // hyper-priors: random intercept SD
  sd_tau_alpha1 ~ normal(0, 0.05)T[0,];
  sd_tau_phi ~ normal(0, 0.05)T[0,];
  
  // hyper-priors: random intercept scale factor, half normal distribution to ensure that it will be positive
  for(i in 1:(K-1)) {
    tau_alpha1_sc[i] ~ normal(0, 1)T[0,];
  }
  for(i in 1:K) {tau_phi_sc[i] ~ normal(0, 1)T[0,];}
  
  // cholesky factor of correlation matrix
  lkj_alpha1_corr ~ lkj_corr_cholesky(2);
  for(r in 1:R) {
    sd_tau_theta[r] ~ normal(0, 0.05)T[0,]; // random intercept SD
    lkj_theta_corr[r] ~ lkj_corr_cholesky(2);
  }
  lkj_phi_corr ~ lkj_corr_cholesky(2);
  
  // multivariate normal distribution 
  alpha1_raw ~ multi_normal_cholesky(zero_r, diag_pre_multiply(tau_alpha1, lkj_alpha1_corr)); // for K-1 parameter estimates due to redundant parametrization
  phi ~ multi_normal_cholesky(zero_rK, diag_pre_multiply(tau_phi, lkj_phi_corr)); 
  
  // multivariate normal random walk
  for(r in 1:R) {
    if(st_id[r] == 0){ // elections outside Synapse Territories with K-1 parties
      for(k in 1:(K-1)){
        tau_theta_sc[r, k] ~ normal(0, 1)T[0,];
      }
      theta_raw[r, 1, 1:(K-1)] ~ multi_normal_cholesky(log_prev_vote[r, 1:(K-1)], diag_pre_multiply(tau_theta[r, 1:(K-1)], lkj_theta_corr[r, 1:(K-1), 1:(K-1)]));  // set initial time step (col) for all groups (rows)
      for(d in 2:max_D+1){
        theta_raw[r, d, 1:(K-1)] ~ multi_normal_cholesky(theta_raw[r, d-1, 1:(K-1)], diag_pre_multiply(tau_theta[r, 1:(K-1)], lkj_theta_corr[r, 1:(K-1), 1:(K-1)]));
      }
    } else { // elections inside Synapse Territories with K parties
      for(k in 1:K){
        tau_theta_sc[r, k] ~ normal(0, 1)T[0,];
      }
      theta_raw[r, 1, :] ~ multi_normal_cholesky(log_prev_vote[r, ], diag_pre_multiply(tau_theta[r, ], lkj_theta_corr[r, ]));  // set initial time step (col) for all groups (rows)
      for(d in 2:max_D+1){
        theta_raw[r, d, :] ~ multi_normal_cholesky(theta_raw[r, d-1, :], diag_pre_multiply(tau_theta[r, ], lkj_theta_corr[r, ]));
      }
    }
  }
  
  // unscaled mean
  for(i in 1:P) {
    if(st_id[r_id[i]] == 0) { // elections outside Synapse Territories with K-1 parties
      log_mu[((i-1)*K+1):((i*K)-1)] = alpha1[r_id[i], 1:(K-1)] +
      log(theta[r_id[i], t_id[i], 1:(K-1)]); 
      log_mu[i*K] = 0;
    } else if(st_id[r_id[i]] == 1) { // elections inside Synapse Territories with K parties
      log_mu[((i-1)*K+1):(i*K)] = alpha1[r_id[i], :] +
      log(theta[r_id[i], t_id[i], :]); 

    }
  }
  
  // scaled mean on original scale to ensure sum to 1 constrain
  for(i in 1:P) {
    if(st_id[r_id[i]] == 0) {
      mu[((i-1)*K+1):((i*K)-1)] = softmax(log_mu[((i-1)*K+1):((i*K)-1)]); 
      mu[i*K] = 0;
    } else if(st_id[r_id[i]] == 1) {
      mu[((i-1)*K+1):(i*K)] = softmax(log_mu[((i-1)*K+1):(i*K)]); 

    }
  }

  // variance 
  for(i in 1:P) {
    if(st_id[r_id[i]] == 0) {
      log_sigma2[((i-1)*K+1):((i*K)-1)] = log(mu[((i-1)*K+1):((i*K)-1)].*(1-mu[((i-1)*K+1):((i*K)-1)])/n[i]) + 
      phi[r_id[p_id[i]], k_id[i]];
      log_sigma2[i*K] = 0; 
      sigma2[((i-1)*K+1):((i*K)-1)] = exp(log_sigma2[((i-1)*K+1):((i*K)-1)]); 
      sigma2[i*K] = 0.00000001; // set sigma2 for SSP outside Synapse Territories to very small value due to >0 constrain
    } else if(st_id[r_id[i]] == 1) {
      log_sigma2[((i-1)*K+1):(i*K)] = log(mu[((i-1)*K+1):((i*K))].*(1-mu[((i-1)*K+1):((i*K))])/n[i]) + 
      phi[r_id[p_id[i]], k_id[i]]; 
      sigma2[((i-1)*K+1):(i*K)] = exp(log_sigma2[((i-1)*K+1):(i*K)]);
    }
  }
  
  // poll estimate
  for(i in 1:N) {
     if(st_id[r_id[p_id[i]]] == 1 || (st_id[r_id[p_id[i]]] == 0 && k_id[i] < 4)) {
      poll[i] ~ normal(mu[i], sqrt(sigma2[i]));
    }
  }
  
}
