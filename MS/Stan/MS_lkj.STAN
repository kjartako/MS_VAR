data {
int<lower=1> T; // Number of observations
int<lower=1> dim; // Dimension of observations
matrix[T,dim] y; // Observations
int<lower=1> nreg; // Number of regimes
int<lower = 0, upper = 1> mean_reg;
int<lower = 0, upper = 1> sigma_reg;
vector[nreg] Q_alpha; // Transition probabilities hyperparameter
real mu_mean; // Mean value of the normal prior for the constant means
real mu_sd; // Standard deviation of the normal prior for the constant means
real eta; // LKJ hyperparameter
real gamma_alpha; // Inverse gamma hyperparameter (shape)
real gamma_beta;  // Inverse gamma hyperparameter (scale)
}

parameters {
simplex[nreg] Q[nreg] ;  // Transition probabilities
vector[dim] mu[mean_reg ? nreg : 1]; // Mean

vector<lower=0>[dim] vari[sigma_reg ? nreg : 1]; // Variance
cholesky_factor_corr[dim] L_corr[sigma_reg ? nreg : 1]; //Lower Cholesky factors of the correlation matrices
}

transformed parameters{
vector[dim] sdev[sigma_reg ? nreg : 1];
matrix[dim,dim] L_sigma[sigma_reg ? nreg : 1];

for(i in 1:(sigma_reg ? nreg : 1)){
	for(j in 1:dim){
		sdev[i,j]=sqrt(vari[i,j]);
	}
	L_sigma[i] = diag_pre_multiply(sdev[i], L_corr[i]);
}
}

model {
// Priors

for (k in 1:nreg){
	Q[k] ~ dirichlet(Q_alpha);
}

for (k in 1:(mean_reg ? nreg : 1)){
	mu[k] ~ normal(mu_mean,mu_sd);	
}

for (k in 1:(sigma_reg ? nreg : 1)){
	vari[k] ~ inv_gamma(gamma_alpha, gamma_beta);
	L_corr[k] ~ lkj_corr_cholesky(eta);
}
{
    //forward algorithm for computing log p(y|...)
    
	real fwd[nreg,nreg];
    real alphas[T,nreg];
	real py[nreg];
	real logQ[nreg,nreg];
	int n_m = mean_reg ? nreg : 1;
	int n_s = sigma_reg ? nreg : 1;
	
	for(i in 1:nreg){
		for (j in 1:nreg){
			logQ[i,j]= log(Q[i,j]);
		}
	}
		
	for(i in 1:nreg){	
		alphas[1,i] = multi_normal_cholesky_lpdf(y[1,] | mu[min(n_m,i)],L_sigma[min(n_s,i)]);
	}
	
    for (t in 2:T){	

		for(i in 1:nreg){		
			py[i] = multi_normal_cholesky_lpdf(y[t,] | mu[min(n_m,i)],L_sigma[min(n_s,i)]);
		}	
		
		for(i in 1:nreg){
		  for (j in 1:nreg){
			fwd[i,j] = alphas[t-1,j] + logQ[i,j] + py[i];
		}
		  alphas[t,i]=log_sum_exp(fwd[i,]);
		}
		  
    }
	
	// adding the marginal log-likelihood to Stan target distribution
	
    target += log_sum_exp(alphas[T,]);	
}
}

generated quantities {
	matrix[dim,dim] sigma[sigma_reg ? nreg : 1];
	int<lower=1,upper=nreg> S[T];
	real log_p_S;
  
	for(i in 1:(sigma_reg ? nreg : 1)){
		sigma[i]= multiply_lower_tri_self_transpose(L_sigma[i]);
	}
	
   // Viterbi algorithm 
  { 
    int back_ptr[T,nreg];
    real best_logp[T+1,nreg];
    real best_total_logp;
	real py[nreg];
	real logQ[nreg,nreg];
	int n_m = mean_reg ? nreg : 1;
	int n_s = sigma_reg ? nreg : 1;
    
    for(i in 1:nreg){
	best_logp[1,i]=0;
	}
	
	for(i in 1:nreg){
		for (j in 1:nreg){
			logQ[i,j]= log(Q[i,j]);
	}
	}
	
    for (t in 1:T) {
				  	  
	    for(i in 1:nreg){
		  py[i] = multi_normal_cholesky_lpdf(y[t,] | mu[min(n_m,i)],L_sigma[min(n_s,i)]);
		  best_logp[t+1,i] = negative_infinity();
	    }	
	  
        for (j in 1:nreg) {
			real logp[nreg];

			for(k in 1:nreg){	
				logp[k] = best_logp[t,j] + logQ[j,k] + py[k];
					
				if (logp[k] > best_logp[t+1,k]){
					back_ptr[t,k] = j;
					best_logp[t+1,k] = logp[k];
				}	
			}
        }
    }
	
    log_p_S = max(best_logp[T+1]);
    for (k in 1:nreg){
      if (best_logp[T+1,k] == log_p_S){
        S[T] = k;
	  }
	}
    for (t in 1:(T-1)){
      S[T-t] = back_ptr[T-t+1, S[T-t+1]];
	}
  }
}
