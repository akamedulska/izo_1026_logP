functions{
  real hplcmodel(real fi, real logkw, real logka, real logS2A){
    
    real logk;					// retention factor
    real S1;						// slope coefficient
    
    S1 = (logkw - logka)*(1+10^logS2A);
    logk = logkw - S1 * fi / (1 + 10^logS2A * fi);
    
    return logk;
  }
}

data{
  int nAnalytes;			  // number of analytes
  int nObs;				    	// number of observations
  int analyte[nObs];		// analytes indexes
  vector[nObs] logkObs;	// observed retention factors
  vector[nObs] fi;			// organic modifier content in the mobile phase
  real logP[nAnalytes];         // molecular descriptor      
  int idxmissing[nAnalytes];    // indexes of logP missing value
  real<lower=0, upper=1> rawlambdai[nAnalytes];  // ratio of ionized to total concentration of a compound
  
  int nfiplot;								  // number of fi for plotting
  vector[nfiplot] fiplot;				// organic modifier content in the mobile phase
  
  int<lower = 0, upper = 1> run_estimation; // 0 for prior predictive, 1 for estimation
}


transformed data {
  real oddsi[nAnalytes];
 
for(j in 1:nAnalytes){
  oddsi[j] = logit(rawlambdai[j]);
}
}

parameters{
  ordered[2] logkwHat;					// mean value of logkw in population
  ordered[2] logkaHat;          // mean value of logka in population
  vector[2] logS2AHat;					// mean curvature coefficient for acetonitrile in population

  real<lower = 0> sigma;					// standard deviation for residuals
  vector<lower = 0>[3] omega1;		// diagonal elements of variance-covariance matrix 
  corr_matrix[3] rho1;					  // correlation matrix
  vector<lower = 0>[3] omega2;	  // diagonal elements of variance-covariance matrix 
  corr_matrix[3] rho2;						// correlation matrix
  real beta[6];									  // coefficient value for molecular descriptor
  real logPmissing[nAnalytes];    // logP Missing values
  vector[3] param1[nAnalytes];
  vector[3] param2[nAnalytes];
  real eta[nAnalytes];
}

transformed parameters{
 
  vector[3] miu1[nAnalytes];	 
  vector[3] miu2[nAnalytes];	 
  real logka[nAnalytes];     // mean value of logka in population
  real logkw[nAnalytes];     // mean value of logkw in population
  real logS2A[nAnalytes];    // mean curvature coefficient for acetonitrile in population
  cov_matrix[3] Omega1;			 // variance-covariance matrix
  cov_matrix[3] Omega2;			 // variance-covariance matrix
  vector[nObs] logkHat;			 // mean value of logk in population
  real<lower=0, upper=1> frac[nAnalytes];

  Omega1 = quad_form_diag(rho1, omega1);	// diag_matrix(omega) * rho * diag_matrix(omega)
  Omega2 = quad_form_diag(rho2, omega2);	// diag_matrix(omega) * rho * diag_matrix(omega)

  for(j in 1:nAnalytes){
    miu1[j,1]  = logkwHat[1]  +  beta[1] * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
    miu2[j,1]  = logkwHat[2]  +  beta[2] * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
    miu1[j,2]  = logkaHat[1]  +  beta[3] * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]); 
    miu2[j,2]  = logkaHat[2]  +  beta[4] * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]); 
    miu1[j,3]  = logS2AHat[1] +  beta[5] * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
    miu2[j,3]  = logS2AHat[2] +  beta[6] * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
    frac[j]  = inv_logit(oddsi[j] + eta[j]);
  }

	for(j in 1:nAnalytes){
	logkw[j]  = param1[j, 1]*frac[j] + param2[j, 1]*(1-frac[j]);
	logka[j]  = param1[j, 2]*frac[j] + param2[j, 2]*(1-frac[j]);
	logS2A[j] = param1[j, 3]*frac[j] + param2[j, 3]*(1-frac[j]);
	}
  

  for(i in 1:nObs){
    logkHat[i] = hplcmodel(fi[i], logkw[analyte[i]], logka[analyte[i]], logS2A[analyte[i]]);
 }
  
}
model{
  
  logkwHat[1]  ~ normal(1.054, 1.136);
  logkwHat[2]  ~ normal(2.053, 1.487);
  logkaHat[1]  ~ normal(-3.437, 1.062);
  logkaHat[2]  ~ normal(-1.885, 1.006);
  logS2AHat[1] ~ normal(log10(2), 0.2);
  logS2AHat[2] ~ normal(log10(2), 0.2);

  beta[1] ~ normal(0.7,0.25);
  beta[2] ~ normal(0.7,0.25);
  beta[3] ~ normal(0.3,0.25);
  beta[4] ~ normal(0.3,0.25);
  beta[5] ~ normal(0,0.25);
  beta[6] ~ normal(0,0.25);

  logPmissing ~ normal(2.56,1.92); // based on logP fit

  eta ~ student_t(7,0,0.23);

  omega1[1] ~ normal(0,1.136);
  omega1[2] ~ normal(0,1.487);
  omega1[3] ~ normal(0,0.2);
  rho1   ~ lkj_corr(1);
  omega2[1] ~ normal(0,1.062);
  omega2[2] ~ normal(0,1.006);
  omega2[3] ~ normal(0,0.2);
  rho2	 ~ lkj_corr(1);
  sigma	 ~ normal(0,0.067);

  for(i in  1:nAnalytes){
  param1[i] ~ multi_student_t(7, miu1[i],Omega1);
  param2[i] ~ multi_student_t(7, miu2[i],Omega2);
  }
  
  if(run_estimation==1){
  logkObs ~ student_t(7,logkHat, sigma);	// observations
  }
}

generated quantities{

  real logkCond[nObs];
  real log_lik[nObs];
  real log_likPred[nObs];
  real logkHatPred[nObs];
  real logkaPred[nAnalytes];
  real logkwPred[nAnalytes];
  real logS2APred[nAnalytes];
  real<lower=0, upper=1> fracPred[nAnalytes];
  
  matrix[nAnalytes,nfiplot] logkHatPlotACond;
  matrix[nAnalytes,nfiplot] logkPlotACond;

  vector[3] paramPred1[nAnalytes];	 
  vector[3] paramPred2[nAnalytes];
  
  for(j in 1:nAnalytes){

  real etaPred[nAnalytes];

  etaPred[j] = student_t_rng(7,0, 0.23);
  fracPred[j]  = inv_logit(oddsi[j] + etaPred[j]);
  
  paramPred1[j] =  multi_student_t_rng(7,miu1[j],Omega1);
  paramPred2[j] =  multi_student_t_rng(7,miu2[j],Omega2);

  logkwPred[j]  = paramPred1[j, 1]*fracPred[j] + paramPred2[j, 1]*(1-fracPred[j]);
  logkaPred[j]  = paramPred1[j, 2]*fracPred[j] + paramPred2[j, 2]*(1-fracPred[j]);
  logS2APred[j] = paramPred1[j, 3]*fracPred[j] + paramPred2[j, 3]*(1-fracPred[j]);
  }
  
  for(i in 1:nObs){
    logkCond[i] = student_t_rng(7,logkHat[i], sigma);
    log_lik[i]	= student_t_lpdf(logkObs[i] | 7,logkHat[i], sigma);
    logkHatPred[i]	= hplcmodel(fi[i], logkwPred[analyte[i]], logkaPred[analyte[i]], logS2APred[analyte[i]]);
    log_likPred[i] = student_t_lpdf(logkObs[i] | 7,logkHatPred[i], sigma); 

  }
}
