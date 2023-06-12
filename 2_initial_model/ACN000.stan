functions{
  real hplcmodel(real fi, real logkw, real logka, real logS2A){
    
    real logk;												// retention factor
    real S1;												// slope coefficient
    
    S1 = (logkw - logka)*(1+exp(logS2A));
    logk					= logkw - S1 * fi / (1 + exp(logS2A) * fi);
    
    return logk;
  }
}
data{
  int nAnalytes;											// number of analytes
  int nObs;												// number of observations
  int analyte[nObs];										// analytes indexes
  vector[nObs] logkObs;									// observed retention factors
  vector[nObs] fi;										// organic modifier content in the mobile phase
  real logP[nAnalytes];      
  int idxmissing[nAnalytes];
  int nfiplot;											// number of fi for plotting
  vector[nfiplot] fiplot;									// organic modifier content in the mobile phase
  
  int<lower = 0, upper = 1> run_estimation; // 0 for prior predictive, 1 for estimation
}

transformed data{
  vector[3] etaHat;
  etaHat					= rep_vector(0.0,3);
}

parameters{
  real logkaHat;										   // retention factor in acetonitrile for  analyte with MLOGP 0
  real logS2AHat;								           // mean curvature coefficient for acetonitrile in population
  real<lower = 0> sigmaadd;							   // standard deviation for residuals
  corr_matrix[3] rho;									   // correlation matrix
  vector<lower = 0>[3] omega;							   // diagonal elements of variance-covariance matrix
  real logkwHat;										   // mean value of logkw in population
  vector[3] eta[nAnalytes];							   // inter-individual variability
  real beta_logkw;									   // coefficient value for molecular descriptor
  real beta_logka;
  real nu;			// normality constant
  real logPmissing[nAnalytes];
}

transformed parameters{
  cov_matrix[3] Omega;									// variance-covariance matrix
  real logka[nAnalytes];									// slope coefficient for acetonitrile
  real logS2A[nAnalytes];									// curvature coefficient for acetonitrile
  real logkw[nAnalytes];									// retention factor in neat water
  vector[nObs] logkHat;									// mean value of logk in population
  
  Omega = quad_form_diag(rho, omega);	// diag_matrix(omega) * rho * diag_matrix(omega)
  
  for(j in 1:nAnalytes){
    logka[j]  = eta[j, 1] + logkaHat + beta_logka * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
    logS2A[j] = eta[j, 2] + logS2AHat;
    logkw[j]  = eta[j, 3] + logkwHat + beta_logkw * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
  }
  
  for(i in 1:nObs){
    logkHat[i] = hplcmodel(fi[i], logkw[analyte[i]], logka[analyte[i]], logS2A[analyte[i]]);
  }
}

model{
  logkwHat				~ normal(2, 5);
  logkaHat			    ~ normal(0, 5);
  logS2AHat				~ normal(log(2), 0.5);
    beta_logkw         	    ~ normal(1,0.5);
    beta_logka         	    ~ normal(1,0.5);
  logPmissing ~ normal(2.56,1.92);
  omega					~ normal(0,5);
  nu                      ~ gamma(2,0.1);
  rho						~ lkj_corr(1);
  sigmaadd				~ normal(0,1);
  eta						~ multi_student_t(nu,etaHat, Omega);	// inter-individual variability

  if(run_estimation==1){
  logkObs					~ normal(logkHat, sigmaadd);	// observations
}

}

generated quantities{
  
  vector[3] etaPred[nAnalytes];							// etas
  real logkaPred[nAnalytes];								// etention factor in acetonitrile
  real logS2APred[nAnalytes];								// curvature coefficient for acetonitrile
  real logkwPred[nAnalytes];								// retention factor in neat water
  vector[nObs] logkHatPred;						    	// predicted logk	
  real logkCond[nObs];
  real logkPred[nObs];
  real log_lik[nObs];
  real log_likPred[nObs];
  matrix[nAnalytes,nfiplot] logkHatPlotACond;
  matrix[nAnalytes,nfiplot] logkPlotACond;
  matrix[nAnalytes,nfiplot] logkHatPlotAPred;
  matrix[nAnalytes,nfiplot] logkPlotAPred;
  
  for(j in 1:nAnalytes){
    etaPred[j] = multi_student_t_rng(nu,etaHat, Omega);	// inter-individual variability
    logkaPred[j] = etaPred[j, 1] + logkaHat + beta_logka * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
    logS2APred[j]	= etaPred[j, 2] + logS2AHat;
    logkwPred[j]	= etaPred[j, 3] + logkwHat + beta_logkw * ((1-idxmissing[j])*logP[j]+idxmissing[j]*logPmissing[j]);
  }
  
  for(i in 1:nObs){
    logkHatPred[i]	= hplcmodel(fi[i], logkwPred[analyte[i]], logkaPred[analyte[i]], logS2APred[analyte[i]]);
    logkCond[i]		= normal_rng(logkHat[i], sigmaadd);
    logkPred[i]		= normal_rng(logkHatPred[i], sigmaadd);
    log_lik[i]		= normal_lpdf(logkObs[i] | logkHat[i], sigmaadd);
    log_likPred[i]	= normal_lpdf(logkObs[i] | logkHatPred[i], sigmaadd);
  }
  
  for(j in 1:nAnalytes){
    for(z in 1:nfiplot){
      logkHatPlotACond[j,z]	= hplcmodel(fiplot[z], logkw[j], logka[j], logS2A[j]);
      logkPlotACond[j,z]		= normal_rng(logkHatPlotACond[j,z], sigmaadd);
      
      logkHatPlotAPred[j,z]	= hplcmodel(fiplot[z], logkwPred[j], logkaPred[j], logS2APred[j]);
      logkPlotAPred[j,z]		= normal_rng(logkHatPlotAPred[j,z], sigmaadd);
      
    }
  }
  
}
