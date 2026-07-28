data {
    int <lower=0> N;       //the number of individuals
    real <lower=0> age[N]; // Age 
    real Y[N];             // log  titer (the new seropositivity)
}





parameters {
    real<lower =0> lambda;       // seroconversion rate for primo-infection
    real<lower =0> lambda_bis;   // seroconversion rate from the recovered state
    real<lower =0> rho;          // seroreversion rate


    real<lower=0> mu_N;          // mean of the distribution of the negatives
    real<lower=0> mu_R;          // difference between mu_N and the mean of the recovered
    real<lower=0> mu_P;          // difference between the mean of the recovered and the mean of the positives
    real<lower=0, upper=10> sigma_N; // std of the different distributions
    real<lower=0, upper=10> sigma_R;
    real<lower=0, upper=10> sigma_P;

}


transformed parameters{
  vector[N] Neg;    // Proportion of sero-negative for each age class
  vector[N] R;      // Proportion of recovered for each age class
  real nu = rho + lambda_bis; // So the equation is easily readable
  
  for(i in 1:N){R[i]= (rho/(nu))*(1 - exp(-(nu)*age[i])) - (rho/(nu - lambda))*(exp(-lambda*age[i]) - exp(-(nu)*age[i]));}
  for(i in 1:N){Neg[i]= exp(-lambda * age[i]);}

   
   
  real muPos;
  real muR;
  muR = mu_N + mu_R;
  muPos = mu_N + mu_R + mu_P;  // Mean of the positives
}








model {
    real LPS[3];

    real I1; 
    real I2; 
    real I3; 
    
    // Priors
    lambda ~ uniform(0,100);
    lambda_bis ~ uniform(0,100);
    rho ~ uniform(0,10);

    mu_N ~ uniform(0,13);
    mu_R ~ normal(log(100)/2,0.1);
    mu_P ~ normal(log(100)/2,0.1);
    sigma_N ~ uniform(0,5); 
    sigma_R ~ uniform(0,5); 
    sigma_P ~ uniform(0,5); 


    // Log Likelihood computation
    for(j in 1:N){
        I3  = R[j];
        I1  = Neg[j];
        I2  = 1 - R[j] - Neg[j];


        if(I3  == 0){LPS[3] = -10000;
        }else{
            LPS[3] = log(I3) + normal_lpdf(Y[j] | muR, sigma_R);} 
      
        
        if(I2 == 0){LPS[2] = -10000;
        }else{
            LPS[2] = log(I2) + normal_lpdf(Y[j] | muPos, sigma_P);}
            
        
        if(I1 == 0){LPS[1] = -10000;
        }else{
            LPS[1] = log(I1) + normal_lpdf(Y[j] | mu_N, sigma_N);}

            target += log(sum(exp(LPS)));
    }
}



