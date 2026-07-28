data {
    int <lower=0> N;       //the number of individuals
    real <lower=0> age[N]; // Age 
    real Y[N];             // log  titer (the new seropositivity)
}





parameters {
    real<lower =0> lambda;     // Seroconversion rate
    real<lower =0> rho;        // Seroreversion rate


    real<lower=0> mu_N;        // mean of the distribution for seronegatives
    real<lower=0> mu_P;        // difference between mu_N and the mean of the seropositives
    real<lower=0, upper=10> sigma_N;  // std of the different distributions
    real<lower=0, upper=10> sigma_P;

}


transformed parameters{
  real muPos;
  muPos = mu_N + mu_P;  // Mean of the positives
  
  // Serocatalytic model
  vector[N] PM; // Proportion of seropositives for each age class

  for(i in 1:N){PM[i]= lambda/(lambda+rho)*(1-exp(-(lambda+rho)*age[i]));}
  
}








model {
    real LPS[2];

    real I1; // infected by both MAYV and CHIKV   
    real I2; 
    
  
  
    // Priors
    lambda ~ uniform(0,100); 
    rho ~ uniform(0,10);     

    mu_N ~ uniform(0,13);
    mu_P ~ normal(log(100), 0.1); 
    sigma_N ~ uniform(0,5); 
    sigma_P ~ uniform(0,5); 


 
    // Log Likelihood computation
    for(j in 1:N){
        I2  = PM[j];
        I1  = 1 - PM[j];


        if(I1  == 0){LPS[1] = -10000;
        }else{
            LPS[1] = log(I1) + normal_lpdf(Y[j] | mu_N, sigma_N);} 
      
        
        if(I2 == 0){LPS[2] = -10000;
        }else{
            LPS[2] = log(I2) + normal_lpdf(Y[j] | muPos, sigma_P);}
            

            target += log(sum(exp(LPS)));
    }
}



