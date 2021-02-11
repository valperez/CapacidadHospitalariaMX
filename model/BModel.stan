data {
  int<lower=1> nestados;                         //Cantidad de estados incluidos en el modelo
  int<lower=1> ndias;                            //Cantidad de días modelados
  int<lower=0> dias_predict;                     //Cantidad de días a predecir (futuro)
  int<lower=0> m;                                //Autoregresive order m
  real<lower=0> sigma_mu_hiper;                     //Sigma del hiperparámetro de la media
  real<lower=0> mu_mu_hiper;                        //Media del hiperparámetro de la media
  real<lower=0> sigma_sigma_hiper;
  real<lower=0> sigma_kappa_hiper;
  real<lower=0> sigma_estado_hiper;
  real<lower=0> mu_mu_time_hiper;
  real<lower=0> sigma_sigma_time_hiper;
  matrix[nestados, ndias] PHosp; //Proporción de hospitalizados para estado i día j
}

parameters {
  //Efectos autoregresivos
  real mu;
  real<lower=0> sigma;
  
  //Para la autoregresión con parámetros lambda
  real lambda[m];
  
  real mu_kappa;
  real kappa_sqrt[nestados];
  
  
  vector[nestados] alpha; //Agregamos un random effect
  real<lower=0> kappa[nestados];
  real<lower=0> sigma_kappa[nestados];
  
  //Efectos de estado
  real mu_estado[nestados];    
  real<lower=0> sigma_estado[nestados];
  
  //Seasonal effect for both 
  real mu_time[nestados];
  real<lower=0> sigma_time[nestados];
  vector[nestados] beta_yearly_cosine;
  vector[nestados] beta_yearly_sine;
  vector[nestados] beta_weekly_sine;
  vector[nestados] beta_weekly_cosine;
  vector[nestados] beta_monthly_sine;
  vector[nestados] beta_monthly_cosine;
  
  
  
}


model {
  
  //Parámetros
  vector[nestados] phi;
  vector[nestados] logit_p_estado;
  
  //Hiperparámetros
  mu           ~ normal(mu_mu_hiper, sigma_mu_hiper);
  sigma        ~ cauchy(0, sigma_mu_hiper);
  
  mu_time      ~ normal(mu_mu_time_hiper, sigma_sigma_time_hiper);
  sigma_time   ~ cauchy(0, sigma_sigma_time_hiper);
  
  mu_estado    ~ normal(mu_mu_hiper, sigma_mu_hiper);
  sigma_estado ~ cauchy(0, sigma_mu_hiper);
  sigma_kappa  ~ normal(0, sigma_kappa_hiper);
  mu_kappa     ~ normal(mu_mu_hiper, sigma_mu_hiper);
  
  //Creamos los parámetros
  lambda               ~ normal(mu, sigma);
  beta_yearly_cosine   ~ normal(mu_time, sigma_time);
  beta_yearly_sine     ~ normal(mu_time, sigma_time);
  beta_weekly_cosine   ~ normal(mu_time, sigma_time);
  beta_weekly_sine     ~ normal(mu_time, sigma_time);
  beta_monthly_sine    ~ normal(mu_time, sigma_time);
  beta_monthly_cosine  ~ normal(mu_time, sigma_time);
  alpha                ~ normal(mu_estado, sigma_estado);
  kappa                ~ normal(mu_kappa, sigma_kappa);
  kappa_sqrt           ~ normal(mu_kappa, sigma_kappa);
  
  
  //Loopeamos
  for (t in (m + 1):ndias){
    logit_p_estado = alpha + cos(2*pi()*t/365) * beta_yearly_cosine + 
                        sin(2*pi()*t/365) * beta_yearly_sine + 
                        cos(2*pi()*t/30)  * beta_monthly_cosine + 
                        sin(2*pi()*t/30)  * beta_monthly_sine +
                        cos(2*pi()*t/7) * beta_weekly_cosine + 
                        sin(2*pi()*t/7) * beta_weekly_sine;
    for (k in 1:m){
      logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*PHosp[1:nestados,t-k];
    }
    phi = to_vector( kappa ) +  to_vector( kappa_sqrt ) .* to_vector( square(PHosp[1:nestados,t-1] - PHosp[1:nestados,t-m]) );
    PHosp[1:nestados,t] ~ beta_proportion(inv_logit(logit_p_estado), phi);
  }
}

//Adapted from https://jwalton.info/Stan-posterior-predictives/
generated quantities {
  // Generate posterior predictives
  matrix[nestados, dias_predict + ndias] HospPred;
  vector[nestados] logit_p_estado;
  vector[nestados] phi;
  
  // First m points are to start model
  HospPred[1:nestados, 1:m] = PHosp[1:nestados, 1:m];

  // Posterior dist for observed
  for (t in (m + 1):ndias){
    logit_p_estado = alpha + cos(2*pi()*t/365) * beta_yearly_cosine + 
                    sin(2*pi()*t/365) * beta_yearly_sine + 
                    cos(2*pi()*t/30)  * beta_monthly_cosine + 
                    sin(2*pi()*t/30)  * beta_monthly_sine +
                    cos(2*pi()*t/7) * beta_weekly_cosine + 
                    sin(2*pi()*t/7) * beta_weekly_sine;
    for (k in 1:m){
      logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*PHosp[1:nestados,t-k];
    }
    phi = to_vector( kappa ) +  to_vector( kappa_sqrt ) .* to_vector( square(PHosp[1:nestados,t-1] - PHosp[1:nestados,t-m]) );
    HospPred[1:nestados, t] = to_vector( beta_proportion_rng(inv_logit(logit_p_estado), phi) );
  }
  
  // Posterior dist for unobserved but still using some observed
  for (t in (ndias + 1):(ndias + m)){
    logit_p_estado = alpha + cos(2*pi()*t/365) * beta_yearly_cosine + 
                    sin(2*pi()*t/365) * beta_yearly_sine + 
                    cos(2*pi()*t/30)  * beta_monthly_cosine + 
                    sin(2*pi()*t/30)  * beta_monthly_sine +
                    cos(2*pi()*t/7) * beta_weekly_cosine + 
                    sin(2*pi()*t/7) * beta_weekly_sine;
    for (k in 1:m){
      if (t - k <= ndias){
        logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*PHosp[1:nestados,t-k];
      } else {
        logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*HospPred[1:nestados,t-k];
      }
    }
    phi = to_vector( kappa ) +  to_vector( kappa_sqrt ) .* to_vector( square(HospPred[1:nestados,t-1] - PHosp[1:nestados,t-m]) );
    HospPred[1:nestados, t] = to_vector( beta_proportion_rng(inv_logit(logit_p_estado), phi) );
  }
  
  // Posterior dist for unobserved 
  for (t in (ndias + m + 1):(ndias +  dias_predict)){
    logit_p_estado = alpha + cos(2*pi()*t/365) * beta_yearly_cosine + 
                    sin(2*pi()*t/365) * beta_yearly_sine + 
                    cos(2*pi()*t/30)  * beta_monthly_cosine + 
                    sin(2*pi()*t/30)  * beta_monthly_sine +
                    cos(2*pi()*t/7) * beta_weekly_cosine + 
                    sin(2*pi()*t/7) * beta_weekly_sine;
    for (k in 1:m){
      logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*HospPred[1:nestados,t-k];
    }
    phi = to_vector( kappa ) +  to_vector( kappa_sqrt ) .* to_vector( square(HospPred[1:nestados,t-1] - HospPred[1:nestados,t-m]) );
    HospPred[1:nestados, t] = to_vector( beta_proportion_rng(inv_logit(logit_p_estado), phi) );
  }
}

