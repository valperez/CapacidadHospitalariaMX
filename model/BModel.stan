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
  matrix[nestados, ndias] PHosp; //Proporción de hospitalizados para estado i día j
}

parameters {
  //Efectos autoregresivos
  real mu;
  real<lower=0> sigma;

  //Para la autoregresión con parámetros lambda
  real lambda[m];

  vector[nestados] alpha; //Agregamos un random effect
  real<lower=0> kappa[nestados];
  real<lower=0> sigma_kappa[nestados];

  //Efectos de estado
  real mu_estado[nestados];    
  real<lower=0> sigma_estado[nestados];
}

model {

  //Parámetros
  real alpha_model[nestados];
  real beta_model[nestados];
  vector[nestados] logit_p_estado;

  //Hiperparámetros
  mu           ~ normal(mu_mu_hiper, sigma_mu_hiper);
  sigma        ~ cauchy(0, sigma_mu_hiper);
  mu_estado    ~ normal(mu_mu_hiper, sigma_mu_hiper);
  sigma_estado ~ cauchy(0, sigma_mu_hiper);
  sigma_kappa  ~ normal(0, sigma_kappa_hiper);

  //Creamos los parámetros
  lambda ~ normal(mu, sigma);
  alpha  ~ normal(mu_estado, sigma_estado);
  kappa  ~ normal(0, sigma_kappa);

  //Loopeamos
  for (t in (m + 1):ndias){
    logit_p_estado = alpha;
    for (k in 1:m){
      logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*PHosp[1:nestados,t-k];
    }
    PHosp[1:nestados,t] ~ beta_proportion(inv_logit(logit_p_estado[1:nestados]), kappa[1:nestados]);
  }
}

//Adapted from https://jwalton.info/Stan-posterior-predictives/
generated quantities {
  // Generate posterior predictives
  matrix[nestados, dias_predict + ndias] HospPred;
  vector[nestados] logit_p_estado;

  // First m points are to start model
  HospPred[1:nestados, 1:m] = PHosp[1:nestados, 1:m];

  // Posterior dist for observed
  for (t in (m + 1):ndias){
    logit_p_estado = alpha;
    for (k in 1:m){
      logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*PHosp[1:nestados,t-k];
    }
    HospPred[1:nestados, t] = to_vector( beta_proportion_rng(inv_logit(logit_p_estado), kappa) );
  }

  // Posterior dist for unobserved but still using some observed
  for (t in (ndias + 1):(ndias + m)){
    logit_p_estado = alpha;
    for (k in 1:m){
      if (t - k <= ndias){
        logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*PHosp[1:nestados,t-k];
      } else {
        logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*HospPred[1:nestados,t-k];
      }
    }
    HospPred[1:nestados, t] = to_vector( beta_proportion_rng(inv_logit(logit_p_estado), kappa) );
  }

  // Posterior dist for unobserved 
  for (t in (ndias + m + 1):(ndias +  dias_predict)){
    logit_p_estado = alpha;
    for (k in 1:m){
      logit_p_estado[1:nestados] = logit_p_estado[1:nestados] + lambda[k]*HospPred[1:nestados,t-k];
    }
    HospPred[1:nestados, t] = to_vector( beta_proportion_rng(inv_logit(logit_p_estado), kappa) );
  }
}