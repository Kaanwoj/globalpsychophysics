// global psychophysics model
// reference parameters are estimated as potentially role-dependent
// The variances are estimated depending on sigindex input

functions {
  real db_lambert(real std_b) {return 10 * log10(std_b * pi() / (10 ^ -6));}
  real db_spl(real std_l) {return 20 * log10(std_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real db_disp(real x_s) { return 40 * log10(x_s / 0.001578822); }
  real db_inv_disp(real db_s) { return 0.001578822 * 10^(db_s / 40); }
  real weigh_fun(real p, real omega_1, real omega) {
      return omega_1 * pow(p, omega);
  }
  real db_convert(real x, int modality) {
    if (modality == 1) return db_spl(x);
    else if (modality == 2) return db_lambert(x);
    else if (modality == 3) return db_disp(x);
    else reject("Unknown modality: ", modality, ". Please choose either 1 (loudnsess), 2 (brightness), 3 (vibration strength");
  }
    real db_inverse(real x, int modality) {
    if (modality == 1) return db_inv_spl(x);
    else if (modality == 2) return db_inv_lambert(x);
    else if (modality == 3) return db_inv_disp(x);
    else reject("Unknown modality: ", modality, ". Please choose either 1 (loudnsess), 2 (brightness), 3 (vibration strength");
  }
  
  real gpm_predict(real std, real alpha_std, real alpha_tgt,
                      real beta_std, real beta_tgt, 
                      real rho_std, real rho_tgt,
                      real omega_1, int p, real omega,
                      int std_modality, int tgt_modality) {
    return db_convert(
      pow(
        inv(alpha_tgt) *
          (weigh_fun(p, omega_1, omega) * alpha_std * 
             pow(db_inverse(std, std_modality), beta_std) -
           weigh_fun(p, omega_1, omega) * alpha_std * 
             pow(db_inverse(rho_std, std_modality), beta_std) +
           alpha_tgt * pow(db_inverse(rho_tgt, tgt_modality), beta_tgt)
          ),
        inv(beta_tgt))
    , tgt_modality);
  }
}

data {
  int<lower=1> ntotal;      // total ntrials
  int<lower=1> ncond;
  int<lower=1> nsig;
  array[ncond] int<lower=1, upper=3> std_modality; // 1 = loudness, 2 = brightness, 3 = vibration strength
  array[ncond] int<lower=1, upper=3> tgt_modality; // 1 = loudness, 2 = brightness, 3 = vibration strength
  array[ncond] int<lower=1, upper=3> p;
  vector<lower=1>[ncond] std;
  array[ntotal] int<lower=1> idx;
  array[ntotal] int<lower=1> sigidx;
  array[ntotal] real<lower=0> tgt;
  int<lower=0, upper=1> onlyprior;
  
  // priors
  real<lower=0> alpha_l_logmu;
  real<lower=0> alpha_l_logsigma;
  real<lower=0> alpha_s_logmu;
  real<lower=0> alpha_s_logsigma;
  real          beta_l_logmu;
  real<lower=0> beta_l_logsigma;
  real          beta_b_logmu;
  real<lower=0> beta_b_logsigma;
  real          beta_s_logmu;
  real<lower=0> beta_s_logsigma;
  real<lower=0> rho_ltol_mu;
  real<lower=0> rho_ltol_sigma;
  real<lower=0> rho_ltob_mu;
  real<lower=0> rho_ltob_sigma;
  real<lower=0> rho_ltos_mu;
  real<lower=0> rho_ltos_sigma;
  real<lower=0> rho_lfroml_mu;
  real<lower=0> rho_lfroml_sigma;
  real<lower=0> rho_lfromb_mu;
  real<lower=0> rho_lfromb_sigma;
  real<lower=0> rho_lfroms_mu;
  real<lower=0> rho_lfroms_sigma;
  real<lower=0> rho_btob_mu;
  real<lower=0> rho_btob_sigma;
  real<lower=0> rho_btol_mu;
  real<lower=0> rho_btol_sigma;
  real<lower=0> rho_btos_mu;
  real<lower=0> rho_btos_sigma;
  real<lower=0> rho_bfromb_mu;
  real<lower=0> rho_bfromb_sigma;
  real<lower=0> rho_bfroml_mu;
  real<lower=0> rho_bfroml_sigma;
  real<lower=0> rho_bfroms_mu;
  real<lower=0> rho_bfroms_sigma;
  real<lower=0> rho_stos_mu;
  real<lower=0> rho_stos_sigma;
  real<lower=0> rho_stol_mu;
  real<lower=0> rho_stol_sigma;
  real<lower=0> rho_stob_mu;
  real<lower=0> rho_stob_sigma;
  real<lower=0> rho_sfroms_mu;
  real<lower=0> rho_sfroms_sigma;
  real<lower=0> rho_sfroml_mu;
  real<lower=0> rho_sfroml_sigma;
  real<lower=0> rho_sfromb_mu;
  real<lower=0> rho_sfromb_sigma;
  real          omega_1_logmu;
  real<lower=0> omega_1_logsigma;
  real          omega_logmu;
  real<lower=0> omega_logsigma;
  real<lower=0> sig_mu;
  real<lower=0> sig_sigma;
}

transformed data {
  real alpha_b = 1.0;
}

parameters {
  real<lower=0> alpha_l;
  real<lower=0> alpha_s;
  real<lower=0> beta_l;
  real<lower=0> beta_b;
  real<lower=0> beta_s;
  real<lower=0> omega_1;
  real<lower=0> omega;
  vector<lower=0>[nsig] sig;
  real rho_ltob;
  real rho_ltos;
  real rho_ltol;
  real rho_lfromb;
  real rho_lfroms;
  real rho_lfroml;
  real rho_btol;
  real rho_btos;
  real rho_btob;
  real rho_bfroml;
  real rho_bfroms;
  real rho_bfromb;
  real rho_stol;
  real rho_stob;
  real rho_stos;
  real rho_sfroml;
  real rho_sfromb;
  real rho_sfroms;
}
                      
transformed parameters {
  array[ncond] real mu; 
  for (i in 1:ncond) {
    real alpha_std_i;
    real alpha_tgt_i;
    real beta_std_i;
    real beta_tgt_i;
    real rho_std_i;
    real rho_tgt_i;
    // set standard parameters for gpm_predict 
    if (std_modality[i] == 1){
      alpha_std_i = alpha_l;
      beta_std_i  = beta_l;
      if (tgt_modality[i] == 1){
        rho_std_i   = rho_ltol;
      }else if(tgt_modality[i] == 2){
        rho_std_i   = rho_ltob;
      }else if(tgt_modality[i] == 3){
        rho_std_i   = rho_ltos;
      }
    } else if (std_modality[i] == 2){
      alpha_std_i = alpha_b;
      beta_std_i  = beta_b;
      if (tgt_modality[i] == 1){
        rho_std_i   = rho_btol;
      }else if(tgt_modality[i] == 2){
        rho_std_i   = rho_btob;
      }else if(tgt_modality[i] == 3){
        rho_std_i   = rho_btos;
      }
    } else if (std_modality[i] == 3){
      alpha_std_i = alpha_s;
      beta_std_i  = beta_s;
      if (tgt_modality[i] == 1){
        rho_std_i   = rho_stol;
      }else if(tgt_modality[i] == 2){
        rho_std_i   = rho_stob;
      }else if(tgt_modality[i] == 3){
        rho_std_i   = rho_stos;
      }
    }
    // set target parameters for gpm_predict 
    if (tgt_modality[i] == 1){
      alpha_tgt_i = alpha_l;
      beta_tgt_i  = beta_l;
      if (std_modality[i] == 1){
        rho_tgt_i   = rho_lfroml;
      }else if(std_modality[i] == 2){
        rho_tgt_i   = rho_lfromb;
      }else if(std_modality[i] == 3){
        rho_tgt_i   = rho_lfroms;
      }
    } else if (tgt_modality[i] == 2){
      alpha_tgt_i = alpha_b;
      beta_tgt_i  = beta_b;
      if (std_modality[i] == 1){
        rho_tgt_i   = rho_bfroml;
      }else if(std_modality[i] == 2){
        rho_tgt_i   = rho_bfromb;
      }else if(std_modality[i] == 3){
        rho_tgt_i   = rho_bfroms;
      }
    } else if (tgt_modality[i] == 3){
      alpha_tgt_i = alpha_s;
      beta_tgt_i  = beta_s;
      if (std_modality[i] == 1){
        rho_tgt_i   = rho_sfroml;
      }else if(std_modality[i] == 2){
        rho_tgt_i   = rho_sfromb;
      }else if(std_modality[i] == 3){
        rho_tgt_i   = rho_sfroms;
      }
    }
    // Fix alpha_s in loud_strong or stron_loud, since one alpha needs to be 
    // fixed and alpha_b is not here
    //if (std_modality[i] != 2 && tgt_modality[i] != 2) {
      // loud<->strong: fix alpha_s side
    //  if (tgt_modality[i] == 3) alpha_tgt_i = 1.0;
    //  else if (std_modality[i] == 3) alpha_std_i = 1.0;
    //}
    // call gpm_predict
    mu[i] = gpm_predict(std[i],
                        alpha_std_i, alpha_tgt_i,
                        beta_std_i,  beta_tgt_i,
                        rho_std_i,   rho_tgt_i,
                        omega_1, p[i], omega,
                        std_modality[i], tgt_modality[i]);
  }
}

model {
  target += lognormal_lpdf(alpha_l | alpha_l_logmu, alpha_l_logsigma);
  target += lognormal_lpdf(alpha_s | alpha_s_logmu, alpha_s_logsigma);
  target += lognormal_lpdf(beta_b | beta_b_logmu, beta_b_logsigma);
  target += lognormal_lpdf(beta_l | beta_l_logmu, beta_l_logsigma);
  target += lognormal_lpdf(beta_s | beta_s_logmu, beta_s_logsigma);
  target += normal_lpdf(rho_ltol | rho_ltol_mu, rho_ltol_sigma);
  target += normal_lpdf(rho_ltob | rho_ltob_mu, rho_ltob_sigma);
  target += normal_lpdf(rho_ltos | rho_ltos_mu, rho_ltos_sigma);
  target += normal_lpdf(rho_lfroml | rho_lfroml_mu, rho_lfroml_sigma);
  target += normal_lpdf(rho_lfromb | rho_lfromb_mu, rho_lfromb_sigma);
  target += normal_lpdf(rho_lfroms | rho_lfroms_mu, rho_lfroms_sigma);
  target += normal_lpdf(rho_btob | rho_btob_mu, rho_btob_sigma);
  target += normal_lpdf(rho_btol | rho_btol_mu, rho_btol_sigma);
  target += normal_lpdf(rho_btos | rho_btos_mu, rho_btos_sigma);
  target += normal_lpdf(rho_bfromb | rho_bfromb_mu, rho_bfromb_sigma);
  target += normal_lpdf(rho_bfroml | rho_bfroml_mu, rho_bfroml_sigma);
  target += normal_lpdf(rho_bfroms | rho_bfroms_mu, rho_bfroms_sigma);
  target += normal_lpdf(rho_stos | rho_stos_mu, rho_stos_sigma);
  target += normal_lpdf(rho_stob | rho_stob_mu, rho_stob_sigma);
  target += normal_lpdf(rho_stol | rho_stol_mu, rho_stol_sigma);
  target += normal_lpdf(rho_sfroms | rho_sfroms_mu, rho_sfroms_sigma);
  target += normal_lpdf(rho_sfromb | rho_sfromb_mu, rho_sfromb_sigma);
  target += normal_lpdf(rho_sfroml | rho_sfroml_mu, rho_sfroml_sigma);
  target += lognormal_lpdf(omega_1 | omega_1_logmu, omega_1_logsigma);
  target += lognormal_lpdf(omega | omega_logmu, omega_logsigma);
  target += normal_lpdf(sig | sig_mu, sig_sigma) -
              normal_lccdf(0 | sig_mu, sig_sigma);
  
  if (!onlyprior) {
    target += normal_lpdf(tgt | to_vector(mu)[idx], to_vector(sig[sigidx]));
  }
}

generated quantities {
  array[ntotal] real tgt_pred = normal_rng(
    to_vector(mu)[idx],
    to_vector(sig[sigidx]));
}
