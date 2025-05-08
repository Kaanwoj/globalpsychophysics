// global psychophysics model for matching p=1
// reference parameters are estimated as role-dependent
functions {
  real db_lambert(real std_b) {return 10 * log10(std_b * pi() / (10 ^ -6));}
  real db_spl(real std_l) {return 20 * log10(std_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real loud_to_bright(real std_loud, real alpha_l, real alpha_b, real beta_l,
                      real beta_b, real w_1, real rho_ltob, real rho_bfroml) {
    return db_lambert(
      inv(alpha_b) *
      pow(w_1 * alpha_l * pow(db_inv_spl(std_loud), beta_l) -
            w_1 * alpha_l * pow(db_inv_spl(rho_ltob), beta_l) +
            alpha_b * pow(db_inv_lambert(rho_bfroml), beta_b),
          inv(beta_b))
    );
  }
  real bright_to_loud(real std_bright, real alpha_l, real alpha_b, real beta_l,
                      real beta_b, real w_1, real rho_btol, real rho_lfromb) {
    return db_spl(
       inv(alpha_l) *
         pow(w_1 * alpha_b * pow(db_inv_lambert(std_bright), beta_b) -
               w_1 * alpha_b * pow(db_inv_lambert(rho_btol), beta_b) +
               alpha_l * pow(db_inv_lambert(rho_lfromb), beta_l),
             inv(beta_l))
    );
  }
}
data {
  int<lower=1> ntotal_lb;
  int<lower=1> ntotal_bl;
  int<lower=1> nstd_lb;
  int<lower=1> nstd_bl;
  vector<lower=1>[nstd_lb] std_lb;
  vector<lower=1>[nstd_bl] std_bl;
  int<lower=1> std_lb_idx[ntotal_lb];
  int<lower=1> std_bl_idx[ntotal_bl];
  vector<lower=0>[ntotal_lb] tgt_lb;
  vector<lower=0>[ntotal_bl] tgt_bl;
  vector<lower=0>[nstd_lb] sig_lb;
  vector<lower=0>[nstd_bl] sig_bl;
}
transformed data {
  int<lower=1, upper=1> alpha_b = 1;
}
parameters {
  real<lower=0> alpha_l;
  real<lower=0> beta_l;
  real<lower=0> beta_b;
  real<lower=0> w_1;
  real          rho_ltob;
  real          rho_bfroml;
  real          rho_btol;
  real          rho_lfromb;
}
transformed parameters {
  vector[nstd_lb] mu_lb;
  vector[nstd_bl] mu_bl;
  for (n in 1:nstd_lb)
    mu_lb[n] = loud_to_bright(std_lb[n], alpha_l, alpha_b, beta_l, beta_b, w_1,
                              rho_ltob, rho_bfroml);
  for (n in 1:nstd_bl)
    mu_bl[n] = bright_to_loud(std_bl[n], alpha_l, alpha_b, beta_l, beta_b, w_1,
                              rho_btol, rho_lfromb);
}
model {
  alpha_l ~ normal(12, 3); #(0.7, 0.1);
  beta_l ~ beta(3, 6);
  beta_b ~ beta(3, 6);
  w_1 ~ normal(1, .3);
  rho_ltob ~ normal(60, 20);
  rho_bfroml ~ normal(60, 20);
  rho_btol ~ normal(60, 20);
  rho_lfromb ~ normal(60, 20);
  target += normal_lpdf(tgt_lb | mu_lb[std_lb_idx], sig_lb[std_lb_idx]);
  target += normal_lpdf(tgt_bl | mu_bl[std_bl_idx], sig_bl[std_bl_idx]);
}
generated quantities {
  array[ntotal_lb] real tgt_lb_pred = normal_rng(mu_lb[std_lb_idx], sig_lb[std_lb_idx]);
  array[ntotal_bl] real tgt_bl_pred = normal_rng(mu_bl[std_bl_idx], sig_bl[std_bl_idx]);
}
