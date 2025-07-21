// global psychophysics model
// reference parameters are estimated as role-dependent
// One variance is estimated for each matching direction

functions {
  real db_lambert(real std_b) {return 10 * log10(std_b * pi() / (10 ^ -6));}
  real db_spl(real std_l) {return 20 * log10(std_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real weigh_fun(real p, real omega1, real omega) {
      return omega1 * pow(p, omega);
  }
  real loud_to_bright(real std_loud, real alpha_l, real alpha_b, real beta_l,
                      real beta_b, real rho_ltob, real rho_bfroml,
                      real omega1, real p, real omega) {
    return db_lambert(
      inv(alpha_b) *
      pow(weigh_fun(p, omega1, omega) * alpha_l * pow(db_inv_spl(std_loud), beta_l) -
            weigh_fun(p, omega1, omega) * alpha_l * pow(db_inv_spl(rho_ltob), beta_l) +
            alpha_b * pow(db_inv_lambert(rho_bfroml), beta_b),
          inv(beta_b))
    );
  }
  real bright_to_loud(real std_bright, real alpha_l, real alpha_b, real beta_l,
                      real beta_b, real rho_btol, real rho_lfromb,
                      real omega1, real p, real omega) {
    return db_spl(
       inv(alpha_l) *
         pow(weigh_fun(p, omega1, omega) * alpha_b * pow(db_inv_lambert(std_bright), beta_b) -
               weigh_fun(p, omega1, omega) * alpha_b * pow(db_inv_lambert(rho_btol), beta_b) +
               alpha_l * pow(db_inv_lambert(rho_lfromb), beta_l),
             inv(beta_l))
    );
  }
}
data {
  int<lower=1> ntotal;      // total ntrials
  int<lower=1> ntotal_lb;   // total ntrials per task
  int<lower=1> ntotal_bl;   // total ntrials per task
  int<lower=1> np;
  array[np] int p;
  int<lower=1> nstd_lb;
  int<lower=1> nstd_bl;
  vector<lower=1>[nstd_lb] std_lb;
  vector<lower=1>[nstd_bl] std_bl;
  array[ntotal_lb] int std_p_lb;
  array[ntotal_bl] int std_p_bl;
  int<lower=1> nstd_p_lb;
  int<lower=1> nstd_p_bl;
  vector<lower=0>[ntotal] tgt;
}
transformed data {
  int<lower=1, upper=1> alpha_b = 1;
}
parameters {
  real<lower=0> alpha_l;
  real<lower=0, upper=1> beta_l;
  real<lower=0, upper=1> beta_b;
  real<lower=0> omega1;
  real<lower=0> omega;
  real<lower=0> sig_lb;
  real<lower=0> sig_bl;
  real rho_ltob;
  real rho_bfroml;
  real rho_btol;
  real rho_lfromb;
}
transformed parameters {
  matrix[nstd_lb, np] mu_lb;
  matrix[nstd_bl, np] mu_bl;
  for (i in 1:np) {         // what if p and std is not fully crossed?
    for (j in 1:nstd_lb)    // column-major order
      mu_lb[j, i] = loud_to_bright(std_lb[j], alpha_l, alpha_b, beta_l, beta_b,
                                   rho_ltob, rho_bfroml, omega1, p[i], omega);
    for (j in 1:nstd_bl)
      mu_bl[j, i] = bright_to_loud(std_bl[j], alpha_l, alpha_b, beta_l, beta_b,
                                rho_btol, rho_lfromb, omega1, p[i], omega);
  };
}
model {
  alpha_l ~ normal(30, 20); #(0.7, 0.1);
  beta_l ~ beta(3, 6);
  beta_b ~ beta(3, 6);
  omega1 ~ normal(1, .3);       //lognormal(0, .1);
  omega ~ normal(0.6, .5);      //lognormal(0, .2);
  sig_lb ~ normal(3, 5);
  sig_bl ~ normal(3, 5);
  rho_ltob ~ normal(50, 20);
  rho_bfroml ~ normal(70, 20);
  rho_btol ~ normal(70, 20);
  rho_lfromb ~ normal(50, 20);

  tgt ~ normal(append_row(to_vector(mu_lb)[std_p_lb],
                          to_vector(mu_bl)[std_p_bl]),
               append_row(rep_vector(sig_lb, ntotal_lb), 
                          rep_vector(sig_bl, ntotal_bl)));
  
}
generated quantities {
  array[ntotal] real tgt_pred = normal_rng(
    append_row(to_vector(mu_lb)[std_p_lb], to_vector(mu_bl)[std_p_bl]),
    append_row(rep_vector(sig_lb, ntotal_lb), rep_vector(sig_bl, ntotal_bl)));
}
