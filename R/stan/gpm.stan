// global psychophysics model
// reference parameters are estimated as potentially role-dependent
// The variances are estimated depending on sigindex input

functions {
  real db_lambert(real std_b) {return 10 * log10(std_b * pi() / (10 ^ -6));}
  real db_spl(real std_l) {return 20 * log10(std_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real weigh_fun(real p, real omega_1, real omega) {
      return omega_1 * pow(p, omega);
  }
  real loud_to_bright(real std_loud, real alpha_l, real alpha_b, real beta_l,
                      real beta_b, real rho_ltob, real rho_bfroml,
                      real omega_1, real p, real omega) {
    return db_lambert(
      pow(
        inv(alpha_b) *
          (weigh_fun(p, omega_1, omega) * alpha_l * 
             pow(db_inv_spl(std_loud), beta_l) -
           weigh_fun(p, omega_1, omega) * alpha_l * 
             pow(db_inv_spl(rho_ltob), beta_l) +
           alpha_b * pow(db_inv_lambert(rho_bfroml), beta_b)
          ),
        inv(beta_b))
    );
  }
  real bright_to_loud(real std_bright, real alpha_l, real alpha_b, real beta_l,
                      real beta_b, real rho_btol, real rho_lfromb,
                      real omega_1, real p, real omega) {
    return db_spl(
      pow(
        inv(alpha_l) *
          (weigh_fun(p, omega_1, omega) * alpha_b *
             pow(db_inv_lambert(std_bright), beta_b) -
           weigh_fun(p, omega_1, omega) * alpha_b *
             pow(db_inv_lambert(rho_btol), beta_b) +
           alpha_l * pow(db_inv_spl(rho_lfromb), beta_l)
          ),
        inv(beta_l))
    );
  }
}
data {
  int<lower=1> ntotal;      // total ntrials
  int<lower=1> ncond;
  int<lower=1> nsig;
  array[ncond] int<lower=1, upper=2> task;   // 1: loud_bright 2: bright_loud
  array[ncond] int<lower=1, upper=3> p;
  vector<lower=1>[ncond] std;
  array[ntotal] int<lower=1> idx;
  array[ntotal] int<lower=1> sigidx;
 // array[ntotal] real<lower=0> tgt;
  array[ntotal] real tgt;
  int<lower=0, upper=1> onlyprior;
  real lowestsoundstd;
  real lowestlightstd;
  
  // priors
  real<lower=0> alpha_l_logmu;
  real<lower=0> alpha_l_logsigma;
  real          beta_b_logmu;
  real<lower=0> beta_b_logsigma;
  real          beta_l_logmu;
  real<lower=0> beta_l_logsigma;
  real<lower=0> rho_btol_mu;
  real<lower=0> rho_btol_sigma;
  real<lower=0> rho_bfroml_mu;
  real<lower=0> rho_bfroml_sigma;
  real<lower=0> rho_ltob_mu;
  real<lower=0> rho_ltob_sigma;
  real<lower=0> rho_lfromb_mu;
  real<lower=0> rho_lfromb_sigma;
  real omega_1_logmu;
  real<lower=0> omega_1_logsigma;
  real          omega_logmu;
  real<lower=0> omega_logsigma;
  real<lower=0> sig_mu;
  real<lower=0> sig_sigma;
}
transformed data {
  int<lower=1, upper=1> alpha_b = 1;
}
parameters {
  real<lower=0> alpha_l;
  real<lower=0> beta_l;
  real<lower=0> beta_b;
  real<lower=0> omega_1;
  real<lower=0> omega;
  vector<lower=0>[nsig] sig;
  real<upper=lowestsoundstd> rho_ltob;
  real rho_bfroml;
  real<upper=lowestlightstd> rho_btol;
  real rho_lfromb;
}
transformed parameters {
  array[ncond] real mu; 
  for (i in 1:ncond) {
    if (task[i] == 1) {
      mu[i] = loud_to_bright(std[i], alpha_l, alpha_b, beta_l, beta_b,
                             rho_ltob, rho_bfroml, omega_1, p[i], omega);
    } else if (task[i] == 2) {
      mu[i] = bright_to_loud(std[i], alpha_l, alpha_b, beta_l, beta_b,
                             rho_btol, rho_lfromb, omega_1, p[i], omega);

    }
  }
}
model {
  target += lognormal_lpdf(alpha_l | alpha_l_logmu, alpha_l_logsigma);
  target += lognormal_lpdf(beta_b | beta_b_logmu, beta_b_logsigma);
  target += lognormal_lpdf(beta_l | beta_l_logmu, beta_l_logsigma);
  target += normal_lpdf(rho_btol | rho_btol_mu, rho_btol_sigma) -
              normal_lcdf(lowestlightstd | rho_btol_mu, rho_btol_sigma);
              // normal_lcdf with a constant does not change the sampling, but
              // defining the variable with upper=lowest... does the trucation
  target += normal_lpdf(rho_bfroml | rho_bfroml_mu, rho_bfroml_sigma);
  target += normal_lpdf(rho_ltob | rho_ltob_mu, rho_ltob_sigma) -
              normal_lcdf(lowestsoundstd | rho_ltob_mu, rho_ltob_sigma);
  target += normal_lpdf(rho_lfromb | rho_lfromb_mu, rho_lfromb_sigma);
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
  array[ntotal] real log_lik;
  for (i in 1:ntotal) {
    log_lik[i] = normal_lpdf(tgt[i] | mu[idx[i]], sig[sigidx[i]]);
  }
}
