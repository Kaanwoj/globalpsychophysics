// global psychophysics model for matching (p = 1)
// reference parameters are estimated as a constant sum
functions {
  real db_lambert(real std_b) {return 10 * log10(std_b * pi() / (10 ^ -6));}
  real db_spl(real std_l) {return 20 * log10(std_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real loud_to_bright(real std_loud, real alpha_l, real beta_l, real beta_b,
                        real w_1, real const_lb) {
    return db_lambert(
      pow(w_1 * alpha_l * pow(db_inv_spl(std_loud), beta_l) -
        const_lb, inv(beta_b))
    );
  }
  real bright_to_loud(real std_bright, real alpha_l, real beta_l,
                        real beta_b, real w_1, real const_bl) {
    return db_spl(
       inv(alpha_l) *
       pow(w_1 * pow(db_inv_lambert(std_bright), beta_b) -
         const_bl, inv(beta_l))
    );
  }
}
data {
  int<lower=1>    ntotal;
  int<lower=1>    ntotal_lb;
  int<lower=1>    ntotal_bl;
  int<lower=1>        nstd_lb;
  int<lower=1>        nstd_bl;
  vector<lower=1>[nstd_lb]      std_lb;
  vector<lower=1>[nstd_bl]      std_bl;
  int<lower=1>       std_lb_idx[ntotal_lb];
  int<lower=1>       std_bl_idx[ntotal_bl];
//vector<lower=0>[ntotal_lb]      tgt_lb;
//vector<lower=0>[ntotal_bl]      tgt_bl;
  vector<lower=0>[ntotal] tgt;
  vector<lower=0>[nstd_lb] sig_lb;
  vector<lower=0>[nstd_bl] sig_bl;
}
parameters {
  real<lower=0> alpha_l;
  real<lower=0, upper=1> beta_l;
  real<lower=0, upper=1> beta_b;
  real<lower=0> w_1;
  real          const_lb;
  real          const_bl;
}
transformed parameters {
  vector[nstd_lb] mu_lb;
  vector[nstd_bl] mu_bl;
  for (n in 1:nstd_lb)
    mu_lb[n] = loud_to_bright(std_lb[n], alpha_l, beta_l, beta_b, w_1, const_lb);
  for (n in 1:nstd_bl)
    mu_bl[n] = bright_to_loud(std_bl[n], alpha_l, beta_l, beta_b, w_1, const_bl);
}
model {
  alpha_l ~ normal(12, 3); #(0.7, 0.1);
  beta_l ~ beta(3, 6);
  beta_b ~ beta(3, 6);
  w_1 ~ normal(1, .3);
  const_lb ~ normal(0, 1);
  const_bl ~ normal(0, 1);
//target += normal_lpdf(tgt_lb | mu_lb[std_lb_idx], sig_lb[std_lb_idx]);
//target += normal_lpdf(tgt_bl | mu_bl[std_bl_idx], sig_bl[std_bl_idx]);
  tgt ~ normal(append_row(mu_lb[std_lb_idx], mu_bl[std_bl_idx])
}
generated quantities {
  array[ntotal_lb] real tgt_lb_pred = normal_rng(mu_lb[std_lb_idx], sig_lb[std_lb_idx]);
  array[ntotal_bl] real tgt_bl_pred = normal_rng(mu_bl[std_bl_idx], sig_bl[std_bl_idx]);
    // Log-likelihood for LOO
  vector[ntotal_lb + ntotal_bl] log_lik;
  // ... for loud to bright observations
  for (n in 1:ntotal_lb) {
    log_lik[n] = normal_lpdf(tgt_lb[n] | mu_lb[std_lb_idx[n]], sig_lb[std_lb_idx[n]]);
  }
  // ... for bright to loud observations
  for (n in 1:ntotal_bl) {
    log_lik[ntotal_lb + n] = normal_lpdf(tgt_bl[n] | mu_bl[std_bl_idx[n]], sig_bl[std_bl_idx[n]]);
  }
}
