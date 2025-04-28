// global psychophysics model
// reference parameters are estimated as a constant sum
functions {
  real db_lambert(real x_b) {return 10 * log10(x_b * pi() / (10 ^ -6));}
  real db_spl(real x_l) {return 20 * log10(x_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real loud_to_bright(real x_loud, real alpha_l, real beta_l, real beta_b,
                        real w_1, real const_lb) {
    return db_lambert(
      pow(w_1 * alpha_l * pow(db_inv_spl(x_loud), beta_l) -
        const_lb, inv(beta_b))
    );
  }
  real bright_to_loud(real x_bright, real alpha_l, real beta_l,
                        real beta_b, real w_1, real const_bl) {
    return db_spl(
       inv(alpha_l) *
       pow(w_1 * pow(db_inv_lambert(x_bright), beta_b) -
         const_bl, inv(beta_l))
    );
  }
}
data {
  int<lower=1>    ntotal_lb;
  int<lower=1>    ntotal_bl;
  int<lower=1>      ntrials;
  int<lower=1>        nx_lb;
  int<lower=1>        nx_bl;
  vector<lower=1>[nx_lb]      x_lb;
  vector<lower=1>[nx_bl]      x_bl;
  int<lower=1>       x_lb_idx[ntotal_lb];
  int<lower=1>       x_bl_idx[ntotal_bl];
  vector<lower=0>[ntotal_lb]      y_lb;
  vector<lower=0>[ntotal_bl]      y_bl;
  vector<lower=0>[nx_lb] sig_lb;
  vector<lower=0>[nx_bl] sig_bl;
}
parameters {
  real<lower=0> alpha_l;
  real<lower=0> beta_l;
  real<lower=0> beta_b;
  real<lower=0> w_1;
  real          const_lb;
  real          const_bl;
}
transformed parameters {
  vector[nx_lb] mu_lb;
  vector[nx_bl] mu_bl;
  for (n in 1:nx_lb)
    mu_lb[n] = loud_to_bright(x_lb[n], alpha_l, beta_l, beta_b, w_1, const_lb);
  for (n in 1:nx_bl)
    mu_bl[n] = bright_to_loud(x_bl[n], alpha_l, beta_l, beta_b, w_1, const_bl);
}
model {
  alpha_l ~ normal(12, 3); #(0.7, 0.1);
  beta_l ~ beta(3, 6);
  beta_b ~ beta(3, 6);
  w_1 ~ normal(1, .3);
  const_lb ~ normal(0, 1);
  const_bl ~ normal(0, 1);
  target += normal_lpdf(y_lb | mu_lb[x_lb_idx], sig_lb[x_lb_idx]);
  target += normal_lpdf(y_bl | mu_bl[x_bl_idx], sig_bl[x_bl_idx]);
}
generated quantities {
  array[ntotal_lb] real y_lb_pred = normal_rng(mu_lb[x_lb_idx], sig_lb[x_lb_idx]);
  array[ntotal_bl] real y_bl_pred = normal_rng(mu_bl[x_bl_idx], sig_bl[x_bl_idx]);
}
