// global psychophysics model
// reference parameters are estimated as constant sum
// The variances are estimated depending on sigindex input

functions {
  real db_lambert(real std_b) {
    return 10 * log10(std_b * pi() / (10 ^ -6));}
  real db_spl(real std_l) {
    return 20 * log10(std_l / (2 * (10 ^ -5)));}
  real db_inv_lambert(real db_b) {
    return (10 ^ (db_b / 10)) * (10 ^ -6) / pi();}
  real db_inv_spl(real db_l) {
    return (10 ^ (db_l / 20)) * 2 * (10 ^ -5);}
  real loud_to_bright(real x_loud, real alpha_l, real alpha_b,
                      real beta_l, real beta_b, real omega_1, real const_lb) {
    return db_lambert(
      pow(
        inv(alpha_b) * 
          (omega_1 * alpha_l * pow(db_inv_spl(x_loud), beta_l) -
           const_lb),
        inv(beta_b))
    );
  }
  real bright_to_loud(real x_bright, real alpha_l, real alpha_b,
                      real beta_l, real beta_b, real omega_1, real const_bl) {
    return db_spl(
       pow(
        inv(alpha_l) *
          (omega_1 * alpha_b * pow(db_inv_lambert(x_bright), beta_b) -
           const_bl), 
         inv(beta_l))
    );
  }
}
data {
  int<lower=1> ntotal; // total ntrials
  int<lower=1> ncond;
  int<lower=1> nsig;
  array[ncond] int<lower=1, upper=2> task; // 1: loud_bright 2: bright_loud
  vector<lower=1>[ncond] std;
  array[ntotal] int<lower=1> idx;
  array[ntotal] int<lower=1> sigidx;
  array[ntotal] real<lower=0> tgt;
  int<lower=0, upper=1> onlyprior;
}
transformed data {
  int<lower=1, upper=1> alpha_b = 1;
}
parameters {
  real<lower=0> alpha_l;
  real<lower=0, upper=1> beta_l;
  real<lower=0, upper=1> beta_b;
  real<lower=0> omega_1;
  vector<lower=0>[nsig] sig;
  real const_lb;
  real const_bl;
}
transformed parameters {
  array[ncond] real mu; 
  for (i in 1:ncond) {
    if (task[i] == 1) {
      mu[i] = loud_to_bright(std[i], alpha_l, alpha_b, beta_l,
                             beta_b, omega_1, const_lb);
    } else if (task[i] == 2) {
      mu[i] = bright_to_loud(std[i], alpha_l, alpha_b, beta_l,
                             beta_b, omega_1, const_bl);

    }
  }
}
model {
  target += lognormal_lpdf(alpha_l | log(10), 1);
  target += beta_lpdf(beta_b | 6, 6);
  target += beta_lpdf(beta_l | 6, 6);
  target += normal_lpdf(const_lb | 0, 1);
  target += normal_lpdf(const_bl | 0, 1);
  target += normal_lpdf(omega_1 | 1, 0.3);
  target += normal_lpdf(sig | 5, 3) - normal_lccdf(0 | 5, 3);
  if (!onlyprior) {
    target += normal_lpdf(tgt | to_vector(mu)[idx],
                                to_vector(sig[sigidx]));
  }
}
generated quantities {
  array[ntotal] real tgt_pred = normal_rng(
    to_vector(mu)[idx],
    to_vector(sig[sigidx]));
}
