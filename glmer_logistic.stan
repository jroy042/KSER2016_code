// Stan model for mixed effects logistic regression
// with random effects for subject and item

// the code is written assuming that both subjects and items
// have at least one slope.  It is not written for models with
// only one grouping factor or intercept-only models.  These
// would require a few simple modifications.  See the online
// supplement for sample R code.

// the model matrix x should have columns with the following
// order and then be converted to CSR format:
//   fixed effects
//   subject effects ordered by subject
//     within-subject, the effects should be in the same order
//   item effects ordered by item
//     within-item, the effects should be in the same order

// the binary response y should be passed as integer 0/1

// the (reparameterized) priors are:
//   fixed effects
//     beta ~ student_t(5,0,sigma_beta)
//     sigma_beta ~ half_normal(0,1)
//   random effects
//     sigma_group ~ half_normal(0,1)
//     omega_group ~ lkj_corr(2)
//     gamma_group ~ mult_normal(0,Sigma_group)
//       where Sigma_group = diag(sigma_group)
//         * omega_group * t(diag(sigma_group))

functions {
  // this function takes a vector and turns it into
  // a matrix with R rows and C columns, loading
  // the values by row
  matrix vec_to_matrix_by_row(vector v, int R, int C){
    matrix[R,C] m;
    for(r in 1:R) m[r] = v[(C*(r-1)+1):(C*r)]';
    return m;
  }
}

data {
  int<lower=2> N;  // number of observations
  int<lower=2> S;  // number of subjects
  int<lower=2> I;  // number of items

  int<lower=1> P;  // number of fixed effects
  int<lower=1,upper=P> QS;  // number of subject effects
  int<lower=1,upper=P> QI;  // number of item effects

  int<lower=0,upper=1> y[N];  // binary response

  // sparse model matrix (CSR)
  int<lower=1> nz;  // number of non-zero elements in x
  vector[nz] x_w;  // non-zero elements in x
  int x_v[nz];  // column indices for x_w
  int x_u[N+1];  // row-start indices for x
}

transformed data {
  int K;  // number of columns in x
  int SF;  // first subject effect column in x
  int SL;  // last subject effect column in x
  int IF;  // first item effect column in x
  int IL;  // last item effect column in x

  K = P + S * QS + I * QI;
  SF = P + 1;
  SL = P + S * QS;
  IF = SL + 1;
  IL = SL + I * QI;
}

parameters {
  vector[P] beta_raw;
  vector<lower=0>[P] tau_beta;
  real<lower=0> sigma_beta;

  matrix[QS,S] gamma_subj_raw;
  vector<lower=0>[QS] sigma_subj;  // subject effect SDs
  cholesky_factor_corr[QS] omega_subj_raw;

  matrix[QI,I] gamma_item_raw;
  vector<lower=0>[QI] sigma_item;  // item effect SDs
  cholesky_factor_corr[QI] omega_item_raw;
}

transformed parameters {
  vector[K] coef;  // all coefficients
  vector[N] y_hat;  // predicted log-odds

  // transform fixed effects
  for(p in 1:P)
    coef[p] = beta_raw[p] * sigma_beta / sqrt(tau_beta[p]);

  // transform subject effects
  coef[SF:SL]
    = to_vector(rep_matrix(sigma_subj,S)
      .* (omega_subj_raw * gamma_subj_raw));

  // transform item effects
  coef[IF:IL]
    = to_vector(rep_matrix(sigma_item,I)
      .* (omega_item_raw * gamma_item_raw));

  // y_hat = x * coef
  y_hat = csr_matrix_times_vector(N,K,x_w,x_v,x_u,coef);
}

model {
  beta_raw ~ normal(0,1);
  tau_beta ~ gamma(2.5,2.5);
  sigma_beta ~ normal(0,1);

  to_vector(gamma_subj_raw) ~ normal(0,1);
  sigma_subj ~ normal(0,1);
  omega_subj_raw ~ lkj_corr_cholesky(2);

  to_vector(gamma_item_raw) ~ normal(0,1);
  sigma_item ~ normal(0,1);
  omega_item_raw ~ lkj_corr_cholesky(2);

  y ~ bernoulli_logit(y_hat);  // logistic model defined
}

generated quantities {
  vector[P] beta;  // fixed effects
  matrix[S,QS] gamma_subj;  // subject effects
  matrix[I,QI] gamma_item;  // item effects
  matrix[QS,QS] omega_subj;  // correlation in subject effects
  matrix[QI,QI] omega_item;  // correlation in item effects

  beta = coef[1:P];
  gamma_subj = vec_to_matrix_by_row(coef[SF:SL],S,QS);
  gamma_item = vec_to_matrix_by_row(coef[IF:IL],I,QI);
  omega_subj = tcrossprod(omega_subj_raw);
  omega_item = tcrossprod(omega_item_raw);
}
