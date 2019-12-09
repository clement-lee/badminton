// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

// 00) prelim
void update(double & par_curr, const double par_prop, double & lpost_curr, const double lpost_prop, double & s, const int i, const int burnin, const double factor = 30.0) {
  const bool accept_reject = log(runif(1)[0]) < lpost_prop - lpost_curr;
  par_curr = accept_reject ? par_prop : par_curr;
  lpost_curr = accept_reject ? lpost_prop : lpost_curr;
  if (i < burnin) {
    s = sqrt(s * s + (accept_reject ? 3.0 : (-1.0)) * s * s / factor / sqrt(i + 1.0));
  }
}

const NumericVector tv(const double x) {
  return NumericVector::create(x);
}



// 01) model
// [[Rcpp::export]]
const double llik_model(const NumericVector y1, const NumericVector y2, const IntegerVector x1, const IntegerVector x2, const vec gamma, const double eta) {
  // log-likelihood of model
  // y1 & y2 are number of games won by players 1 & 2, respectively
  // x1 & x2 are the indices of players 1 & 2 in gamma
  // gamma is the vector of log-"strengths"
  const int n = y1.size();
  if (y2.size() != n || x1.size() != n || x2.size() != n) {
    stop("llik_model: y1, y2, x1 and x2 have to be of the same length.");
  }
  double llik;
  NumericVector gamma1 = wrap(gamma), gamma2 = wrap(gamma);
  gamma1 = gamma1[x1];
  gamma2 = gamma2[x2];
  const NumericVector
    r1 = 1.0 / (1.0 + exp(gamma2 - gamma1)),
    r2 = 1.0 / (1.0 + exp(gamma1 - gamma2));
  NumericVector p1(n, 0.0), p2(n, 0.0);
  if (eta > 1.0 || eta < 0.0) {
    llik = -INFINITY;
  }
  else {
    if (eta == 1.0) {
      p1 = log(r1);
      p2 = log(r2);
    }
    else if (eta == 0.0) {
      p1 = pnorm(gamma1 - gamma2, 0.0, 1.0, true, true);
      p2 = pnorm(gamma2 - gamma1, 0.0, 1.0, true, true);
    }
    else {
      const double beta = 1.0 / eta;
      p1 = pbeta(r1, beta, beta, true, true);
      p2 = pbeta(r2, beta, beta, true, true);
    }
    const NumericVector l = y1 * p1 + y2 * p2;
    llik = sum(l);
  }
  if (llik != llik) {
    llik = -INFINITY;
  }
  return llik;
}

// [[Rcpp::export]]
const double lprior_model(const vec gamma, const double eta, const double mu, const double sigma, const double a_eta = 1.0, const double b_eta = 1.0, const double m_mu = 0.0, const double s_mu = 100.0, const double a_sigma = 1.0, const double b_sigma = 0.001) {
  const NumericVector gamma0 = wrap(gamma);
  double lprior =
    sum(dnorm(gamma0, mu, sigma, true)) +
    dbeta(tv(eta), a_eta, b_eta, true)[0] +
    dnorm(tv(mu), m_mu, s_mu, true)[0] +
    dgamma(tv(sigma), a_sigma, 1.0 / b_sigma, true)[0];
  if (lprior != lprior) {
    lprior = -INFINITY;
  }
  return lprior;
}

// [[Rcpp::export]]
List mh_model(const NumericVector y1,
              const NumericVector y2,
              const IntegerVector x1,
              const IntegerVector x2,
              const int m,
              const mat M0_par,
              const double eta = 0.5,
              const double mu = 0.0,
              const double sigma = 1.0,
              const double a_eta = 1.0,
              const double b_eta = 1.0,
              const double m_mu = 0.0,
              const double s_mu = 100.0,
              const double a_sigma = 1.0,
              const double b_sigma = 0.001,
              const int N = 40000,
              const int thin = 1,
              const int burnin = 10000,
              const int print_freq = 100) {
  if (M0_par.n_rows != 3 || M0_par.n_cols != 3) {
    stop("mh_model: M0_par has to be 3x3 matrix.");
  }
  const NumericVector hypers = NumericVector::create(eta, sigma, a_eta, b_eta, s_mu, a_sigma, b_sigma);
  if (is_true(any(hypers <= 0.0))) {
    stop("mh_model: initial value of eta & other hyperparameters must be positive.");
  }
  NumericVector gamma_init = rnorm(m, mu, sigma);
  gamma_init[0] = 0.0; // for identifiability?
  vec gamma_curr = as<vec>(gamma_init), gamma_prop = as<vec>(gamma_init), par_curr(3), par_prop(3);
  par_curr[0] = eta;
  par_curr[1] = mu;
  par_curr[2] = sigma;
  auto lpost = [y1, y2, x1, x2, a_eta, b_eta, m_mu, s_mu, a_sigma, b_sigma](const vec gamma, const vec par) {
    double lpost = llik_model(y1, y2, x1, x2, gamma, par[0]) + lprior_model(gamma, par[0], par[1], par[2], a_eta, b_eta, m_mu, s_mu, a_sigma, b_sigma);
    if (lpost != lpost) {
      lpost = -INFINITY;
    }
    return lpost;
  };
  double lpost_curr = lpost(gamma_curr, par_curr), lpost_prop;
  mat gamma_mat(N, m);
  NumericVector eta_vec(N), mu_vec(N), sigma_vec(N), lpost_vec(N);
  running_stat_vec<vec> stats_par(true);
  mat M1_par(3, 3);
  running_stat_vec<vec> stats_gamma(true);
  mat M0_gamma(m, m, fill::eye), M1_gamma(m, m);
  M0_gamma(0, 0) = 0.0;
  // run
  int i, j;
  for (i = 0; i < N * thin + burnin; i++) {
    // update gamma
    if (i < 0.1 * burnin || runif(1)[0] < 0.05) {
      gamma_prop = mvnrnd(gamma_curr, 0.01 / (m + 0.0) * M0_gamma);
    }
    else {
      gamma_prop = mvnrnd(gamma_curr, (2.34 * 2.34 / (m + 0.0)) * M1_gamma);
    }
    lpost_prop = lpost(gamma_prop, par_curr);
    if (log(runif(1)[0]) < lpost_prop - lpost_curr) {
      gamma_curr = gamma_prop;
      lpost_curr = lpost_prop;
    }
    // update eta, mu, sigma
    if (i < 0.1 * burnin || runif(1)[0] < 0.05) {
      par_prop = mvnrnd(par_curr, 0.01 / 3.0 * M0_par);
    }
    else {
      par_prop = mvnrnd(par_curr, (2.34 * 2.34 / 3.0) * M1_par);
    }
    lpost_prop = lpost(gamma_curr, par_prop);
    if (log(runif(1)[0]) < lpost_prop - lpost_curr) {
      par_curr = par_prop;
      lpost_curr = lpost_prop;
    }
    // adaptive
    if (i < burnin) {
      stats_gamma(gamma_curr);
      M1_gamma = stats_gamma.cov();
      stats_par(par_curr);
      M1_par = stats_par.cov();
    }
    // print & save
    if ((i + 1) % print_freq == 0) {
      Rcout << "Iteration " << i + 1 << ": lpost = " << lpost_curr << endl;
      Rcout << "eta: " << par_curr[0] << endl;
      Rcout << "mu: " << par_curr[1] << endl;
      Rcout << "sigma: " << par_curr[2] << endl;
      Rcout << "M1_par: " << M1_par << endl;
      Rcout << endl;
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      gamma_mat.row(j) = gamma_curr.t();
      eta_vec[j] = par_curr[0];
      mu_vec[j] = par_curr[1];
      sigma_vec[j] = par_curr[2];
      lpost_vec[j] = lpost_curr;
    }
  }
  // output
  DataFrame lpost_df = DataFrame::create(Named("iteration") = seq_len(N), Named("lpost") = lpost_vec),
    par = DataFrame::create(Named("eta") = eta_vec,
                            Named("mu") = mu_vec,
                            Named("sigma") = sigma_vec),
    initial = DataFrame::create(Named("eta") = tv(eta),
                                Named("mu") = tv(mu),
                                Named("sigma") = tv(sigma)),
    hyper = DataFrame::create(Named("a_eta") = tv(a_eta),
                              Named("b_eta") = tv(b_eta),
                              Named("m_mu") = tv(m_mu),
                              Named("s_mu") = tv(s_mu),
                              Named("a_sigma") = tv(a_sigma),
                              Named("b_sigma") = tv(b_sigma)),
    sds = DataFrame::create(Named("sd_eta") = tv(M1_par(0,0)),
                            Named("sd_mu") = tv(M1_par(1,1)),
                            Named("sd_sigma") = tv(M1_par(2,2))),
    scalars = DataFrame::create(Named("N") = IntegerVector::create(N),
                                Named("thin") = IntegerVector::create(thin),
                                Named("burnin") = IntegerVector::create(burnin),
                                Named("print_freq") = IntegerVector::create(print_freq));
  List output = List::create(Named("gamma_par") = gamma_mat, Named("gamma_initial") = gamma_init, Named("lpost") = lpost_df, Named("par") = par, Named("initial") = initial, Named("hyper") = hyper, Named("sds") = sds, Named("scalars") = scalars);
  return output;
}
