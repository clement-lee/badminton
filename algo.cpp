#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

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
const double llik_model(const NumericVector y1, const NumericVector y2, const IntegerVector x1, const IntegerVector x2, const NumericVector gamma, const double eta) {
  // log-likelihood of model
  // y1 & y2 are number of games won by players 1 & 2, respectively
  // x1 & x2 are the indices of players 1 & 2 in gamma
  // gamma is the vector of log-"strengths"
  const int n = y1.size();
  if (y2.size() != n || x1.size() != n || x2.size() != n) {
    stop("llik_model: y1, y2, x1 and x2 have to be of the same length.");
  }
  double llik;
  const NumericVector
    gamma1 = gamma[x1],
    gamma2 = gamma[x2],
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
const double lprior_model(const NumericVector gamma, const double eta, const double mu, const double sigma, const double a_eta = 1.0, const double b_eta = 1.0, const double m_mu = 0.0, const double s_mu = 100.0, const double a_sigma = 1.0, const double b_sigma = 0.001) {
  double lprior =
    sum(dnorm(gamma, mu, sigma, true)) +
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
              const double eta = 0.5,
              const double mu = 0.0,
              const double sigma = 1.0,
              const double a_eta = 1.0,
              const double b_eta = 1.0,
              const double m_mu = 0.0,
              const double s_mu = 100.0,
              const double a_sigma = 1.0,
              const double b_sigma = 0.001,
              double sd_gamma = 0.001,
              double sd_eta = 0.001,
              double sd_mu = 0.001,
              double sd_sigma = 0.001,
              const int N = 40000,
              const int thin = 1,
              const int burnin = 10000,
              const int print_freq = 100,
              const int update_eta = true) {
  const NumericVector hypers = NumericVector::create(eta, sigma, a_eta, b_eta, s_mu, a_sigma, b_sigma);
  if (is_true(any(hypers <= 0.0))) {
    stop("mh_model: initial value of eta & other hyperparameters must be positive.");
  }
  NumericVector gamma_init = rnorm(m, mu, sigma), gamma_curr = clone(gamma_init), gamma_prop = clone(gamma_curr), sds_gamma(m, sd_gamma);
  gamma_init[0] = 0.0;
  gamma_curr[0] = 0.0; // for identifiability?
  gamma_prop[0] = 0.0;
  double eta_curr = eta, eta_prop,
    mu_curr = mu, mu_prop,
    sigma_curr = sigma, sigma_prop;
  auto lpost = [y1, y2, x1, x2, a_eta, b_eta, m_mu, s_mu, a_sigma, b_sigma](const NumericVector gamma, const double eta, const double mu, const double sigma) {
    double lpost = llik_model(y1, y2, x1, x2, gamma, eta) + lprior_model(gamma, eta, mu, sigma, a_eta, b_eta, m_mu, s_mu, a_sigma, b_sigma);
    if (lpost != lpost) {
      lpost = -INFINITY;
    }
    return lpost;
  };
  double lpost_curr = lpost(gamma_curr, eta_curr, mu_curr, sigma_curr), lpost_prop;
  NumericMatrix gamma_mat(N, m);
  NumericVector eta_vec(N), mu_vec(N), sigma_vec(N), lpost_vec(N);
  // run
  int i, j;
  for (i = 0; i < N * thin + burnin; i++) {
    // update gamma
    for (j = 1; j < m; j++) { // [0] not updated
      gamma_prop = clone(gamma_curr);
      gamma_prop[j] = rnorm(1, gamma_curr[j], sds_gamma[j])[0];
      lpost_prop = lpost(gamma_prop, eta_curr, mu_curr, sigma_curr);
      update(gamma_curr[j], gamma_prop[j], lpost_curr, lpost_prop, sds_gamma[j], i, burnin);
    }
    // update eta
    if (update_eta) {
      eta_prop = rnorm(1, eta_curr, sd_eta)[0];
      lpost_prop = lpost(gamma_curr, eta_prop, mu_curr, sigma_curr);
      update(eta_curr, eta_prop, lpost_curr, lpost_prop, sd_eta, i, burnin);
    }
    // update mu
    mu_prop = rnorm(1, mu_curr, sd_mu)[0];
    lpost_prop = lpost(gamma_curr, eta_curr, mu_prop, sigma_curr);
    update(mu_curr, mu_prop, lpost_curr, lpost_prop, sd_mu, i, burnin);
    // update sigma
    sigma_prop = rnorm(1, sigma_curr, sd_sigma)[0];
    lpost_prop = lpost(gamma_curr, eta_curr, mu_curr, sigma_prop);
    update(sigma_curr, sigma_prop, lpost_curr, lpost_prop, sd_sigma, i, burnin);
    // print & save
    if ((i + 1) % print_freq == 0) {
      Rcout << "Iteration " << i + 1 << ": lpost = " << lpost_curr << endl;
      Rcout << "eta: " << eta_curr << " (" << sd_eta << ")" << endl;
      Rcout << "mu: " << mu_curr << " (" << sd_mu << ")" << endl;
      Rcout << "sigma: " << sigma_curr << " (" << sd_sigma << ")" << endl;
      Rcout << endl;
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      gamma_mat(j, _) = gamma_curr;
      eta_vec[j] = eta_curr;
      mu_vec[j] = mu_curr;
      sigma_vec[j] = sigma_curr;
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
    sds = DataFrame::create(Named("sd_eta") = tv(sd_eta),
                            Named("sd_mu") = tv(sd_mu),
                            Named("sd_sigma") = tv(sd_sigma)),
    scalars = DataFrame::create(Named("N") = IntegerVector::create(N),
                                Named("thin") = IntegerVector::create(thin),
                                Named("burnin") = IntegerVector::create(burnin),
                                Named("print_freq") = IntegerVector::create(print_freq));
  List output = List::create(Named("gamma_par") = gamma_mat, Named("gamma_initial") = gamma_init, Named("gamma_sds") = sds_gamma, Named("lpost") = lpost_df, Named("par") = par, Named("initial") = initial, Named("hyper") = hyper, Named("sds") = sds, Named("scalars") = scalars);
  return output;
}
