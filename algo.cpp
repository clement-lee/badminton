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
const double llik_model(const NumericVector y1, const NumericVector y2, const IntegerVector x1, const IntegerVector x2, const NumericVector gamma, const double beta) {
  // log-likelihood of model
  // y1 & y2 are number of games won by players 1 & 2, respectively
  // x1 & x2 are the indices of players 1 & 2 in gamma
  // gamma is the vector of log-"strengths"
  const int n = y1.size();
  if (y2.size() != n || x1.size() != n || x2.size() != n) {
    stop("llik_model: y1, y2, x1 and x2 have to be of the same length.");
  }
  const NumericVector
    gamma1 = gamma[x1],
    gamma2 = gamma[x2],
    r1 = 1.0 / (1.0 + exp(gamma2 - gamma1)),
    r2 = 1.0 / (1.0 + exp(gamma1 - gamma2)),
    p1 = pbeta(r1, beta, beta, true, true),
    p2 = pbeta(r2, beta, beta, true, true),
    l = y1 * p1 + y2 * p2;
  double llik = sum(l);
  if (llik != llik) {
    llik = -INFINITY;
  }
  return llik;
}

// [[Rcpp::export]]
const double lprior_model(const NumericVector gamma, const double beta, const double mu, const double sigma, const double a_beta = 1.0, const double b_beta = 0.001, double m_mu = 0.0, const double s_mu = 100.0, const double a_sigma = 1.0, const double b_sigma = 0.001) {
  double lprior =
    sum(dnorm(gamma, mu, sigma, true)) +
    dgamma(tv(beta), a_beta, 1.0 / b_beta, true)[0] +
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
              const double beta,
              const double mu,
              const double sigma,
              const double a_beta = 1.0,
              const double b_beta = 0.001,
              const double m_mu = 0.0,
              const double s_mu = 100.0,
              const double a_sigma = 1.0,
              const double b_sigma = 0.001,
              double sd_gamma = 0.001,
              double sd_beta = 0.001,
              double sd_mu = 0.001,
              double sd_sigma = 0.001,
              const int N = 40000,
              const int thin = 1,
              const int burnin = 10000,
              const int print_freq = 100) {
  const NumericVector hypers = NumericVector::create(beta, sigma, a_beta, b_beta, s_mu, a_sigma, b_sigma);
  if (is_true(any(hypers <= 0.0))) {
    stop("mh_model: initial value of beta & other hyperparameters must be positive.");
  }
  NumericVector gamma_init = rnorm(m, mu, sigma), gamma_curr = clone(gamma_init), gamma_prop = clone(gamma_curr), sds_gamma(m, sd_gamma);
  double beta_curr = beta, beta_prop,
    mu_curr = mu, mu_prop,
    sigma_curr = sigma, sigma_prop;
  auto lpost = [y1, y2, x1, x2, a_beta, b_beta, m_mu, s_mu, a_sigma, b_sigma](const NumericVector gamma, const double beta, const double mu, const double sigma) {
    double lpost = llik_model(y1, y2, x1, x2, gamma, beta) + lprior_model(gamma, beta, mu, sigma, a_beta, b_beta, m_mu, s_mu, a_sigma, b_sigma);
    if (lpost != lpost) {
      lpost = -INFINITY;
    }
    return lpost;
  };
  double lpost_curr = lpost(gamma_curr, beta_curr, mu_curr, sigma_curr), lpost_prop;
  NumericMatrix gamma_mat(N, m);
  NumericVector beta_vec(N), mu_vec(N), sigma_vec(N);
  // run
  int i, j;
  for (i = 0; i < N * thin + burnin; i++) {
    // update gamma
    for (j = 0; j < m; j++) {
      gamma_prop = clone(gamma_curr);
      gamma_prop[j] = rnorm(1, gamma_curr[j], sds_gamma[j])[0];
      lpost_prop = lpost(gamma_prop, beta_curr, mu_curr, sigma_curr);
      update(gamma_curr[j], gamma_prop[j], lpost_curr, lpost_prop, sds_gamma[j], i, burnin);
    }
    // update beta
    beta_prop = rnorm(1, beta_curr, sd_beta)[0];
    lpost_prop = lpost(gamma_curr, beta_prop, mu_curr, sigma_curr);
    update(beta_curr, beta_prop, lpost_curr, lpost_prop, sd_beta, i, burnin);
    // update mu
    mu_prop = rnorm(1, mu_curr, sd_mu)[0];
    lpost_prop = lpost(gamma_curr, beta_curr, mu_prop, sigma_curr);
    update(mu_curr, mu_prop, lpost_curr, lpost_prop, sd_mu, i, burnin);
    // update sigma
    sigma_prop = rnorm(1, sigma_curr, sd_sigma)[0];
    lpost_prop = lpost(gamma_curr, beta_curr, mu_curr, sigma_prop);
    update(sigma_curr, sigma_prop, lpost_curr, lpost_prop, sd_sigma, i, burnin);
    // print & save
    if ((i + 1) % print_freq == 0) {
      Rcout << "Iteration " << i + 1 << endl;
    }
    if (i >= burnin && (i - burnin + 1) % thin == 0) {
      j = (i - burnin + 1) / thin - 1;
      gamma_mat(j, _) = gamma_curr;
      beta_vec[j] = beta_curr;
      mu_vec[j] = mu_curr;
      sigma_vec[j] = sigma_curr;
    }
  }
  // output
  DataFrame par = DataFrame::create(Named("beta") = beta_vec,
                                    Named("mu") = mu_vec,
                                    Named("sigma") = sigma_vec),
    initial = DataFrame::create(Named("beta") = tv(beta),
                                Named("mu") = tv(mu),
                                Named("sigma") = tv(sigma)),
    hyper = DataFrame::create(Named("a_beta") = tv(a_beta),
                              Named("b_beta") = tv(b_beta),
                              Named("m_mu") = tv(m_mu),
                              Named("s_mu") = tv(s_mu),
                              Named("a_sigma") = tv(a_sigma),
                              Named("b_sigma") = tv(b_sigma)),
    sds = DataFrame::create(Named("sd_beta") = tv(sd_beta),
                            Named("sd_mu") = tv(sd_mu),
                            Named("sd_sigma") = tv(sd_sigma)),
    scalars = DataFrame::create(Named("N") = IntegerVector::create(N),
                                Named("thin") = IntegerVector::create(thin),
                                Named("burnin") = IntegerVector::create(burnin),
                                Named("print_freq") = IntegerVector::create(print_freq));
  List output = List::create(Named("gamma_par") = gamma_mat, Named("gamma_initial") = gamma_init, Named("gamma_sds") = sds_gamma, Named("par") = par, Named("initial") = initial, Named("hyper") = hyper, Named("sds") = sds, Named("scalars") = scalars);
  return output;
}
