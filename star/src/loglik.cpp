#include <Rcpp.h>
using namespace Rcpp;

// Bernstein polynomial of degree K with coefficients coef (length K+1) at w in (0,1)
double bern_eval(double w, const NumericVector& coef) {
  int K = coef.size() - 1;
  double out = 0.0;
  for (int k = 0; k <= K; k++) {
    out += coef[k] * R::choose(K, k) * std::pow(w, k) * std::pow(1.0 - w, K - k);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector bern_eval_vec(NumericVector w, NumericVector coef) {
  int n = w.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = bern_eval(w[i], coef);
  return out;
}

// tau(w) = min(w,1-w)/max(w,1-w)
inline double tau_fun(double w) {
  double mn = std::min(w, 1.0 - w), mx = std::max(w, 1.0 - w);
  return mn / mx;
}

// sigma(w) = 1 + chibar * tau(w)^(1+nu(w)),  nu = exp(Bernstein(z))
// [[Rcpp::export]]
NumericVector sigma_cpp(NumericVector w, double chibar, NumericVector z) {
  int n = w.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    double nu = std::exp(bern_eval(w[i], z));
    out[i] = 1.0 + chibar * std::pow(tau_fun(w[i]), 1.0 + nu);
  }
  return out;
}

// rho(w) = ((w/u1)^p + ((1-w)/u2)^p)^(-1/p)
// [[Rcpp::export]]
NumericVector rho_cpp(NumericVector w, double u1, double u2, double p) {
  int n = w.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = std::pow(std::pow(w[i] / u1, p) + std::pow((1.0 - w[i]) / u2, p), -1.0 / p);
  }
  return out;
}

// log f_W(w) for the Bernstein (Beta mixture) angular density with weights phi
inline double log_fW(double w, const NumericVector& phi) {
  int K = phi.size();
  double dens = 0.0;
  for (int k = 1; k <= K; k++) {
    dens += phi[k - 1] * R::dbeta(w, (double)k, (double)(K - k + 1), 0);
  }
  return std::log(dens);
}

// [[Rcpp::export]]
NumericVector fW_cpp(NumericVector w, NumericVector phi) {
  int n = w.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = std::exp(log_fW(w[i], phi));
  return out;
}

// Exact spliced log-likelihood, d = 2.
// theta components: phi (K_W), alpha, b (K_lam+1), u1, u2, p, chibar, z (K_nu+1)
// data: radii r, angles w.
// log f(e_i) = log f_W(w_i) + log f_{R|W}(r_i|w_i) - log r_i
// [[Rcpp::export]]
double loglik_cpp(NumericVector r, NumericVector w,
                  NumericVector phi, double alpha, NumericVector b,
                  double u1, double u2, double p,
                  double chibar, NumericVector z) {
  int n = r.size();
  double ll = 0.0;
  for (int i = 0; i < n; i++) {
    double lam = std::exp(bern_eval(w[i], b));
    double rho = std::pow(std::pow(w[i] / u1, p) + std::pow((1.0 - w[i]) / u2, p), -1.0 / p);
    ll += log_fW(w[i], phi) - std::log(r[i]);
    if (r[i] <= rho) {
      ll += R::dgamma(r[i], alpha, 1.0 / lam, 1);  // body: Gamma(alpha, rate lam)
    } else {
      double nu = std::exp(bern_eval(w[i], z));
      double sig = 1.0 + chibar * std::pow(tau_fun(w[i]), 1.0 + nu);
      double A = R::pgamma(rho, alpha, 1.0 / lam, 0, 1);  // log upper tail of body at rho
      ll += A - std::log(sig) - (r[i] - rho) / sig;
    }
    if (!std::isfinite(ll)) return -1e300;
  }
  return ll;
}

// Joint survivor P(E1>t1, E2>t2) by Gauss-Legendre quadrature of
// int f_W(w) * Sbar( max(t1/w, t2/(1-w)) | w ) dw, with the full spliced survivor.
// nodes and weights are supplied on (0,1).
// [[Rcpp::export]]
double joint_surv_cpp(double t1, double t2,
                      NumericVector nodes, NumericVector wts,
                      NumericVector phi, double alpha, NumericVector b,
                      double u1, double u2, double p,
                      double chibar, NumericVector z) {
  int M = nodes.size();
  double out = 0.0;
  for (int m = 0; m < M; m++) {
    double w = nodes[m];
    double lev = std::max(t1 / w, t2 / (1.0 - w));
    double lam = std::exp(bern_eval(w, b));
    double rho = std::pow(std::pow(w / u1, p) + std::pow((1.0 - w) / u2, p), -1.0 / p);
    double S;
    if (lev <= rho) {
      S = R::pgamma(lev, alpha, 1.0 / lam, 0, 0);
    } else {
      double nu = std::exp(bern_eval(w, z));
      double sig = 1.0 + chibar * std::pow(tau_fun(w), 1.0 + nu);
      S = R::pgamma(rho, alpha, 1.0 / lam, 0, 0) * std::exp(-(lev - rho) / sig);
    }
    out += wts[m] * std::exp(log_fW(w, phi)) * S;
  }
  return out;
}
