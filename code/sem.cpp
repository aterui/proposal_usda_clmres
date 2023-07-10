#include <TMB.hpp>

enum valid_family {
  gaussian_family = 0,
  poisson_family = 1,
};

enum valid_link {
  identity_link = 0,
  log_link = 1
};

template<class Type>
Type inv_link(Type eta, int link) {
  Type ans;
  switch (link) {
  case identity_link:
    ans = eta;
    break;
  case log_link:
    ans = exp(eta);
    break;
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(D);
  DATA_INTEGER(link);
  DATA_INTEGER(family);
  
  // parameters:
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(u);
  PARAMETER(log_sigma0);
  PARAMETER(log_sigma);
  PARAMETER(log_theta);
  
  // transformed parameters
  Type sigma0 = exp(log_sigma0);
  Type sigma = exp(log_sigma);
  Type theta = exp(log_theta);
  
  // report
  ADREPORT(b);
  ADREPORT(sigma0);
  ADREPORT(sigma);
  ADREPORT(theta);
  
  // negative log likelihood
  Type jnll = 0;
  
  matrix<Type> R = exp(- theta * D.array());
  density::MVNORM_t<Type> dmnorm(R);
  density::SCALE_t<density::MVNORM_t<Type>> scl_dmnorm = density::SCALE(dmnorm, sigma);
  
  vector<Type> eta = X * b + u;
  jnll += scl_dmnorm(u);
  
  vector<Type> mu(eta.size());
  for (int i = 0; i < eta.size(); i++)
    mu(i) = inv_link(eta(i), link);
  
  Type nll = 0;
  for (int i = 0; i < eta.size(); i++) {
    switch (family) {
      case gaussian_family:
        nll = dnorm(y(i), mu(i), sigma0, true);
        break;
      case poisson_family:
        nll = dpois(y(i), mu(i), true);
        break;
      default:
        error("Family not implemented!");
    }//switch end
    
    jnll -= nll;
  }
  
  return jnll;
}