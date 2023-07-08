
// linear regression
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_VECTOR(y);
DATA_MATRIX(X);
DATA_MATRIX(D);

// parameters:
PARAMETER_VECTOR(b);
PARAMETER(log_sigma);
PARAMETER(log_theta);

using namespace density;

// transformed parameters
Type sigma = exp(log_sigma);
Type theta = exp(log_theta);

// report
ADREPORT(b);
ADREPORT(sigma);
ADREPORT(theta);

// negative log likelihood
matrix<Type> R = exp(- D.array() / theta); 
matrix<Type> S = sigma * R.array(); 
MVNORM_t<Type> dmnorm(S);
vector<Type> u = y - X * b;

parallel_accumulator<Type> nll(this);
nll += dmnorm(u);

return nll;
}
