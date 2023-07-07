
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
PARAMETER(lambda);

// procedures: (transformed parameters)
Type sigma = exp(log_sigma);
Type theta = exp(log_theta);

// report
ADREPORT(b);
ADREPORT(sigma);
ADREPORT(theta);
ADREPORT(lambda);

// negative log likelihood
Type nll = 0;
matrix<Type> W = exp(-theta * D.array()); 
//vector<Type> sum_W = W.rowwise().sum();
//matrix<Type> Q(W);

//for (int i = 0; i < W.rows(); i++)
//  Q.row(i) = W.row(i) / sum_W(i);

vector<Type> u = W * y;
nll = -sum(dnorm(y, X * b + lambda * u, sigma, true));

return nll;
}
