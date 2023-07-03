
tmb_model <- "
// linear regression
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_MATRIX(m_x);
DATA_VECTOR(v_y);

// parameters:
PARAMETER_VECTOR(v_b); // slope
PARAMETER(log_sigma); // log(residual SD)
// we fit sigma on a log scale to keep it > 0

// procedures: (transformed parameters)
Type sigma = exp(log_sigma);

Type nll = 0.0; // initialize negative log likelihood

nll = -sum(dnorm(v_y, m_x * v_b, sigma, true));

return nll;
}"
write(tmb_model, file = "code/regression.cpp")

library(TMB)
compile("code/regression.cpp")
dyn.load(dynlib("code/regression"))

set.seed(122)
df0 <- tibble(x1 = runif(20, 1, 10), x2 = runif(20, 1, 10)) %>% 
  mutate(y = rnorm(20, mean = 1.8 + 2.4 * x1 - 3 * x2, sd = exp(0.3)))
fm <- y ~ .

m_x <- data.matrix(model.matrix(fm, data = df0))

obj <- MakeADFun(
  data = list(m_x = m_x, v_y = y), 
  parameters = list(v_b = rep(0, 3), log_sigma = 0),
  DLL = "regression")
opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
rep <- sdreport(obj)
