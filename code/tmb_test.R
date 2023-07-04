
tmb_model <- "
// linear regression
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_MATRIX(m_x);
DATA_MATRIX(m_d);
DATA_VECTOR(v_y);

// parameters:
PARAMETER_VECTOR(v_b); // slope
PARAMETER_VECTOR(v_z);
PARAMETER(log_sigma); // log(residual SD)
PARAMETER(log_sigma_z); // log(residual SD)

// we fit sigma on a log scale to keep it > 0

// procedures: (transformed parameters)
Type sigma = exp(log_sigma);
Type sigma_z = exp(log_sigma_z);
Type nll = 0;
matrix<Type> m_r = exp(-m_d.array()); 

vector<Type> eps = m_r * v_z;
vector<Type> z = v_z;
vector<Type> v_x = m_x * v_b;
nll -= -sum(dnorm(v_z, z, sigma_z, true));
nll -= -sum(dnorm(v_y, v_x + eps, sigma, true));

return nll;
}"
write(tmb_model, file = "code/regression.cpp")

library(TMB)
compile("code/regression.cpp")
dyn.load(dynlib("code/regression"))

set.seed(122)
N <- 50
m_d <- matrix(runif(N*N), N, N)
diag(m_d) <- 0
m_r <- exp(-m_d)
z <- rnorm(N, sd = exp(0.01))

df0 <- tibble(x1 = runif(N, 1, 10),
              x2 = runif(N, 1, 10)) %>% 
  mutate(y = rnorm(N,
                   mean = 1.8 + 2.4 * x1 - 3 * x2 + m_r %*% z,
                   sd = exp(0.1)))

fm <- y ~ .
m_x <- data.matrix(model.matrix(fm, data = df0))

obj <- MakeADFun(
  data = list(m_x = m_x, v_y = df0 %>% pull(y), m_d = m_d), 
  parameters = list(v_b = c(1.8, 2.4, -3),
                    v_z = z,
                    log_sigma = 0.1,
                    log_sigma_z = 0.01),
  random = "v_z",
  DLL = "regression")
optim(par = obj$par, fn = obj$fn, method = "BFGS")

opt <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)
(rep <- sdreport(obj))
