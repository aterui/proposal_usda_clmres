
# setup -------------------------------------------------------------------

source("code/library.R")
library(TMB)

# model -------------------------------------------------------------------
list.files("code", pattern = c("\\.dll$|\\.cpp$|\\.o$"), full.names = T) %>% 
  file.remove()

tmb_model <- "
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
}"
write(tmb_model, file = "code/sem.cpp")
compile("code/sem.cpp")


# simulated data ----------------------------------------------------------

set.seed(121)
N <- 200
theta <- 10
eps <- rnorm(N, sd = 0.1)
X <- cbind(1, rnorm(N), rnorm(N))
beta <- c(0.01, 0.8, 0.2)

D <- cbind(runif(N, 0, 1000), runif(N, 0, 1000)) %>% 
  dist(diag = TRUE, upper = TRUE) %>% 
  data.matrix()
S <- 0.1 * exp(-D / theta)

y_hat <- X %*% beta
y <- mvtnorm::rmvnorm(1, mean = y_hat, sigma = S) %>% c()

# fitting -----------------------------------------------------------------

dyn.load(dynlib("code/sem"))

obj <- MakeADFun(data = list(y = y,
                             D = D,
                             X = X), 
                 parameters = list(b = rep(0, ncol(X)),
                                   log_sigma = log(0.1),
                                   log_theta = log(10)),
                 DLL = "sem")

opt <- nlminb(start = obj$par,
              obj = obj$fn,
              gr = obj$gr)

rep <- sdreport(obj)

summary(rep, "report", p.value = TRUE) %>% 
  round(3)
