
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
DATA_SCALAR(lambda);

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
}"
write(tmb_model, file = "code/regression.cpp")
compile("code/regression.cpp")


# simulated data ----------------------------------------------------------

set.seed(123)
N <- 300
theta <- 0.1
lambda <- 1
eps <- rnorm(N, sd = 0.1)
X <- cbind(1, rnorm(N), rnorm(N))
beta <- c(0.01, 1.8, 2.4)

D <- cbind(runif(N, 0, 1), runif(N, 0, 1)) %>% 
  dist(diag = TRUE, upper = TRUE) %>% 
  data.matrix()
W <- exp(-theta * D)
diag(W) <- 0
#Q <- W / rowSums(W)

y <- solve(diag(N) - lambda * W) %*% (X %*% beta + eps) %>% 
  c()


# fitting -----------------------------------------------------------------

dyn.load(dynlib("code/regression"))

obj <- MakeADFun(data = list(y = y,
                             D = D,
                             X = X), 
                 parameters = list(b = rep(0, ncol(X)),
                                   log_sigma = log(1),
                                   log_theta = log(0.1)),
                 DLL = "regression")

opt <- nlminb(start = obj$par,
              obj = obj$fn,
              gr = obj$gr)

rep <- sdreport(obj)

summary(rep, "report", p.value = TRUE)
