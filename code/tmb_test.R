
# setup -------------------------------------------------------------------

#source("code/library.R")
library(TMB)
library(dplyr)

# model -------------------------------------------------------------------
list.files("code", pattern = c("\\.dll$|\\.o$"), full.names = T) %>% 
  file.remove()

compile("code/sem.cpp")


# simulated data ----------------------------------------------------------

set.seed(11)
N <- 40
theta <- 0.1
sigma <- 0.1
eps <- rnorm(N, sd = 0.1)
X <- cbind(1, rnorm(N), rnorm(N))
beta <- c(0.01, 0.8, 0.2)
coord <- cbind(runif(N, 0, 100), runif(N, 0, 100))

D <- coord %>% 
  dist(diag = TRUE, upper = TRUE) %>% 
  data.matrix()

S <- sigma * exp(- theta * D)
y <- mvtnorm::rmvnorm(1, mean = X %*% beta, sigma = S) %>% c()


# fitting -----------------------------------------------------------------

.valid_link <- c(identity = 0,
                 log = 1)

.valid_family <- c(gaussian = 0,
                   poisson = 1)

dyn.load(dynlib("code/sem"))

obj <- MakeADFun(data = list(y = y,
                             D = D,
                             X = X,
                             link = .valid_link[poisson()$link],
                             family = .valid_family[poisson()$family]), 
                 parameters = list(b = beta,
                                   u = rep(1, length(y)),
                                   log_sigma0 = log(0.1),
                                   log_sigma = log(0.1),
                                   log_theta = log(0.1)),
                 random = "u",
                 DLL = "sem")

opt <- nlminb(start = obj$par,
              obj = obj$fn,
              gr = obj$gr, control = list(iter.max = 300,
                                          eval.max = 400))

rep <- sdreport(obj)

summary(rep, "report", p.value = TRUE) %>% 
  round(3)


# tmb ---------------------------------------------------------------------

library(glmmTMB)
pos <- numFactor(coord[,1], coord[,2])

df0 <- tibble(y = y, pos = pos, group = factor(rep(1, length(y))))
fit <- glmmTMB(y ~ 1 + X[, 2] + X[,3] + exp(pos + 0 | group),
               data = df0, family = "gaussian",
               dispformula=~0)

summary(fit)
