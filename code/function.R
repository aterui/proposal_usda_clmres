

# arc2d8 ------------------------------------------------------------------

arc2d8 <- function(x) {
  ## whitebox uses flow pointer D8
  # 64, 128,  1
  # 32,   0,  2 
  # 16,   8,  4
  
  ## ArcGIS uses D8 algorithm
  # 32, 64, 128
  # 16,  0,   1 
  #  8,  4,   2
  
  # check class()
  if(class(x) != "SpatRaster") stop("Provide data in class 'SpatRaster'")
  
  # begin with northeast through north
  fdir_arc <- c(0, 2^(0:7), 247, 255)
  fdir_d8 <- c(0, fdir_arc[3:9], fdir_arc[2], NA, NA)
  y <- terra::subst(x, from = fdir_arc, to = fdir_d8)
  
  return(y)  
}

# from matrix to vector ---------------------------------------------------

m2v <- function(x) {
  
  if (is.null(colnames(x))) {
    index_from <- rep(1:dim(x)[1], times = dim(x)[2])
    index_to <- rep(1:dim(x)[2], each = dim(x)[1])
  } else {
    index_from <- rep(colnames(x), times = dim(x)[2])
    index_to <- rep(rownames(x), each = dim(x)[1])
  }
  
  return(dplyr::tibble(value = c(x),
                       from = index_from,
                       to = index_to))
  
}


# wrapper jags ------------------------------------------------------------

geojags <- function(formula,
                    data,
                    dm,
                    init,
                    n_chain = 3,
                    adapt = 1000,
                    sample = 1000,
                    burnin = floor(sample / 2),
                    thin = 3) {
  
  ## parameters ####
  para <- c("b",
            "sigma",
            "theta",
            "lambda")
  
  ## model file ####
  model_text <- "
  model {
    ninfo <- 0.1
    df0 <- 6
    
    # prior -------------------------------------------------------------------
    
    for (k in 1:Nb) {
      b[k] ~ dnorm(0, ninfo)
    }
    
    for (k in 1) {
      tau[k] ~ dscaled.gamma(2.5, df0)  
      sigma[k] <- pow(tau[k], -0.5)
    }
    
    theta ~ dunif(0, 10)
    lambda ~ dnorm(0, ninfo)T(0, 1)
    
    # likelihood --------------------------------------------------------------
    
    for(i in 1:Nsample) {
      Y[i] ~ dnorm(mu[i], tau[1])
      mu[i] <- inprod(X[i, ], b[]) + lambda * inprod(Q[i, ], y[])
    }
    
    for (i in 1:Nsample) {
      S[i, i] <- 0
      Q[i, 1:Nsample] <- S[i, ] / (sum(S[i, ]) + step(z[i]))
      z[i] <- -sum(S[i, ])
      for (j in (i + 1):Nsample) {
        S[i, j] <- exp(-theta * m_dist[i, j])# * m_w[i, j]
        S[j, i] <- S[i, j]
      }
    }
    
  }
  
  data {
    
    y <- Y
    
    for (n in 1:Ncomb) {
      m_dist[From[n], To[n]] <- Distance[n]
      m_w[From[n], To[n]] <- W[n]
    }
    
  }"
  
  model_path <- paste(tempdir(), "model.R", sep = "\\")
  write(model_text, model_path)
  m <- runjags::read.jagsfile(model_path)
  
  ## data ####
  
  if(missing(init)) {
    inits <- replicate(n_chain,
                       list(.RNG.name = "base::Mersenne-Twister",
                            .RNG.seed = NA),
                       simplify = FALSE)
    
    for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1
  }
  
  Y <- model.frame(formula, data)[,1]
  X <- model.matrix(formula, data)
  
  
  # jags --------------------------------------------------------------------
  d_jags <- list(
    # regression
    Y = Y,
    X = X,
    Nsample = length(Y),
    Nb = ncol(X),
    # spatial component
    From = dm$from,
    To = dm$to,
    Distance = dm$distance, # network distance
    W = dm$c_flow, # flow connectance
    Ncomb = nrow(dm))
  
  ## run jags ####
  post <- runjags::run.jags(m$model,
                            monitor = para,
                            data = d_jags,
                            n.chains = n_chain,
                            inits = inits,
                            method = "parallel",
                            burnin = burnin,
                            sample = sample,
                            adapt = adapt,
                            thin = thin,
                            n.sims = n_chain,
                            module = "glm")
  
  return(MCMCvis::MCMCsummary(post$mcmc))
}

