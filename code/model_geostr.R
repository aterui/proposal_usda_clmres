model {
  
  ninfo <- 0.1
  df0 <- 3
  
  # prior -------------------------------------------------------------------
  
  b0 ~ dnorm(0, ninfo)
  
  theta[1] ~ dunif(0, theta[2])
  theta[2] ~ dunif(0, 10)
  
  for (k in 1:3) {
    tau[k] ~ dscaled.gamma(2.5, df0)  
    sigma[k] <- pow(tau[k], -0.5)
  }
  
  # likelihood --------------------------------------------------------------
  
  for(i in 1:Nsample) {
    log_y[i] ~ dnorm(mu[i], tau[1])
    mu[i] <- b0 + 
      inprod(rho_down[i, ], z_down[]) + 
      inprod(rho_dist[i, ], z_dist[])
    
    z_down[i] ~ dnorm(0, tau[2])
    z_dist[i] ~ dnorm(0, tau[3])
  }
  
  for (i in 1:Nsample) {
    for (j in 1:Nsample) {
      rho_down[i, j] <- exp(-theta[1] * m_down[i, j]) * m_w[i, j]
      rho_dist[i, j] <- exp(-theta[2] * m_dist[i, j])
    }
  }
  
}

data {
  
  for (i in 1:Nsample) {
    log_y[i] <- log(Y[i])
  }
  
  for (n in 1:Ncomb) {
    m_down[From[n], To[n]] <- Down[n]
    m_dist[From[n], To[n]] <- Distance[n]
    m_w[From[n], To[n]] <- W[n]
  }
  
}