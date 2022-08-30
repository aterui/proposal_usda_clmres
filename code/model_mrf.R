model {
  
  # likelihood --------------------------------------------------------------
  for (n in 1:Nsample) {
    Y[n] ~ dbern(p[Site[n], Species[n]])
  }
  
  for (s in 1:Nsite) {
    for (i in 1:Nspecies) {
      
      logit(p[s, i]) <- beta0[i] + inprod(beta[, i], x[s, ])
      
    }
  }
  
  
  # prior -------------------------------------------------------------------
  tau0 <- 0.1
  phi ~ dunif(0, 1)
  
  for (i in 1:Nspecies) {
    
    beta0[i] ~ dnorm(0, tau0)
    
    for (j in (i + 1):Nspecies) {
      ## symmetric assumption
      b[i, j] ~ dnorm(0, tau_b[i, j])
      tau_b[i, j] <- z[i, j] * 0.01 + (1 - z[i, j]) * 100
      z[i, j] ~ dbern(phi)
      
      beta[i, j] <- b[i, j]
      beta[j, i] <- b[i, j]
    }
    
    for (k in i) {
      ## diagonal elements are zero
      beta[i, k] <- 0
    }
    
  }
  
}

data {
  
  for (n in 1:Nsample) {
    x[Site[n], Species[n]] <- Y[n]
  }
  
}