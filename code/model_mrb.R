model {
  
  # likelihood --------------------------------------------------------------
  
  for (n in 1:Nsample) {
    Y[n] ~ dbern(p[Site[n], Species[n]])
  }
  
  for (s in 1:Nsite) {
    for (i in 1:Nspecies) {
      logit(p[s, i]) <- inprod(alpha[, i], X[s, ]) + inprod(beta[, i], Q[s, ])
    }
  }
  
  
  # prior -------------------------------------------------------------------
  
  tau0 <- 0.1
  phi ~ dunif(0, 1)
  
  for (i in 1:Nspecies) {
    
    for (k in 1:Nk) {
      alpha[k, i] ~ dnorm(0, tau0)
    } #k
    
    for (j in (i + 1):Nspecies) {
      ## sparse priors
      b[i, j] ~ dnorm(0, tau_b[i, j])
      tau_b[i, j] <- z[i, j] * tau0 + (1 - z[i, j]) * 100
      z[i, j] ~ dbern(phi)
      
      ## symmetric assumption
      beta[i, j] <- b[i, j]
      beta[j, i] <- b[i, j]
    } #j
    
    for (k in i) {
      ## diagonal elements are zero
      beta[i, k] <- 0
    } #k
    
  } # i
  
}

data {
  
  for (n in 1:Nsample) {
    # site X species matrix
    Q[Site[n], Species[n]] <- Y[n]
  }
  
}