model {
  
  # likelihood --------------------------------------------------------------
  
  for (n in 1:Nsample) {
    Y[n] ~ dbern(p[Site[n], Species[n]])
  }
  
  for (s in 1:Nsite) {
    for (i in 1:Nspecies) {
      
      logit(p[s, i]) <- beta0[group[s], i] + inprod(beta[, i], x[s, ])
      
    }
  }
  
  
  # prior -------------------------------------------------------------------
  
  tau0 <- 0.1
  phi ~ dunif(0, 1)
  
  for (i in 1:Nspecies) {
    
    for (g in 1:Ng) {
      beta0[g, i] ~ dnorm(mu_beta0, tau_beta0)
    } #g
    
    for (j in (i + 1):Nspecies) {
      ## symmetric assumption
      b[i, j] ~ dnorm(0, tau_b[i, j])
      tau_b[i, j] <- z[i, j] * tau0 + (1 - z[i, j]) * 100
      z[i, j] ~ dbern(phi)
      
      beta[i, j] <- b[i, j]
      beta[j, i] <- b[i, j]
    } #j
    
    for (k in i) {
      ## diagonal elements are zero
      beta[i, k] <- 0
    } #k
    
  } # i

  mu_beta0 ~ dnorm(0, tau0)
  tau_beta0 ~ dscaled.gamma(2.5, 3)
  
}

data {
  
  for (s in 1:Nsite) {
    group[s] <- group0[s, 1]
  }
    
  for (n in 1:Nsample) {
    group0[Site[n], Species[n]] <- Group[n]
    x[Site[n], Species[n]] <- Y[n]
  }
  
}