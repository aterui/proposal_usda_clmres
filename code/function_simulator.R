pacman::p_load(tidyverse,
               foreach)

f_gibbs <- function(n_step = 100,
                    n_site = 100,
                    n_species = 20,
                    n_burn = 0,
                    alpha = NULL,
                    alpha_mu = 0,
                    alpha_sd = 1) {
  
  ## row id for final output
  row_id <- seq(1, n_site * n_species * n_step + 1, by = n_site * n_species)
  
  ## alpha
  if (is.null(alpha)) {
    m_alpha <- rnorm(n_species * n_species,
                     mean = alpha_mu,
                     sd = alpha_sd) %>% 
      matrix(ncol = n_species,
             nrow = n_species)
    
    m_alpha[upper.tri(m_alpha)] <- 0
    m_alpha <- m_alpha + t(m_alpha)
    diag(m_alpha) <- 0
  } else {
    m_alpha <- alpha
  }
  
  v_alpha0 <- boot::logit(runif(n_species, 0.1, 0.8))
  
  m_y <- rbinom(n_site * n_species,
                1,
                boot::inv.logit(rep(v_alpha0, each = n_site))) %>% 
    matrix(nrow = n_site, ncol = n_species)
  
  
  ## Y frame for output with site & species ID  
  Y <- foreach(i = 1:n_step,
               .combine = bind_rows) %do% {
                 Y <- tibble(site = rep(1:n_site, n_species),
                             species = rep(1:n_species, each = n_site),
                             step = NA,
                             occupancy = NA)
               } %>% 
    data.matrix()
  
  m_p <- matrix(NA, nrow = n_site, ncol = n_species)
  
  for (i in 1:n_step) {
    for (j in 1:n_species) {
      m_p[, j] <- boot::inv.logit(v_alpha0[j] + m_y %*% m_alpha[, j])
      m_y[, j] <- rbinom(n_site, 1, m_p[, j])
    }
    Y[row_id[i]:(row_id[i + 1] - 1), "occupancy"] <- c(m_y)  
    Y[row_id[i]:(row_id[i + 1] - 1), "step"] <- i
  }
  
  if (n_burn > 0) {
    return(list(Y = Y[row_id[n_burn + 1]:(row_id[n_step + 1] - 1),],
                alpha = m_alpha,
                alpha0 = v_alpha0))
  } else {
    return(list(Y = Y,
                alpha = m_alpha,
                alpha0 = v_alpha0))
  }
  
}
