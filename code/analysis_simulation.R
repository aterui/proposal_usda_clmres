
# setup -------------------------------------------------------------------

rm(list = ls())

pacman::p_load(tidyverse,
               foreach)

source(here::here("code/function_simulator.R"))


# simulated data ----------------------------------------------------------

n_step <- 100
n_species <- 20
n_site <- 10
n_group <- 20
phi <- 0.5

m_alpha <- rnorm(n_species * n_species,
                 mean = 0,
                 sd = 1) * rbinom(n_species * n_species, 1, phi) %>% 
  matrix(ncol = n_species,
         nrow = n_species)

m_alpha[upper.tri(m_alpha)] <- 0
m_alpha <- m_alpha + t(m_alpha)
diag(m_alpha) <- 0

df_sim <- foreach(i = seq_len(n_group),
                  .combine = bind_rows) %do% {
                    list_gibbs <- f_gibbs(n_species = n_species,
                                          n_site = n_site,
                                          n_step = n_step,
                                          n_burn = n_step - 1,
                                          alpha = m_alpha)
                    
                    df0 <- list_gibbs$Y %>% 
                      as_tibble() %>% 
                      mutate(group = i,
                             alpha0 = rep(list_gibbs$alpha0,
                                          each = n_site))
                    
                    return(df0)
                  }

# export
saveRDS(list(df_sim, m_alpha), file = "output/df_sim.rds")



