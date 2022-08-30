
# setup -------------------------------------------------------------------

rm(list = ls())
pacman::p_load(runjags,
               tidyverse)


# data --------------------------------------------------------------------

list_df <- readRDS("output/df_sim.rds")

df0 <- list_df[[1]] %>% 
  filter(group == 1)


# jags setup --------------------------------------------------------------

## parameters ####
para <- c("phi",
          "beta0",
          "beta")

## model file ####
m <- runjags::read.jagsfile("code/model_mrf.R")

## mcmc setup ####
n_ad <- 1000
n_iter <- 2.0E+3
n_thin <- max(3, ceiling(n_iter / 500))
n_burn <- ceiling(max(10, n_iter/2))
n_sample <- ceiling(n_iter / n_thin)
n_chain <- 4

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1


# jags --------------------------------------------------------------------

d_jags <- list(Site = df0$site,
               Species = df0$species,
               Y = df0$occupancy,
               Nsample = nrow(df0),
               Nspecies = n_distinct(df0$species),
               Nsite = n_distinct(df0$site))

## run jags ####
post <- runjags::run.jags(m$model,
                          monitor = para,
                          data = d_jags,
                          n.chains = n_chain,
                          inits = inits,
                          method = "parallel",
                          burnin = n_burn,
                          sample = n_sample,
                          adapt = n_ad,
                          thin = n_thin,
                          n.sims = n_chain,
                          module = "glm")

mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)


# validation --------------------------------------------------------------

df_est <- mcmc_summary %>% 
  as_tibble(rownames = "param") %>% 
  mutate(param_id = str_remove(param, "\\[.{1,}\\]")) %>% 
  select(param,
         param_id,
         low = `2.5%`,
         high = `97.5%`,
         median = `50%`)


df_beta0 <- df_est %>% 
  filter(param_id == "beta0") %>% 
  mutate(true = distinct(df0, alpha0) %>% pull(alpha0))

df_beta <- df_est %>% 
  filter(param_id == "beta") %>% 
  mutate(true = c(list_df[[2]]))

df_plot <- bind_rows(df_beta0, df_beta)
  
  
# visualize ---------------------------------------------------------------

ggplot(df_plot) +
  geom_abline(intercept = 0,
              slope = 1,
              color = grey(0.75)) +
  geom_segment(aes(x = true,
                   xend = true,
                   y = low,
                   yend = high),
               alpha = 0.1) +
  geom_point(aes(x = true,
                 y = median),
             alpha = 0.5) +
  facet_wrap(facets = ~param_id) +
  theme_bw()
