# setup -------------------------------------------------------------------

rm(list = ls())
source(here::here("code/library.R"))


# data --------------------------------------------------------------------

df0 <- readRDS("data_fmt/data_tp_sub.rds")
df_dist <- readRDS("data_fmt/data_distance.rds") %>% 
  mutate(across(.cols = c("distance", "up", "down"),
                \(x) as.numeric(x)),
         fcon = ifelse(up > 0, 0, 1))

m <- df_dist %>% 
  pivot_wider(id_cols = from,
              names_from = to,
              values_from = distance) %>% 
  select(c(1, (order(.$from) + 1))) %>% 
  arrange(from) %>% 
  dplyr::select(-1) %>% 
  data.matrix()


# jags setup --------------------------------------------------------------

## parameters ####
para <- c("b0",
          "sigma",
          "theta")

## model file ####
m <- runjags::read.jagsfile("code/model_geostr.R")

## mcmc setup ####
n_ad <- 1000
n_iter <- 1.0E+3
n_thin <- max(3, ceiling(n_iter / 250))
n_burn <- ceiling(max(10, n_iter/2))
n_sample <- ceiling(n_iter / n_thin)
n_chain <- 4

inits <- replicate(n_chain,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1


# jags --------------------------------------------------------------------

d_jags <- list(Y = df0$value,
               From = df_dist$from,
               To = df_dist$to,
               Down = df_dist$down, # downstream distance
               Distance = df_dist$distance, # network distance
               W = df_dist$fcon, # flow connectance
               Nsample = nrow(df0),
               Ncomb = nrow(df_dist))

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

(mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc))

# while(max(mcmc_summary$Rhat, na.rm = T) >= 1.1) {
#   post <- runjags::extend.jags(post,
#                                burnin = 0,
#                                sample = n_sample,
#                                adapt = n_ad,
#                                thin = n_thin,
#                                n.sims = n_chain,
#                                combine = TRUE)
#   
#   mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)
# }
# 
# mcmc_sample <- post$mcmc
# n_total_mcmc <- (post$sample / n_sample) * n_iter + n_burn
