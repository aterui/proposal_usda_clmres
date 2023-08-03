# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

df_fish <- readRDS("data_fmt/data_mrb_fish_binary.rds")
df_site <- readRDS("data_fmt/data_mrb_site.rds")

## subsampling
df_fish <- df_fish %>% 
  left_join(df_site %>% dplyr::select(siteid, index),
            by = "siteid") %>% 
  filter(index < 51)
  
sp <- df_fish %>% 
  group_by(species) %>% 
  summarize(n = sum(presence)) %>% 
  filter(n > 20) %>% 
  pull(species)

df_fish <- df_fish %>% 
  filter(species %in% sp) %>% 
  mutate(sp_index = as.numeric(factor(species))) %>% 
  arrange(species)

df_site <- df_site %>% 
  filter(index < 51)

# jags setup --------------------------------------------------------------

## parameters ####
para <- c("phi",
          "alpha",
          "beta")

inits <- replicate(3,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (j in 1:3) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1

## data ####
X <- model.matrix(~ log(area) + boot::logit(frac_agri) + btw, df_site)

d_jags <- with(df_fish,
               list(Y = presence,
                    X = X,
                    Site = index,
                    Species = sp_index,
                    Nsite  = n_distinct(siteid),
                    Nspecies = n_distinct(species),
                    Nsample = nrow(df_fish),
                    Nk = ncol(X))
               )

m <- runjags::read.jagsfile("code/model_mrb")

## mcmc ####
sample <- 500
adapt <- 100
thin <- 3
burnin <- floor(sample * thin / 2)

# run ---------------------------------------------------------------------

post <- runjags::run.jags(m$model,
                          monitor = para,
                          data = d_jags,
                          n.chains = 3,
                          inits = inits,
                          method = "parallel",
                          burnin = burnin,
                          sample = sample,
                          adapt = adapt,
                          thin = thin,
                          n.sims = 3,
                          module = "glm")

MCMCvis::MCMCsummary(post$mcmc)
