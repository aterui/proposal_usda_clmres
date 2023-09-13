# setup -------------------------------------------------------------------

rm(list = ls()); gc()
source("code/library.R")
source("code/function.R")

set.seed(122)
df_site <- readRDS("data_fmt/data_mrb_site.rds") %>% 
  #sample_n(size = 100) %>% 
  mutate(x_log_area = c(scale(log(area))),
         x_logit_agri = c(scale(boot::logit(frac_agri))),
         x_btw = c(scale(btw)))

df_fish <- readRDS("data_fmt/data_mrb_fish_binary.rds") %>% 
  filter(siteid %in% unique(df_site$siteid))

## subsampling
df_fish <- df_fish %>% 
  left_join(df_site %>% dplyr::select(siteid, index),
            by = "siteid")

sp <- df_fish %>% 
  group_by(species) %>% 
  summarize(n = sum(presence)) %>% 
  filter(n > floor(nrow(df_site) * 0.25)) %>% 
  pull(species)

Y <- df_fish %>% 
  filter(species %in% sp) %>% 
  arrange(species) %>% 
  pivot_wider(id_cols = c(siteid, index),
              names_from = species,
              values_from = presence)

X <- df_site %>% 
  dplyr::select(siteid, index, starts_with("x_"))

df_m <- Y %>% 
  left_join(X,
            by = c("siteid", "index")) %>% 
  dplyr::select(-c(siteid, index))

# glm ---------------------------------------------------------------------

library(glmnet)

nx <- 1 + sum(str_detect(colnames(df_m), "x_"))
ns <- ncol(df_m) - (nx - 1)
A <- matrix(NA, ns, ns)
B <- matrix(NA, nx, ns)

for (i in 1:length(sp)) {
  print(i)  
  y <- df_m[, i] %>% pull()
  X <- df_m[, -i] %>% data.matrix()
  
  cv <- cv.glmnet(X, y, family = "binomial", type.measure = "class")
  m <- glmnet(x = X, y = y,
              family = "binomial", 
              lambda = cv$lambda.min)
  
  beta <- coef(m)
  b <- as.vector(beta)
  id <- rownames(beta)
  
  spindex <- !str_detect(id, "(Intercept)|x_.{1,}")
  xindex <- str_detect(id, "(Intercept)|x_.{1,}")
  
  A[i, -i] <- b[spindex]
  B[,i] <- b[xindex]
  
}

Beta <- as_tibble(t(B)) %>% 
  mutate(commonname = colnames(Y[,-c(1,2)]))

# fish trait --------------------------------------------------------------

df_trait <- read_csv(here::here("data_raw/fish_trait.csv")) %>% 
  rename_with(.fn = str_to_lower) %>% 
  mutate(commonname = str_replace_all(commonname,
                                      "\\s",
                                      "_"),
         commonname = str_remove_all(commonname,
                                     "-"),
         commonname = str_to_lower(commonname),
         commonname = ifelse(commonname == "eastern_blacknose_dace",
                             "blacknose_dace",
                             commonname)) %>% 
  filter(commonname %in% sp)

df_plot <- Beta %>% 
  left_join(df_trait, by = "commonname") %>% 
  mutate(primary_consumer = ifelse(detritus + algphyto + macvascu > 0, 1, 0),
         predator = ifelse(invlvfsh + fshcrcrb, 1, 0),
         substrate = ifelse(claysilt + sand > 1, 1, 0),
         lithophils = ifelse(a_1_3a + a_2_3a + a_1_3b + a_2_3b +
                               b_1_3a + b_2_3a + b_2_3b > 0,
                             1, 0)) %>% 
  drop_na(lithophils)


# figure ------------------------------------------------------------------

fm <- lm(log(fecundity) ~ log(maxtl),
         df_plot)
df_plot <- df_plot %>% 
  mutate(fec = resid(fm))

MASS::rlm(V3 ~ fec + log(maxtl) + lithophils + substrate + primary_consumer,
          df_plot) %>% 
  summary()

df_plot %>% 
  ggplot(aes(y = V3,
             x = log(maxtl))) +
  geom_point()
