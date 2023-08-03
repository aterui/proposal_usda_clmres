# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

df_fish <- readRDS("data_fmt/data_mrb_fish_binary.rds")
df_site <- readRDS("data_fmt/data_mrb_site.rds")

## subsampling
df_fish <- df_fish %>% 
  left_join(df_site %>% dplyr::select(siteid, index),
            by = "siteid")

sp <- df_fish %>% 
  group_by(species) %>% 
  summarize(n = sum(presence)) %>% 
  filter(n > 200) %>% 
  pull(species)

df_fish <- df_fish %>% 
  filter(species %in% sp) %>% 
  mutate(sp_index = as.numeric(factor(species))) %>% 
  arrange(species)

# glm ---------------------------------------------------------------------

library(glmnet)
data(BinomialExample)
x <- BinomialExample$x
y <- BinomialExample$y
cvfit <- cv.glmnet(x, y, family = "binomial", type.measure = "class")
s <- cvfit$lambda.min
fit <- glmnet(x, y, family = "binomial")
coef(fit, s = s)
