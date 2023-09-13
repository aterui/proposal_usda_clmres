
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))
source(here::here("code/function.R"))
beta <- readRDS("output/data_flow_beta.rds") %>% 
  as_tibble(rownames = "date") %>% 
  mutate(date = as.Date(date)) %>% 
  rename(const = `(Intercept)`,
         z = `log(area)`,
         b = `frac_agri`,
         zb = `log(area):frac_agri`)

# data --------------------------------------------------------------------

set.seed(123)
df_obs <- readRDS("data_fmt/data_tpsub_master.rds") %>% 
  drop_na(site_id) %>% 
  left_join(beta, by = "date") %>% 
  mutate(log_flow = const + z * log(area) + b * frac_agri + zb * log(area) * frac_agri,
         load = value * exp(log_flow))

df0 <- df_obs %>% 
  #sample_n(50) %>% 
  arrange(site_id)
  
v_omega <- df0 %>% pull(area)

m_weight <- outer(v_omega, v_omega, FUN = "/") %>% 
  sqrt()

colnames(m_weight) <- rownames(m_weight) <- df0 %>% pull(site_id)
df_w <- m2v(m_weight) %>% 
  mutate(across(.cols = c(from, to),
                \(x) as.numeric(x)))

df_dist <- readRDS("data_fmt/data_distance.rds") %>% 
  filter(type == "observation") %>% 
  mutate(across(.cols = c("distance", "up", "down"),
                \(x) as.numeric(x)),
         from = str_remove(from, "o_") %>% as.numeric(),
         to = str_remove(to, "o_") %>% as.numeric()) %>% 
  right_join(df_w,
             by = c("from", "to")) %>% 
  mutate(c_flow = ifelse(up > 0, 0, value)) %>% 
  filter(from <= to)

## distance matrix
D <- df_dist %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = from,
              values_from = distance,
              names_from = to, values_fill = 0) %>% 
  arrange(from) %>% 
  select(-from) %>% 
  data.matrix()

D <- D + t(D)

## weight matrix
W <- df_dist %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = from,
              values_from = c_flow,
              names_from = to,
              values_fill = 0) %>% 
  arrange(from) %>% 
  select(-from) %>% 
  data.matrix()

W <- W + t(W)
diag(W) <- 1


# fit ---------------------------------------------------------------------
library(snglmm)
m0 <- lm(log(load) ~ frac_agri + frac_wet + log_flow,
         data = df0)

m <- snglmm(log(load) ~ frac_agri + frac_wet + log_flow,
            data = df0,
            spatial = "exp",
            control = list(eval.max = 2000,
                           iter.max = 1000),
            inits = list(log_sigma = log(0.5),
                         log_phi = log(1),
                         log_lambda = log(0.5)),
            D = D,
            W = W)


# # prediction spatial ------------------------------------------------------
# 
# df_pred <- df_obs %>% 
#   filter(!(site_id %in% pull(df0, site_id))) %>% 
#   arrange(site_id)
# 
# cD <- readRDS("data_fmt/data_distance.rds") %>% 
#   filter(type == "observation") %>% 
#   mutate(across(.cols = c("distance", "up", "down"),
#                 \(x) as.numeric(x)),
#          from = str_remove(from, "o_") %>% as.numeric(),
#          to = str_remove(to, "o_") %>% as.numeric()) %>% 
#   filter(from %in% df0$site_id,
#          to %in% df_pred$site_id) %>% 
#   arrange(to) %>% 
#   pivot_wider(id_cols = from,
#               values_from = distance,
#               names_from = to, values_fill = 0) %>% 
#   arrange(from) %>% 
#   select(-from) %>% 
#   data.matrix()
# 
# C <- readRDS("data_fmt/data_distance.rds") %>% 
#   filter(type == "observation") %>% 
#   mutate(across(.cols = c("distance", "up", "down"),
#                 \(x) as.numeric(x)),
#          from = str_remove(from, "o_") %>% as.numeric(),
#          to = str_remove(to, "o_") %>% as.numeric()) %>% 
#   filter(from %in% df0$site_id,
#          to %in% df_pred$site_id) %>% 
#   mutate(across(.cols = c("distance", "up", "down"),
#                 \(x) as.numeric(x)),
#          from = str_remove(from, "o_") %>% as.numeric(),
#          to = str_remove(to, "o_") %>% as.numeric(),
#          c_flow = ifelse(up > 0, 0, 1)) %>% 
#   arrange(to) %>% 
#   pivot_wider(id_cols = from,
#               values_from = c_flow,
#               names_from = to, values_fill = 0) %>% 
#   arrange(from) %>% 
#   select(-from) %>% 
#   data.matrix()
# 
# v_omega <- readRDS("data_fmt/data_tpsub_master.rds") %>%
#   filter(site_id %in% df0$site_id) %>%
#   arrange(site_id) %>% 
#   pull(area)
# 
# v_omega_p <- readRDS("data_fmt/data_tpsub_master.rds") %>%
#   filter(!(site_id %in% df0$site_id)) %>%
#   arrange(site_id) %>% 
#   pull(area)
# 
# cW0 <- outer(v_omega, v_omega_p, FUN = "/")
# cW1 <- apply(cW0, 2, FUN = function(x) ifelse(x > 1, 1/x, x))
# cW <- cW1 * C
# 
# X0 <- model.matrix(~ frac_agri + frac_wet + log_flow, df_pred)
# 
# y_hat <- kriging(m,
#                  newdata = X0,
#                  cD = cD,
#                  cW = C) %>% exp()
# 
# df_out <- tibble(y_hat, site_id = df_pred$site_id)
# 
# df_out <- df_out %>% 
#   left_join(df_pred, by = "site_id")
# 
# 
# # prediction non spatial --------------------------------------------------
# 
# df_pred <- df_obs %>% 
#   filter(!(site_id %in% pull(df0, site_id))) %>% 
#   arrange(site_id)
# 
# X0 <- model.matrix(~ frac_agri + frac_wet + log_flow, df_pred)
# y_hat0 <- exp(drop(X0 %*% coef(m0)))
# 
# df_pred <- df_out %>% 
#   mutate(y_hat0 = y_hat0) %>% 
#   pivot_longer(cols = c(y_hat, y_hat0),
#                names_to = "type",
#                values_to = "y") %>% 
#   mutate(type = ifelse(type == "y_hat", "spatial", "non_spatial"))
# 
# 
# # plot --------------------------------------------------------------------

fit0 <- predict(m0)
fit <- attr(m$fit, "report")$eta + attr(m$fit, "report")$u

df_fit <- tibble(type = rep(c("Non-spatial", "Spatial"),
                            each = length(fit)),
                 flow = rep(exp(df0$log_flow), 2),
                 cons_y = exp(c(fit0, fit)) / flow,
                 cons_y_obs = rep(df0$value, 2))

#options(scipen = 100, digits = 1)

g_fit <- ggplot(df_fit,
                aes(y = cons_y,
                    x = cons_y_obs,
                    color = type)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed") +
  #geom_smooth(method = "gam", se = FALSE) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(y = "Predicted TP (mg / L)",
       x = "Observed TP (mg / L)",
       color = "Model") +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text =  element_text(size = 20))

ggsave(g_fit,
       filename = "output/fig_tpfit.pdf",
       width = 8, height = 5)
