
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
m0 <- lm(log(load + 1) ~ frac_agri + frac_wet * log(area),
         data = df0)

m <- snglmm(log(load + 1) ~ frac_agri + log_flow + frac_wet,
            data = df0,
            spatial = "exp",
            control = list(eval.max = 2000,
                           iter.max = 1000),
            inits = list(log_lambda = log(1)),
            D = D,
            W = W)


# prediction spatial ------------------------------------------------------

df_pred <- df_obs %>% 
  filter(!(site_id %in% pull(df0, site_id))) %>% 
  arrange(site_id)

cD <- readRDS("data_fmt/data_distance.rds") %>% 
  filter(from %in% df0$index,
         to %in% df_pred$index,
         str_detect(from, "o_")) %>% 
  mutate(across(.cols = c("distance", "up", "down"),
                \(x) as.numeric(x)),
         from = str_remove(from, "o_") %>% as.numeric(),
         to = str_remove(to, "o_") %>% as.numeric()) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = from,
              values_from = distance,
              names_from = to, values_fill = 0) %>% 
  arrange(from) %>% 
  select(-from) %>% 
  data.matrix()

C <- readRDS("data_fmt/data_distance.rds") %>% 
  filter(from %in% df0$index,
         to %in% df_pred$index,
         str_detect(from, "o_")) %>% 
  mutate(across(.cols = c("distance", "up", "down"),
                \(x) as.numeric(x)),
         from = str_remove(from, "o_") %>% as.numeric(),
         to = str_remove(to, "o_") %>% as.numeric(),
         c_flow = ifelse(up > 0, 0, 1)) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = from,
              values_from = c_flow,
              names_from = to, values_fill = 0) %>% 
  arrange(from) %>% 
  select(-from) %>% 
  data.matrix()

v_omega <- readRDS("data_fmt/data_tpsub_master.rds") %>%
  filter(!is.na(site_id),
         index %in% df0$index) %>%
  arrange(site_id) %>% 
  pull(area)

v_omega_p <- readRDS("data_fmt/data_tpsub_master.rds") %>%
  filter(!is.na(site_id),
         index != "p_35",
         !(index %in% df0$index)) %>%
  arrange(site_id) %>% 
  pull(area)

cW0 <- outer(v_omega, v_omega_p, FUN = "/")
cW1 <- apply(cW0, 2, FUN = function(x) ifelse(x > 1, 1/x, x))
cW <- cW1 * C

X0 <- model.matrix(~ frac_agri + frac_wet + log(area), df_pred)

y_hat <- kriging(m,
                 newdata = X0,
                 cD = cD,
                 cW = cW) %>% exp()
names(y_hat) <- df_pred$index
df_out <- tibble(y_hat, index = names(y_hat))

df_out <- df_out %>% 
  left_join(df_pred, by = "index")


# prediction non spatial --------------------------------------------------

df_pred <- df_obs %>% 
  filter(!(site_id %in% pull(df0, site_id))) %>% 
  arrange(site_id)

X0 <- model.matrix(~ frac_agri + frac_wet + log(area), df_pred)
y_hat0 <- exp(drop(X0 %*% coef(m0)))

df_plot <- df_out %>% 
  mutate(y_hat0 = y_hat0) %>% 
  pivot_longer(cols = c(y_hat, y_hat0),
               names_to = "type",
               values_to = "y") %>% 
  mutate(type = ifelse(type == "y_hat", "spatial", "non_spatial"))


# plot --------------------------------------------------------------------

ggplot(df_plot,
       aes(y = y,
           x = value,
           color = type)) +
  geom_point() +
  geom_abline(intercept = 0,
              slope = 1, linetype = "dashed") +
  geom_smooth(method = "lm")

