
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))
source(here::here("code/function.R"))


# data --------------------------------------------------------------------

df0 <- readRDS("data_fmt/data_tpsub_master.rds") %>% 
  drop_na(site_id) %>% 
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
  left_join(df_w,
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
m0 <- snglmm(log(value) ~ boot::logit(frac_agri) + log(area),
             data = df0)

m <- snglmm(log(value) ~ boot::logit(frac_agri) + log(area),
            data = df0,
            spatial = "exp",
            control = list(eval.max = 2000,
                           iter.max = 1000),
            D = D,
            W = W)


# prediction --------------------------------------------------------------

df_pred <- readRDS("data_fmt/data_tpsub_master.rds") %>% 
  filter(is.na(site_id),
         index != "p_35")

cD <- readRDS("data_fmt/data_distance.rds") %>% 
  filter(to %in% df_pred$index,
         str_detect(from, "o_")) %>% 
  mutate(across(.cols = c("distance", "up", "down"),
                \(x) as.numeric(x))) %>% 
  arrange(as.numeric(str_remove(to, "p_"))) %>% 
  pivot_wider(id_cols = from,
              values_from = distance,
              names_from = to, values_fill = 0) %>% 
  arrange(as.numeric(str_remove(from, "o_"))) %>% 
  select(-from) %>% 
  data.matrix()

C <- readRDS("data_fmt/data_distance.rds") %>% 
  filter(to %in% df_pred$index,
         str_detect(from, "o_")) %>% 
  mutate(across(.cols = c("distance", "up", "down"),
                \(x) as.numeric(x)),
         c_flow = ifelse(up > 0, 0, 1)) %>% 
  arrange(as.numeric(str_remove(to, "p_"))) %>% 
  pivot_wider(id_cols = from,
              values_from = c_flow,
              names_from = to, values_fill = 0) %>% 
  arrange(as.numeric(str_remove(from, "o_"))) %>% 
  select(-from) %>% 
  data.matrix()

v_omega <- readRDS("data_fmt/data_tpsub_master.rds") %>%
  filter(!is.na(site_id)) %>%
  pull(area)

v_omega_p <- readRDS("data_fmt/data_tpsub_master.rds") %>%
  filter(is.na(site_id),
         index != "p_35") %>%
  pull(area)

cW0 <- outer(v_omega, v_omega_p, FUN = "/")
cW1 <- apply(cW0, 2, FUN = function(x) ifelse(x > 1, 1/x, x))
cW <- cW1 * C

X0 <- model.matrix(~ boot::logit(frac_agri) + log(area), df_pred)

y_hat <- kriging(m, newdata = X0, cD = cD, cW = cW) %>% exp()
names(y_hat) <- df_pred$index
df_out <- tibble(y_hat, index = names(y_hat))

df_y_hat <- df_pred %>% 
  left_join(df_out, by = "index") %>% 
  dplyr::select(y_hat, node_id)


# map ---------------------------------------------------------------------

sf_str <- st_read(here::here("data_fmt/epsg4326_str_mrb_raw.gpkg"))

sf_str <- sf_str %>% 
  left_join(df_y_hat, by = c("from" = "node_id"))

g_ex <- ggplot(sf_str) +
  geom_sf(aes(color = y_hat),
          linewidth = 0.8) +
  #scale_color_viridis_c() +
  MetBrewer::scale_color_met_c("Hiroshige", direction = -1) +
  geom_sf(data = st_as_sf(df0),
          aes(color = value),
          alpha = 0.5, shape = 21) +
  theme_bw() +
  labs(color = "TP (mg/L)")

ggsave(g_ex,
       filename = "output/figure_tp_ex.pdf", 
       height = 5, width = 6)
