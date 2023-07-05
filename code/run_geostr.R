# setup -------------------------------------------------------------------

rm(list = ls())
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


# jags setup --------------------------------------------------------------

(fit <- geojags(log(value) ~ boot::logit(frac_agri) + log(area),
                data = df0,
                dm = df_dist,
                adapt = 100,
                sample = 1000,
                thin = 4))