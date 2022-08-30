list_df <- readRDS("output/df_sim.rds")

m <- list_df[[1]] %>% 
  filter(group == 1) %>% 
  pivot_wider(id_cols = site,
              names_from = species,
              values_from = occupancy,
              names_prefix = "sp") %>% 
  select(starts_with("sp")) %>% 
  data.matrix() %>% 
  t()

m_co <- tcrossprod(m)

tibble(y = c(m_co),
       sp1 = rep(1:nrow(m_co), ncol(m_co)),
       sp2 = rep(1:ncol(m_co), each = nrow(m_co)),
       alpha = c(list_df[[2]])) %>% 
  ggplot(aes(x = sp1,
             y = rev(sp2),
             fill = alpha)) +
  geom_raster() +
  geom_text(aes(label = round(y, 2))) +
  MetBrewer::scale_fill_met_c("Hiroshige",
                              direction = -1)
