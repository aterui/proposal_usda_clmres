# setup -------------------------------------------------------------------

rm(list = ls()); gc()
source("code/library.R")
source("code/function.R")

set.seed(122)
df_site <- readRDS("data_fmt/data_mrb_site.rds") %>% 
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

sA <- symmetrize(A, "min")

Xp <- c(1, 0, 0, 0)
o <- drop(Xp %*% B)
diag(sA) <- o


# energy landscape --------------------------------------------------------
tictoc::tic()
log_energy <- local_energy(N = length(sp), A = sA)
tictoc::toc()

tictoc::tic()
m <- local_minima(N = length(sp), e = log_energy)
tictoc::toc()

tictoc::tic()
graph <- graph_from_data_frame(attr(m, "neighbor"), directed = FALSE)
tictoc::toc()

saveRDS(list(log_energy = log_energy, minima = m, graph = graph),
        "output/data_ela.rds")

# plot --------------------------------------------------------------------
# 
# g <- graph_from_adjacency_matrix(sA,
#                                  mode = "lower",
#                                  weighted = "weight")
# 
# E(g)$sign <- ifelse(E(g)$weight > 0, "Plus", "Minus")
# 
# gnet <- ggraph::ggraph(g, layout = "circle") +
#   geom_edge_arc(aes(alpha = abs(weight),
#                     color = sign),
#                 width = 1) +
#   coord_fixed() +
#   geom_node_point(size = 5) +
#   scale_edge_color_manual(values = c(`Plus` = "steelblue",
#                                      `Minus` = "salmon")) +
#   theme_void() +
#   theme(legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.key.size = unit(1, "cm"),
#         plot.margin = margin(t = 1, r = 1, b = 1, l = 1,
#                              unit = "cm")) +
#   guides(edge_color = "none",
#          edge_alpha = "none")
# 
# ggsave(gnet, filename = "output/fig_fish_network.pdf",
#        width = 10,
#        height = 8)
# 
# 
