# setup -------------------------------------------------------------------

rm(list = ls()); gc()
source("code/library.R")
source("code/function.R")


# data --------------------------------------------------------------------

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


# tictoc::tic()
# h <- local_energy(N = length(sp), A = sA)
# tictoc::toc()
# 
# tictoc::tic()
# m <- local_minima(N = length(sp), h = h)
# tictoc::toc()

# tictoc::tic()
# graph <- graph_from_data_frame(attr(m, "neighbor"), directed = FALSE)
# tictoc::toc()

# saveRDS(list(log_energy = log_energy, minima = m, graph = graph),
#         "output/data_ela.rds")


## energy landscape: agriculture impact

v_x <- seq(0.05, 0.95, length = 10)
xq <- (boot::logit(v_x) - mean(boot::logit(df_site$frac_agri))) / sd(boot::logit(df_site$frac_agri))

cout <- foreach(u = seq_len(length(xq)), .combine = bind_rows) %do% {
  print(u)
  Xp <- c(1, 0, xq[u], 0)
  o <- drop(Xp %*% B)
  diag(sA) <- o
  
  message("calculating local energy...")
  h <- local_energy(N = length(sp), A = sA, ncore = 10)
  
  message("calculating local minima...")
  m <- local_minima(N = length(sp), h = h, ncore = 10)
  mu_gap <- mean(abs(h - min(h[m])))
  
  message("calculating tipping points...")
  dt_tip <- tipping(N = length(sp), h = h, ncore = 10)
  dt_tip <- dt_tip[n_to > 1,]
  dt_tip$gap <- dt_tip$h_ridge - dt_tip$h_minima
  
  combo <- combn(m, 2)
  height <- sapply(seq_len(ncol(combo)), function(v) {
    x <- combo[, v]
    dt_sub <- dt_tip[to %in% x,]
    ridges <- as.numeric(names(which(table(dt_sub$node) == 2)))
    dt_ridge <- dt_sub[node %in% ridges,]
    tapply(dt_ridge$gap, dt_ridge$to, min)
  })
  
  phi <- mean(height)
  
  list_out <- list(id = u,
                   mu_gap = mu_gap,
                   phi = phi,
                   n_minima = length(m),
                   frac_agri = v_x[u],
                   minima = list(c(m)),
                   minma_energy = list(h[m]))
  
  # message("gibbs sampling...")
  # trans <- gibbs(s = m,
  #                h = h,
  #                neighbor = attributes(m)$neighbor,
  #                attempt = 10,
  #                freq = 50,
  #                magnitude = 1000)
  # 
  # trans$z <- ifelse(trans$minima_from == trans$minima_to, 1, 0)
  # v_phi <- tapply(trans$z, trans$initial_state, mean)
  # 
  # list_out <- list(id = u,
  #                  phi = v_phi,
  #                  n_minima = length(m),
  #                  frac_agri = v_x[u],
  #                  minima = list(c(m)),
  #                  minma_energy = list(h[m]))
  
  return(list_out)
}

saveRDS(cout, "output/data_phi.rds")

g_phi <- cout %>%
  group_by(id) %>%
  summarize(phi = mean(phi),
            frac_agri = unique(frac_agri)) %>%
  ggplot(aes(y = phi,
             x = frac_agri * 100)) +
  geom_point(size = 4) +
  geom_line() +
  labs(y = expression("Resilience"~~phi),
       x = "% Agriculture") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

ggsave(g_phi,
       filename = "output/fig_phi.pdf",
       width = 6, height = 5)

# subgraph ----------------------------------------------------------------
# 
# G <- readRDS("output/data_ela.rds")
# 
# minima <- G$minima
# log_e <- G$log_energy
# e <- exp(-log_e)
# dt_nei <- attributes(minima)$neighbor
# 
# subv <- sapply(minima,
#                function(x) {
#                  v <- shortest_paths(G$graph,
#                                      from = x,
#                                      to = minima,
#                                      mode = "out",
#                                      output = "vpath")
# 
#                  unique(unlist(v$vpath))
#                })
# 
# subv <- sort(unique(unlist(subv)))
# 
# subg <- subgraph(G$graph, vids = subv)
# V(subg)$log_e <- log_e[subv]
# V(subg)$minima <- as.numeric(names(V(subg))) %in% minima
# V(subg)$richness <- sapply(as.numeric(names(V(subg))),
#                            function(x) sum(as.integer(intToBits(x - 1))))
# 
# lo <- create_layout(subg, layout = "linear")
# min_id <- which(as.numeric(names(V(subg))) %in% minima)
# 
# lo$x <- V(subg)$richness
# lo$y <- log_e[subv]
# 
# g_ela <- ggraph(lo) +
#   geom_edge_link(alpha = 0.25) +
#   geom_node_point(aes(color = richness,
#                       shape = minima),
#                   size = 7,
#                   alpha = 0.75) +
#   MetBrewer::scale_color_met_c("Hiroshige",
#                                direction = -1) +
#   labs(color = "Richness",
#        y = "Community Energy E") +
#   guides(shape = "none") +
#   theme_classic() +
#   theme(axis.title.y = element_text(size = 24),
#         axis.text.y = element_text(size = 20),
#         axis.line.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.line.x.bottom = element_blank(),
#         axis.line.x.top = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 24),
#         legend.key.size = unit(2, 'cm'),
#         legend.position = c(0.2, 0.3))
# 
# ggsave(g_ela,
#        filename = "output/fig_ela_mrb.pdf",
#        height = 8, width = 9)
# 
# plot --------------------------------------------------------------------

g <- graph_from_adjacency_matrix(sA,
                                 mode = "lower",
                                 weighted = "weight")

E(g)$sign <- ifelse(E(g)$weight > 0, "Plus", "Minus")

lo <- create_layout(g, layout = "linear", circular = TRUE)
gnet <- ggraph::ggraph(lo) +
  geom_edge_arc(aes(alpha = abs(weight),
                    color = sign),
                width = 1) +
  coord_fixed() +
  geom_node_point(size = 5) +
  scale_edge_color_manual(values = c(`Plus` = "steelblue",
                                     `Minus` = "salmon")) +
  theme_void() +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "cm"),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1,
                             unit = "cm")) +
  guides(edge_color = "none",
         edge_alpha = "none")

ggsave(gnet, filename = "output/fig_fish_network.pdf",
       width = 10,
       height = 8)


