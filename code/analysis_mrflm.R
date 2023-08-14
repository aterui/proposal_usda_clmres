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
  filter(n > floor(nrow(df_site) * 0.2)) %>% 
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


# plot --------------------------------------------------------------------

g <- graph_from_adjacency_matrix(sA,
                                 mode = "lower",
                                 weighted = "weight")

E(g)$sign <- ifelse(E(g)$weight > 0, "Plus", "Minus")

gnet <- ggraph::ggraph(g, layout = "circle") +
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


# energy landscape --------------------------------------------------------

## binary matrix of presence absence
list_spm <- lapply(c(0, seq_len(length(sp))),
                   function(x) {
                     
                     if (x == 0) {
                       m <- rep(0, length(sp))
                     } else {
                       cbn <- combn(length(sp), x)
                       m <- apply(cbn, MARGIN = 2,
                                  function(i) {
                                    v <- rep(0, length(sp))
                                    v[i] <- 1
                                    return(v)
                                  })
                     }
                     
                     return(as.matrix(t(m)))
                   })

list_one <- lapply(list_spm,
                   function(X) {
                     X[X == 0] <- -1
                     return(X)
                   })

## diagonal elements
## non-species-association factors
Xp <- c(1, 0, 0, 0)
o <- drop(Xp %*% B)
diag(sA) <- o

energy <- lapply(list_spm,
                 function(X) exp(-rowSums(X %*% sA))) %>% 
  unlist()


## neighbor analysis
tictoc::tic()
ncore <- parallel::detectCores() - 1
cl <- parallel::makeCluster(ncore)
doSNOW::registerDoSNOW(cl)

nid0 <- cumsum(lapply(list_spm, nrow))
df_nb <- foreach(i = 1:5,
                 .combine = rbind) %dopar% {
                   
                   m <- list_one[[i + 1]] %*% t(list_one[[i]])
                   list_nb <- apply(m, 2,
                                    FUN = function(x) which(x == length(sp) - 2) + nid0[i],
                                    simplify = FALSE)
                   
                   nnb <- length(list_nb[[1]]) # number of neiboughs to each node
                   from <- rep(if (i == 1) 1 else (nid0[i - 1] + 1):nid0[i],
                               each = nnb)
                   cout <- data.table::data.table(from, to = unlist(list_nb), k = i)
                   return(cout)
                 }

parallel::stopCluster(cl)
tictoc::toc()
gc()

graph <- graph_from_data_frame(df_nb[,c("from", "to")], directed = FALSE)
V(graph)$energy <- energy[seq_len(length(V(graph)))]

minima <- sapply(seq_len(length(V(graph))),
                 function(x) {
                   gap <- energy[x] - energy[neighbors(graph, x)]
                   y <- all(gap < 0)
                   return(y)
                 })
# 
# # adj <- igraph::graph.adjlist(nb, mode = "total")
# # V(adj)$energy <- energy
# # 
# # 
# # glay <- create_layout(adj, layout = "lgl")
# # glay$y <- log(glay$energy)
# # 
# # ggraph(glay) +
# #   geom_edge_link(alpha = 0.5) +
# #   geom_node_point(aes(color = energy,
# #                       size = energy))
