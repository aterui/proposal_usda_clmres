
# delineate watersheds ----------------------------------------------------

watershed <- function(str_grid,
                      f_dir,
                      f_acc,
                      outlet,
                      snap_dist = 5,
                      output_dir = "data_fmt",
                      filename = "watershed",
                      file_ext = "gpkg",
                      keep_outlet = FALSE) {
  
  ## temporary files
  message("Saving temporary files...")
  v_name <- tempdir() %>% 
    paste(c("strg.tif",
            "outlet.shp",
            "outlet_snap.shp",
            "upa.tif",
            "dir.tif",
            "wsd.tif"),
          sep = "\\")
  
  terra::writeRaster(str_grid,
                     filename = v_name[str_detect(v_name, "strg")],
                     overwrite = TRUE)
  
  terra::writeRaster(f_dir,
                     filename = v_name[str_detect(v_name, "dir")],
                     overwrite = TRUE)
  
  terra::writeRaster(f_acc,
                     filename = v_name[str_detect(v_name, "upa")],
                     overwrite = TRUE)
  
  sf::st_write(outlet,
               dsn = v_name[str_detect(v_name, "outlet.shp")],
               append = FALSE)
  
  ## snapping
  message("Snap outlet points to the nearest stream grid...")
  whitebox::wbt_jenson_snap_pour_points(pour_pts = v_name[str_detect(v_name, "outlet\\.")],
                                        streams = v_name[str_detect(v_name, "strg")],
                                        output = v_name[str_detect(v_name, "outlet_snap")],
                                        snap_dist = snap_dist)
  
  ## delineation
  message("Delineate watersheds...")
  whitebox::wbt_unnest_basins(d8_pntr = v_name[str_detect(v_name, "dir")],
                              pour_pts = v_name[str_detect(v_name, "outlet_snap")],
                              output = v_name[str_detect(v_name, "wsd")])
  
  ## vectorize
  message("Vectorize raster watersheds...")
  
  sf_wsd <- list.files(path = tempdir(),
                       pattern = "wsd",
                       full.names = TRUE) %>%
    lapply(terra::rast) %>%
    lapply(stars::st_as_stars) %>%
    lapply(sf::st_as_sf,
           merge = TRUE,
           as_points = FALSE) %>%
    dplyr::bind_rows() %>%
    sf::st_transform(crs = 5070) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(site_id = sum(dplyr::c_across(cols = ends_with("tif")),
                                na.rm = TRUE)) %>%
    dplyr::select(site_id) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(area = units::set_units(st_area(.), "km^2")) %>%
    dplyr::group_by(site_id) %>%
    dplyr::slice(which.max(area)) %>% # remove duplicates by outlet
    dplyr::ungroup() %>%
    dplyr::relocate(site_id, area) %>%
    dplyr::arrange(site_id) %>%
    sf::st_transform(crs = 4326)
  
  if (!(output_dir %in% list.files(".")))
    dir.create(output_dir)
  
  sf::st_write(sf_wsd,
               dsn = paste0(output_dir,
                            "/",
                            filename,
                            ".",
                            file_ext),
               append = FALSE)
  
  if (keep_outlet) {
    outlet_snap <- sf::st_read(dsn = v_name[str_detect(v_name, "outlet_snap")])
    sf::st_write(outlet_snap,
                 dsn = paste0(output_dir,
                              "/",
                              "outlet_snap",
                              ".",
                              file_ext),
                 append = FALSE)    
  }
  
  ## remove temporary files
  message("Removing temporary files...")
  files <- list.files(tempdir(), full.names = T)
  cl <- call("file.remove", files)
  bools <- suppressWarnings(eval(cl, envir = parent.frame()))
  # message(paste0("Following files were removed from a temporary dirctory: ",
  #                files[bools]))
}

# arc2d8 ------------------------------------------------------------------

arc2d8 <- function(x) {
  ## whitebox uses flow pointer D8
  # 64, 128,  1
  # 32,   0,  2 
  # 16,   8,  4
  
  ## ArcGIS uses D8 algorithm
  # 32, 64, 128
  # 16,  0,   1 
  #  8,  4,   2
  
  # check class()
  if(class(x) != "SpatRaster") stop("Provide data in class 'SpatRaster'")
  
  # begin with northeast through north
  fdir_arc <- c(0, 2^(0:7), 247, 255)
  fdir_d8 <- c(0, fdir_arc[3:9], fdir_arc[2], NA, NA)
  y <- terra::subst(x, from = fdir_arc, to = fdir_d8)
  
  return(y)  
}

# from matrix to vector ---------------------------------------------------

m2v <- function(x) {
  
  if (is.null(colnames(x))) {
    index_from <- rep(1:dim(x)[1], times = dim(x)[2])
    index_to <- rep(1:dim(x)[2], each = dim(x)[1])
  } else {
    index_from <- rep(colnames(x), times = dim(x)[2])
    index_to <- rep(rownames(x), each = dim(x)[1])
  }
  
  return(dplyr::tibble(value = c(x),
                       from = index_from,
                       to = index_to))
  
}


# wrapper jags ------------------------------------------------------------

geojags <- function(formula,
                    data,
                    dm,
                    init,
                    n_chain = 3,
                    adapt = 1000,
                    sample = 1000,
                    burnin = floor(sample / 2),
                    thin = 3) {
  
  ## parameters ####
  para <- c("b",
            "sigma",
            "theta",
            "lambda")
  
  ## model file ####
  model_text <- "
  model {
    ninfo <- 0.1
    df0 <- 6
    
    # prior -------------------------------------------------------------------
    
    for (k in 1:Nb) {
      b[k] ~ dnorm(0, ninfo)
    }
    
    for (k in 1) {
      tau[k] ~ dscaled.gamma(2.5, df0)  
      sigma[k] <- pow(tau[k], -0.5)
    }
    
    theta ~ dunif(0, 10)
    lambda ~ dnorm(0, ninfo)T(0, 1)
    
    # likelihood --------------------------------------------------------------
    
    for(i in 1:Nsample) {
      Y[i] ~ dnorm(mu[i], tau[1])
      mu[i] <- inprod(X[i, ], b[]) + lambda * inprod(Q[i, ], y[])
    }
    
    for (i in 1:Nsample) {
      S[i, i] <- 0
      Q[i, 1:Nsample] <- S[i, ] / (sum(S[i, ]) + step(z[i]))
      z[i] <- -sum(S[i, ])
      for (j in (i + 1):Nsample) {
        S[i, j] <- exp(-theta * m_dist[i, j])# * m_w[i, j]
        S[j, i] <- S[i, j]
      }
    }
    
  }
  
  data {
    
    y <- Y
    
    for (n in 1:Ncomb) {
      m_dist[From[n], To[n]] <- Distance[n]
      m_w[From[n], To[n]] <- W[n]
    }
    
  }"
  
  model_path <- paste(tempdir(), "model.R", sep = "\\")
  write(model_text, model_path)
  m <- runjags::read.jagsfile(model_path)
  
  ## data ####
  
  if(missing(init)) {
    inits <- replicate(n_chain,
                       list(.RNG.name = "base::Mersenne-Twister",
                            .RNG.seed = NA),
                       simplify = FALSE)
    
    for (j in 1:n_chain) inits[[j]]$.RNG.seed <- (j - 1) * 10 + 1
  }
  
  Y <- model.frame(formula, data)[,1]
  X <- model.matrix(formula, data)
  
  
  # jags --------------------------------------------------------------------
  d_jags <- list(
    # regression
    Y = Y,
    X = X,
    Nsample = length(Y),
    Nb = ncol(X),
    # spatial component
    From = dm$from,
    To = dm$to,
    Distance = dm$distance, # network distance
    W = dm$c_flow, # flow connectance
    Ncomb = nrow(dm))
  
  ## run jags ####
  post <- runjags::run.jags(m$model,
                            monitor = para,
                            data = d_jags,
                            n.chains = n_chain,
                            inits = inits,
                            method = "parallel",
                            burnin = burnin,
                            sample = sample,
                            adapt = adapt,
                            thin = thin,
                            n.sims = n_chain,
                            module = "glm")
  
  return(MCMCvis::MCMCsummary(post$mcmc))
}


# min matrix --------------------------------------------------------------

symmetrize <- function(X, method = "min") {
  tX <- t(X)
  l <- X[lower.tri(X)]
  u <- tX[lower.tri(tX)]
  
  if (method == "mean")
    y <- apply(cbind(l, u), MARGIN = 1,
               FUN = mean)
  
  if (method == "min")
    y <- apply(cbind(l, u), MARGIN = 1,
               FUN = function(x) x[which.min(abs(x))])
  
  if (method == "max")
    y <- apply(cbind(l, u), MARGIN = 1,
               FUN = function(x) x[which.max(abs(x))])
  
  M <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
  M[lower.tri(M)] <- y
  M <- M + t(M)
  return(M)
}


# ela function ------------------------------------------------------------

local_energy <- function(N, A, ncore = 4) {
  
  library(foreach)
  # possible number of community states
  K <- 2L ^ N
  pw2 <- as.integer(2L ^ seq(0, N - 1, by = 1))
  
  n_per_gr <- floor(K / ncore)
  gr <- seq_len(K) %/% n_per_gr + 1
  gr <- sort(ifelse(gr > ncore, sample(ncore), gr))
  
  # parallel region
  cl <- parallel::makeCluster(ncore)
  doSNOW::registerDoSNOW(cl)
  
  cout <- foreach(g = 1L:max(gr),
                  .combine = c) %dopar% {
                    
                    int <- as.integer(which(gr == g) - 1)
                    sapply(int,
                           function(i) {
                             x <- as.integer(intToBits(i)[1L:N])
                             y <- -sum(drop((A * x) %*% x))
                             return(y)
                           })
                  }
  
  parallel::stopCluster(cl)
  gc(); gc()
  
  return(cout)
}

local_minima <- function(N, h, ncore = 4) {
  
  library(foreach)
  # possible number of community states
  K <- 2L ^ N
  pw2 <- as.integer(2L ^ seq(0, N - 1, by = 1))
  
  n_per_gr <- floor(K / ncore)
  gr <- seq_len(K) %/% n_per_gr + 1
  gr <- sort(ifelse(gr > ncore, sample(ncore), gr))
  
  # parallel region
  cl <- parallel::makeCluster(ncore)
  doSNOW::registerDoSNOW(cl)
  
  bl_min <- foreach(g = 1L:max(gr),
                    .combine = c) %dopar% {
                      
                      int <- as.integer(which(gr == g) - 1)
                      
                      sapply(int, function(i) {
                        x <- as.integer(intToBits(i)[1L:N]) * (-1L)
                        x[x == 0L] <- 1L
                        nei <- i + pw2 * x
                        all(h[nei + 1L] > h[i + 1L])
                      })                        
                    }
  
  dt_nei <- foreach(g = 1L:max(gr),
                    .combine = rbind) %dopar% {
                      
                      int <- as.integer(which(gr == g) - 1)
                      
                      v_nei <- c(sapply(int, function(i) {
                        x <- as.integer(intToBits(i)[1L:N]) * (-1L)
                        x[x == 0L] <- 1L
                        nei <- i + pw2 * x + 1L
                        return(nei)
                      }))
                      
                      dt0 <- data.table::data.table(from = rep(int, each = N) + 1L,
                                                    to = v_nei,
                                                    int = rep(int, each = N))
                      return(dt0[to > from, ])
                    }
  
  parallel::stopCluster(cl)
  gc(); gc()
  
  index_min <- which(bl_min == 1)
  attributes(index_min)$neighbor <- dt_nei
  
  return(index_min)
}

tipping <- function(N, h, ncore = 4) {
  # possible number of community states
  K <- 2L ^ N
  pw2 <- as.integer(2L ^ seq(0, N - 1, by = 1))
  
  n_per_gr <- floor(K / ncore)
  gr <- seq_len(K) %/% n_per_gr + 1
  gr <- sort(ifelse(gr > ncore, sample(ncore), gr))
  
  # parallel region
  cl <- parallel::makeCluster(ncore)
  doSNOW::registerDoSNOW(cl)
  
  ## neg2: T/F node with >= 2 negative gap neighbors
  neg2 <- foreach(g = 1L:max(gr),
                  .combine = c) %dopar% {
                    
                    int <- as.integer(which(gr == g) - 1)
                    
                    sapply(int, function(i) {
                      x <- as.integer(intToBits(i)[1L:N]) * (-1L)
                      x[x == 0L] <- 1L
                      nei <- i + pw2 * x
                      sum(h[nei + 1L] < h[i + 1L]) > 1
                    })                        
                  }  
  
  # tipping points
  ## convert neg2 index to integer
  int <- which(neg2) - 1
  dt_to <- foreach(k = int,
                   .combine = rbind) %dopar% {
    
    ### identify nodes with negative gradients
    x <- as.integer(intToBits(k)[1L:N]) * (-1L)
    x[x == 0L] <- 1L
    nb0 <- k + pw2 * x
    int_neg <- nb0[h[nb0 + 1] < h[k + 1]]
    
    ### identify ultimate local minima
    to <- sapply(int_neg, function(i) {
      repeat {
        x <- as.integer(intToBits(i)[1L:N]) * (-1L)
        x[x == 0L] <- 1L
        nei0 <- i + pw2 * x
        h_nei <- h[nei0 + 1]
        h_i <- h[i + 1]
        if (all(h_nei > h_i)) break
        i <- nei0[which.min(h_nei)]
      }
      return(i + 1)
    })
    
    return(data.table::data.table(node = k + 1,
                                  n_to = length(unique(to)),
                                  to = unique(to),
                                  h_ridge = h[k + 1],
                                  h_minima = h[unique(to)])
           )
  }
  
  parallel::stopCluster(cl)
  gc(); gc()
    
  return(dt_to)
}


# gibbs sampling ----------------------------------------------------------

gibbs <- function(s,
                  h,
                  neighbor,
                  attempt = 10,
                  magnitude = 1,
                  freq = 10) {
  
  dt_nei <- neighbor
  e <- exp(-h)
  
  list_out <- foreach(s = s, .combine = rbind) %do% {
    
    # initial basin of attraction
    m <- s
    repeat {
      v_nei <- c(unlist(dt_nei[from == m, "to"],
                        use.names = FALSE), 
                 unlist(dt_nei[to == m, "from"],
                        use.names = FALSE)
                 )
      if (all(h[v_nei] > h[m])) break
      m <- v_nei[which.min(h[v_nei] - h[m])]
    }
    minima <- m
    
    # perturbation attempts    
    org <- s
    foreach (j = seq_len(attempt),
             .combine = rbind) %do% {
               
               i <- 1
               while(i <= freq) {
                 v_nei <- c(unlist(dt_nei[from == org, "to"],
                                   use.names = FALSE), 
                            unlist(dt_nei[to == org, "from"],
                                   use.names = FALSE)
                            )
                 
                 p <- e[v_nei] / (e[v_nei] + e[org])
                 
                 index <- sample(x = length(p), size = magnitude, replace = T)
                 z <- rbinom(n = magnitude, size = 1, prob = p[index])
                 
                 if (any(z == 1)) {
                   trans_index <- index[min(which(z == 1))]
                   org <- v_nei[trans_index]
                 } else {
                   org <- org
                 }
                 
                 i <- i + 1
               }
               
               x <- org
               repeat {
                 v_nei <- c(unlist(dt_nei[from == x, "to"]), 
                            unlist(dt_nei[to == x, "from"]))
                 names(v_nei) <- NULL 
                 
                 if (all(h[v_nei] > h[x])) break
                 
                 x <- v_nei[which.min(h[v_nei] - h[x])]
               }
               
               return(data.table(initial_state = s,
                                 minima_from = minima,
                                 minima_to = x))
             }
  }
  
  return(list_out)   
}

