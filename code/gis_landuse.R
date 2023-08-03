
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))
source(here::here("code/function.R"))


# read data ---------------------------------------------------------------

## temporary directory
root <- here::here() %>% 
  str_remove("public-proj_markov-network")

albers_rs_lu <- terra::rast("E:/github/priv-proj_midwest-gis/data_gis/albers_lu.tif")
albers_sf_wsd <- st_read("data_fmt/epsg4326_wsd.gpkg") %>% 
  st_transform(5070)
df_tp <- st_read("data_fmt/epsg4326_point_snap.gpkg") %>% 
  as_tibble()

# raster classification ---------------------------------------------------

## forest: 111-126
albers_rs_forest <- terra::classify(albers_rs_lu,
                                    rcl = rbind(c(110.9, 126.1, 1),
                                                c(2, 500, 0)))

## urban: 50
albers_rs_urban <- terra::classify(albers_rs_lu,
                                   rcl = rbind(c(49.9, 50.1, 1),
                                               c(2, 500, 0)))

## agri: 40
albers_rs_agri <- terra::classify(albers_rs_lu,
                                  rcl = rbind(c(39.9, 40.1, 1),
                                              c(2, 500, 0)))

## multi-layer stake
albers_rs_fua <- terra::rast(list(albers_rs_forest,
                                  albers_rs_urban,
                                  albers_rs_agri))

names(albers_rs_fua) <- c("forest", "urban", "agri")


# extraction --------------------------------------------------------------

df_lu <- exact_extract(albers_rs_fua,
                       albers_sf_wsd,
                       "mean",
                       append_cols = TRUE) %>% 
  rename(frac_forest = mean.forest,
         frac_urban = mean.urban,
         frac_agri = mean.agri) %>% 
  as_tibble()


# join --------------------------------------------------------------------

df_m <- df_tp %>% 
  left_join(df_lu, by = "node_id")

saveRDS(df_m, "data_fmt/data_tpsub_master.rds")

ggplot(df_m,
       aes(x = area,
           y = value)) +
  geom_point() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")
  