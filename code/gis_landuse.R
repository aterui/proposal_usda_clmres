
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))
source(here::here("code/function.R"))


# read data ---------------------------------------------------------------

## temporary directory
root <- here::here() %>% 
  str_remove("public-proj_markov-network")

albers_rs_lu <- terra::rast("F:/gis/nlcd_mrb/NLCD_2011_Land_Cover_L48_20210604_6xx9CG8DGXK8ycYmOKID.tiff")

albers_sf_wsd <- st_read("data_fmt/epsg4326_wsd.gpkg") %>% 
  st_transform(st_crs(albers_rs_lu))

df_tp <- st_read("data_fmt/outlet_snap.gpkg") %>% 
  as_tibble()

# raster classification ---------------------------------------------------

## forest: 41 - 43
albers_rs_forest <- terra::classify(albers_rs_lu,
                                    rcl = rbind(c(40.9, 43.1, 1),
                                                c(0, 100, 0)), )

## agri: 81 - 82
albers_rs_agri <- terra::classify(albers_rs_lu,
                                  rcl = rbind(c(80.9, 82.1, 1),
                                              c(0, 100, 0)))

## wetland: 90, 95
albers_rs_wet <- terra::classify(albers_rs_lu,
                                  rcl = rbind(c(89.9, 95.1, 1),
                                              c(0, 100, 0)))

## multi-layer stake
albers_rs_multi <- terra::rast(list(albers_rs_forest,
                                    albers_rs_agri,
                                    albers_rs_wet))

names(albers_rs_multi) <- c("forest", "agri", "wet")


# extraction --------------------------------------------------------------

df_lu <- exact_extract(albers_rs_multi,
                       albers_sf_wsd,
                       "mean",
                       append_cols = TRUE) %>% 
  rename(frac_forest = mean.forest,
         frac_agri = mean.agri,
         frac_wet = mean.wet) %>% 
  as_tibble()


# join --------------------------------------------------------------------

df_m <- df_tp %>% 
  left_join(df_lu, by = "site_id")

saveRDS(df_m, "data_fmt/data_tpsub_master.rds")

ggplot(df_m,
       aes(x = frac_agri,
           y = value,
           color = frac_wet)) +
  geom_point() +
  #scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")
  