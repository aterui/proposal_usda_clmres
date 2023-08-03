
# setup -------------------------------------------------------------------

rm(list = ls())
source("code/library.R")
source("code/function.R")

## fish data in MW
df_fish <- list.files("E:/github/priv-proj_midwest-data/data_fmt",
                      full.names = TRUE,
                      pattern = "data_fmt_") %>% 
  lapply(read_csv) %>% 
  bind_rows() %>% 
  rename_with(.fn = str_to_lower) %>% 
  mutate(across(.cols = where(is.character),
                .fn = str_to_lower)) %>% 
  dplyr::select(-1)

wgs84_sf_fish <- df_fish %>% 
  st_as_sf(coords = c("lon", "lat")) %>% 
  st_set_crs(4269)

## mask layer
wgs84_sf_mask <- st_read("data_raw/epsg4269_huc4_mrb.gpkg") %>% 
  dplyr::select(NULL) %>% 
  mutate(index = 1)

## vector stream
wgs84_sf_str <- st_read("data_fmt/epsg4326_str_mrb_raw.gpkg")

## land use data
### forest: 111-126, urban: 50, agri: 40
albers_rs_lu <- terra::rast("E:/github/priv-proj_midwest-gis/data_gis/albers_lu.tif")

albers_rs_forest <- terra::classify(albers_rs_lu,
                                    rcl = rbind(c(110.9, 126.1, 1),
                                                c(2, 500, 0)))

albers_rs_urban <- terra::classify(albers_rs_lu,
                                   rcl = rbind(c(49.9, 50.1, 1),
                                               c(2, 500, 0)))

albers_rs_agri <- terra::classify(albers_rs_lu,
                                  rcl = rbind(c(39.9, 40.1, 1),
                                              c(2, 500, 0)))

albers_rs_fua <- terra::rast(list(albers_rs_forest,
                                  albers_rs_urban,
                                  albers_rs_agri))

names(albers_rs_fua) <- c("forest", "urban", "agri")


# format fish data --------------------------------------------------------

df_mrbfish <- st_join(wgs84_sf_fish, wgs84_sf_mask) %>% 
  filter(!is.na(index)) %>% 
  mutate(presence = 1)

wgs84_sf_site0 <- df_mrbfish %>% 
  distinct(siteid, geometry) %>% 
  st_transform(4326) %>% 
  st_as_sf()

wgs84_sf_str_buff <- wgs84_sf_str %>% 
  st_transform(5070) %>% 
  st_buffer(dist = 500) %>% 
  st_union() %>% 
  st_transform(4326) %>% 
  st_as_sf() %>% 
  st_make_valid() %>% 
  mutate(index = 1)

wgs84_sf_site <- wgs84_sf_site0 %>% 
  st_join(wgs84_sf_str_buff) %>% 
  filter(!is.na(index)) %>% 
  mutate(index = row_number())

df_binary <- df_mrbfish %>% 
  as_tibble() %>% 
  arrange(siteid, species) %>% 
  filter(siteid %in% unique(wgs84_sf_site$siteid)) %>% 
  pivot_wider(id_cols = siteid,
              names_from = species,
              names_prefix = "f_",
              values_from = presence,
              values_fill = 0) %>% 
  pivot_longer(cols = starts_with("f_"),
               names_to = "species",
               names_transform = list(species = function(x) str_remove(x, "f_")),
               values_to = "presence")


# delineate watersheds at fish sites --------------------------------------

## strg: grid stream layer
## dir: flow direction
## upa: upstream drainage area
# wgs84_sr_strg <- terra::rast("data_raw/epsg4326_stream_grid_50sqkm.tif")
# wgs84_sr_dir <- terra::rast("E:/github/priv-proj_midwest-gis/data_gis/epsg4326_dir.tif") %>% 
#   arc2d8()
# wgs84_sr_upa <- terra::rast("E:/github/priv-proj_midwest-gis/data_gis/epsg4326_upa.tif")

# watershed(str_grid = wgs84_sr_strg,
#           f_dir = wgs84_sr_dir,
#           f_acc = wgs84_sr_upa, 
#           outlet = wgs84_sf_site,
#           output_dir = "data_fmt",
#           filename = "epsg4326_wsd_fish")

# format river data -------------------------------------------------------

## assign betweeness values
graph_str <- wgs84_sf_str %>% 
  mutate(seg_index = row_number()) %>% 
  st_touches() %>% 
  graph.adjlist(mode = "all") # undirected graph

btw <- betweenness(graph_str,
                   directed = FALSE,
                   cutoff = -1,
                   normalized = TRUE)

str_nearest <- st_nearest_feature(wgs84_sf_site, wgs84_sf_str)

# wgs84_sf_str %>% 
#   ggplot() +
#   geom_sf(aes(color = btw)) +
#   MetBrewer::scale_color_met_c("Hiroshige",
#                                direction = -1)


# environments at fish sites ----------------------------------------------

wgs84_sf_wsd_fish <- st_read("data_fmt/epsg4326_wsd_fish.gpkg")
albers_sf_wsd_fish <- wgs84_sf_wsd_fish %>% 
  st_transform(5070)

df_lu <- exact_extract(albers_rs_fua,
                       albers_sf_wsd_fish,
                       "mean",
                       append_cols = TRUE) %>% 
  rename(frac_forest = mean.forest,
         frac_urban = mean.urban,
         frac_agri = mean.agri) %>% 
  as_tibble()

df_site <- wgs84_sf_site %>% 
  left_join(df_lu,
            by = c("index" = "site_id")) %>% 
  as_tibble() %>% 
  mutate(btw = btw[str_nearest])


# export ------------------------------------------------------------------

saveRDS(df_binary, "data_fmt/data_mrb_fish_binary.rds")
saveRDS(df_site, "data_fmt/data_mrb_site.rds")
