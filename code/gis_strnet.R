
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))
source(here::here("code/function.R"))


# read data ---------------------------------------------------------------

## temporary directory
root <- here::here() %>% 
  str_remove("public-proj_markov-network")

## read layers
### strg: grid stream layer
wgs84_sr_strg <- terra::rast("data_raw/epsg4326_stream_grid_50sqkm.tif")

### dir: flow direction
wgs84_sr_dir <- terra::rast("E:/github/priv-proj_midwest-gis/data_gis/epsg4326_dir.tif") %>% 
  arc2d8()

### upa: upstream drainage area
wgs84_sr_upa <- terra::rast("E:/github/priv-proj_midwest-gis/data_gis/epsg4326_upa.tif")

### sf_mask: mask layer
wgs84_sf_mask <- st_read(here::here("data_raw/epsg4269_huc4_mrb.gpkg")) %>% 
  st_transform(crs = 5070) %>% 
  st_buffer(dist = 10000) %>% 
  dplyr::select(NULL) %>% 
  st_transform(crs = 4326) %>% 
  mutate(id = 1)

wgs84_sf_mask0 <- st_read(here::here("data_raw/epsg4269_huc4_mrb.gpkg")) %>% 
  st_transform(crs = 4326) %>% 
  dplyr::select(NULL) %>% 
  mutate(id = 1)
  
### sf_str: vector stream
wgs84_sf_str <- st_read(dsn = "data_raw/epsg4326_channel_50sqkm.shp") %>% 
  st_set_crs(4326) %>% 
  st_join(wgs84_sf_mask) %>% 
  drop_na(id) %>% 
  dplyr::select(NULL)

comp <- wgs84_sf_str %>% 
  st_touches() %>% 
  graph.adjlist() %>% 
  components()

wgs84_sf_str <- wgs84_sf_str %>% 
  mutate(group = comp$membership) %>% 
  filter(group == which.max(comp$csize))

albers_sf_str_buff <- wgs84_sf_str %>% 
  st_transform(5070) %>% 
  st_buffer(dist = 500) %>% 
  st_union() %>% 
  st_as_sf() %>% 
  st_transform(5070) %>%
  mutate(id = row_number())
  
### wgs84_sf_tp: vector point for TP sampling
set.seed(111)
wgs84_sf_tp <- readRDS(here::here("data_raw/data_np.rds")) %>% 
  filter(characteristic %in% c("Phosphorus", "Total Phosphorus, mixed forms")) %>% 
  arrange(date) %>% 
  dplyr::select(-c(activity_id, value_raw, value_raw_unit)) %>% 
  mutate(unique_site = paste(round(lon, 2), round(lat, 2))) %>% 
  group_by(unique_site) %>% 
  slice(which.min(date)) %>% 
  st_as_sf(coords = c("lon", "lat"),
           crs = 4269) %>% 
  st_transform(crs = 5070) %>% 
  st_join(albers_sf_str_buff) %>% 
  st_join(st_transform(wgs84_sf_mask0, crs = 5070)) %>% 
  drop_na(starts_with("id.")) %>% 
  dplyr::select(-starts_with("id.")) %>% 
  sample_n(size = 200) %>% 
  mutate(site_id = row_number()) %>% 
  st_as_sf() %>% 
  st_transform(4326)

saveRDS(wgs84_sf_tp %>% as_tibble(),
        here::here("data_fmt/data_tp_sub.rds"))

# save to tempdir() -------------------------------------------------------

v_name <- tempdir() %>% 
  paste(c("strg.tif",
          "outlet.shp",
          "outlet_snap.shp",
          "upa.tif",
          "dir.tif",
          "wsd.tif"),
        sep = "\\")

terra::writeRaster(wgs84_sr_strg,
                   filename = v_name[str_detect(v_name, "strg")],
                   overwrite = TRUE)

terra::writeRaster(wgs84_sr_dir,
                   filename = v_name[str_detect(v_name, "dir")],
                   overwrite = TRUE)

terra::writeRaster(wgs84_sr_upa,
                   filename = v_name[str_detect(v_name, "upa")],
                   overwrite = TRUE)

# wbt_extract_streams(flow_accum = v_name[str_detect(v_name, "upa")],
#                     output = v_name[str_detect(v_name, "strg")],
#                     threshold = 50)
# 
# x <- terra::rast(v_name[str_detect(v_name, "strg")])
# terra::writeRaster(x, filename = here::here("data_raw/epsg4326_stream_grid_50sqkm.tif"))
# 
# wbt_raster_streams_to_vector(streams = v_name[str_detect(v_name, "strg")],
#                              d8_pntr = v_name[str_detect(v_name, "dir")],
#                              output = "data_raw/epsg4326_channel_50sqkm.shp")

st_write(wgs84_sf_tp,
         dsn = v_name[str_detect(v_name, "outlet.shp")],
         append = FALSE)

list.files(tempdir(), full.names = T)


# gis ---------------------------------------------------------------------

## pour point snap
wbt_jenson_snap_pour_points(pour_pts = v_name[str_detect(v_name, "outlet\\.")],
                            streams = v_name[str_detect(v_name, "strg")],
                            output = v_name[str_detect(v_name, "outlet_snap")],
                            snap_dist = 2)

wgs84_sf_tp_snap <- st_read(dsn = v_name[str_detect(v_name, "outlet_snap")])

wgs84_sfnet <- as_sfnetwork(wgs84_sf_str, directed = FALSE)
wgs84_sfnet %>%
  activate(edges) %>% 
  st_as_sf() %>% 
  st_write("data_fmt/epsg4326_str_mrb_raw.gpkg")

wgs84_sfnetb <- st_network_blend(wgs84_sfnet, wgs84_sf_tp_snap) %>%
  activate(edges) %>%
  mutate(w = units::set_units(edge_length(), "km")) %>% 
  activate(nodes) %>% 
  mutate(node_id = row_number(),
         prefix = ifelse(is.na(site_id),
                         "p_",
                         "o_")) %>% 
  group_by(prefix) %>% 
  mutate(index = paste0(prefix, ifelse(is.na(site_id),
                                       row_number(),
                                       site_id))) %>% 
  ungroup() %>% 
  arrange(node_id)

## network
wgs84_sfnetb %>%
  activate(edges) %>%
  as_tibble() %>%
  st_as_sf() %>%
  st_write(here::here("data_fmt/epsg4326_str_mrb.gpkg"),
           append = FALSE)

wgs84_sf_tp_snap <- wgs84_sfnetb %>% 
  activate(nodes) %>% 
  as_tibble() %>% 
  mutate(node_id = row_number())

st_write(wgs84_sf_tp_snap,
         dsn = v_name[str_detect(v_name, "outlet_snap")],
         append = FALSE)

st_write(wgs84_sf_tp_snap,
         dsn = here::here("data_fmt/epsg4326_point_snap.gpkg"),
         append = FALSE)

## watershed delineation
wbt_unnest_basins(d8_pntr = v_name[str_detect(v_name, "dir")],
                  pour_pts = v_name[str_detect(v_name, "outlet_snap")],
                  output = v_name[str_detect(v_name, "wsd")])

wgs84_sf_wsd <- list.files(path = tempdir(),
                           pattern = "wsd",
                           full.names = TRUE) %>% 
  lapply(terra::rast) %>% 
  lapply(stars::st_as_stars) %>% 
  lapply(sf::st_as_sf,
         merge = TRUE,
         as_points = FALSE) %>%
  bind_rows() %>% 
  st_transform(crs = 5070) %>% 
  rowwise() %>% 
  mutate(node_id = sum(c_across(cols = ends_with("tif")),
                       na.rm = TRUE)) %>% 
  dplyr::select(node_id) %>% 
  ungroup() %>% 
  mutate(area = units::set_units(st_area(.), "km^2")) %>% 
  group_by(node_id) %>% 
  slice(which.max(area)) %>% # remove duplicates by outlet
  ungroup() %>% 
  relocate(node_id, area) %>% 
  arrange(node_id) %>% 
  st_transform(crs = 4326)

st_write(wgs84_sf_wsd,
         dsn = here::here("data_fmt/epsg4326_wsd.gpkg"),
         append = FALSE)

## distance
root_id <- 35 # root node id; currently manual check with QGIS
m_dist <- st_network_cost(wgs84_sfnetb,
                          weights = "w")

site_id <- wgs84_sfnetb %>% 
  activate(nodes) %>% 
  as_tibble() %>%
  pull(index)

v_u <- m_dist[root_id, -root_id]
m_u <- outer(v_u, v_u, FUN = "-")

m_d <- m_dist[-root_id, -root_id]
colnames(m_d) <- site_id[-root_id]; rownames(m_d) <- site_id[-root_id]

m_up <- abs(m_u - m_d) / 2
colnames(m_up) <- site_id[-root_id]; rownames(m_up) <- site_id[-root_id]

m_down <- abs(m_u + m_d) / 2
colnames(m_down) <- site_id[-root_id]; rownames(m_down) <- site_id[-root_id]

x <- m2v(m_d) %>% 
  rename(distance = value)

y <- m2v(m_up) %>% 
  rename(up = value)

z <- m2v(m_down) %>% 
  rename(down = value)

df_dist <- purrr::reduce(list(x, y, z),
                         dplyr::left_join,
                         by = c("from", "to")) %>% 
  mutate(across(.cols = c("distance", "up", "down"),
                \(x) round(x, 2)),
         type = ifelse(str_detect(from, "o_") & str_detect(to, "o_"),
                       yes = "observation",
                       no = "prediction")
         ) %>% 
  arrange(up)

saveRDS(df_dist, here::here("data_fmt/data_distance.rds"))

# mapping -----------------------------------------------------------------
# 
# ggplot() +
#   geom_sf(data = graph %>%
#             activate(edges) %>% 
#             as_tibble() %>% 
#             st_as_sf() %>% 
#             filter(edge_id %in% 1:100)) +
#   geom_sf(data = graph %>%
#             activate(nodes) %>% 
#             as_tibble() %>% 
#             st_as_sf() %>% 
#             filter(edge_id %in% 1:100),
#           color = "red")
  

# clean tempdir() ---------------------------------------------------------

file.remove(list.files(tempdir(), full.names = T))
