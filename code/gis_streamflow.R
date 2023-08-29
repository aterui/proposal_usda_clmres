
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))
source(here::here("code/function.R"))


# stream gauge site -------------------------------------------------------

gauge_mrb <- c("05290000",
               "05291000",
               "05292000",
               "05293000",
               "05294000",
               "05300000",
               "05301000",
               "05304500",
               "05304995",
               "05305000",
               "05311000",
               "05311150",
               "05313500",
               "05315000",
               "05316500",
               "05316580",
               "05316770",
               "05317000",
               "05317200",
               "05319500",
               "05320000",
               "05320500",
               "05325000",
               "05327000",
               "05330000")

sf_gauge <- dataRetrieval::readNWISsite(gauge_mrb) %>% 
  as_tibble() %>% 
  st_as_sf(coords = c("dec_long_va", "dec_lat_va")) %>% 
  st_set_crs(4269) %>% 
  st_transform(4326) %>% 
  dplyr::select(site_no)


# -------------------------------------------------------------------------

wgs84_sr_strg <- terra::rast("data_raw/epsg4326_stream_grid_50sqkm.tif")
wgs84_sr_dir <- terra::rast("F:/github/priv-proj_midwest-gis/data_gis/epsg4326_dir.tif") %>% 
  arc2d8()
wgs84_sr_upa <- terra::rast("F:/github/priv-proj_midwest-gis/data_gis/epsg4326_upa.tif")

watershed(str_grid = wgs84_sr_strg,
          f_dir = wgs84_sr_dir,
          f_acc = wgs84_sr_upa,
          outlet = sf_gauge,
          filename = "watershed_gauge")

wgs84_gauge_wsd <- st_read("data_fmt/watershed_gauge.gpkg") %>% 
  mutate(area = units::set_units(st_area(.), "km^2"),
         site_no = sf_gauge$site_no)

albers_gauge_wsd <- wgs84_gauge_wsd %>% 
  st_transform(5070)

# landuse
albers_rs_lu <- terra::rast("F:/gis/nlcd_mrb/NLCD_2011_Land_Cover_L48_20210604_6xx9CG8DGXK8ycYmOKID.tiff")

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

df_lu <- exact_extract(albers_rs_multi,
                       albers_gauge_wsd,
                       "mean",
                       append_cols = TRUE) %>% 
  rename(frac_forest = mean.forest,
         frac_agri = mean.agri,
         frac_wet = mean.wet) %>% 
  as_tibble() %>% 
  mutate(area = as.numeric(area))

# flow --------------------------------------------------------------------

df_flow <- foreach(i = seq_len(length(gauge_mrb)),
                   .combine = bind_rows) %do% {

                     df_x <- readNWISdata(siteNumbers = gauge_mrb[i],
                                          parameterCd = "00060",
                                          startDate = "2000-01-01",
                                          endDate = "2015-12-31")
                     
                   }

df_m <- df_flow %>% 
  left_join(df_lu,
            by = "site_no") %>% 
  rename(flow = X_00060_00003,
         unit = X_00060_00003_cd,
         date = dateTime) %>% 
  mutate(unit = "cubic_feet",
         flow_l = flow * 28.3168,
         unit_l = "L")


# analysis ----------------------------------------------------------------

fit <- lme4::lmer(log(flow_l) ~ log(area) * frac_agri + (1 | date),
                  df_m)

beta <- coef(fit)$date
saveRDS(beta, "output/data_flow_beta.rds")
