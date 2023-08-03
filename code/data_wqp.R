
# setup -------------------------------------------------------------------

source(here::here("code/library.R"))


# sample mask layer -------------------------------------------------------

wgs84_mrb <- st_read(here::here("data_raw/epsg4269_huc4_mrb.gpkg")) %>% 
  st_transform(crs = 4326)

# nutrient data -----------------------------------------------------------

## Characterirstics
# [1] "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)"                      
# [2] "Organic Nitrogen"                                                                  
# [3] "Ammonia and ammonium"                                                              
# [4] "Nitrite"                                                                           
# [5] "Nitrate"                                                                           
# [6] "Kjeldahl nitrogen"                                                                 
# [7] "Inorganic nitrogen (nitrate and nitrite)"                                          
# [8] "Orthophosphate"                                                                    
# [9] "Phosphorus"                                                                        
# [10] "Nitrogen"                                                                          
# [11] "Inorganic nitrogen (nitrate and nitrite) ***retired***use Nitrate + Nitrite"       
# [12] "Ammonia-nitrogen"                                                                  
# [13] "Total Kjeldahl nitrogen (Organic N & NH3)"                                         
# [14] "Ammonium"                                                                          
# [15] "Ammonia"                                                                           
# [16] "Nitrate + Nitrite"                                                                 
# [17] "Total Nitrogen, mixed forms"                                                       
# [18] "Total Kjeldahl nitrogen"                                                           
# [19] "Total Phosphorus, mixed forms"                                                     
# [20] "Phosphate-phosphorus***retired***use Total Phosphorus, mixed forms"                
# [21] "Nutrient-nitrogen***retired***use TOTAL NITROGEN, MIXED FORMS with speciation AS N"

# nutrient data
df_nu0 <- readWQPdata(characteristicType = "Nutrient",
                      bBox = as.vector(st_bbox(wgs84_mrb)),
                      startDateLo = "2010-01-01",
                      startDateHigh = "2010-12-31",
                      convertType = FALSE)

df_nu <- df_nu0 %>% 
  filter(ActivityMediaName == "Water",
         ResultValueTypeName == "Actual") %>% 
  mutate(value_raw = as.numeric(ResultMeasureValue)) %>% 
  dplyr::select(activity_id = ActivityIdentifier,
                location_id = MonitoringLocationIdentifier,
                characteristic = CharacteristicName,
                sample_fraction = ResultSampleFractionText,
                media = ActivityMediaName,
                date = ActivityStartDate,
                hydro = HydrologicCondition,
                value_raw,
                value_raw_unit = ResultMeasure.MeasureUnitCode) %>% 
  mutate(date = as.Date(date))

df_site <- attributes(df_nu)$siteInfo %>% 
  dplyr::select(location_id = MonitoringLocationIdentifier,
                lat = dec_lat_va,
                lon = dec_lon_va) %>% 
  mutate(lat = as.numeric(lat),
         lon = as.numeric(lon))


# filter ------------------------------------------------------------------

df_no3 <- df_nu %>% 
  filter(characteristic %in% c("Nitrate"),
         value_raw_unit %in% c("mg/L", "mg/l asNO3", "mg/l NO3"),
         sample_fraction == "Dissolved") %>% 
  drop_na(value_raw) %>% 
  mutate(value = value_raw,
         value_unit = "mg/l") %>% 
  left_join(df_site,
            by = "location_id") %>% 
  group_by(lon, lat) %>% 
  slice(which.min(date)) %>% 
  ungroup()
  
df_tp <- df_nu %>%
  filter(characteristic %in% c("Phosphorus",
                               "Total Phosphorus, mixed forms"),
         sample_fraction == "Total",
         value_raw_unit != "% by wt") %>%
  drop_na(value_raw) %>%
  mutate(characteristic = "TP",
         value = case_when(value_raw_unit == "mg/l as P" ~ value_raw,
                           value_raw_unit == "mg/L" ~ value_raw,
                           value_raw_unit == "ug/L" ~ 0.001 * value_raw,
                           value_raw_unit == "ppb" ~ 0.001 * value_raw),
         value_unit = "mg/l") %>%
  left_join(df_site,
            by = "location_id") %>%
  group_by(lon, lat) %>% 
  slice(which.min(date)) %>% 
  ungroup()

df_np <- bind_rows(df_no3, df_tp)

saveRDS(df_np,
        here::here("data_raw/data_np.rds"))
