date <- seq(as.Date("2011-01-01"), as.Date("2023-11-16"), by = "month")
date_df = data.frame(date = date,
                     num = seq(1, length(date)),
                     month = month(date))

raster_sst <- rast("~/stage_M2/1.data/copernicus_data/SST.tif")
coord_df <- crds(raster_sst)

nb_of_month <- dim(raster_sst)[3]

sst_mean <- coord_df

for (i in 1:12){
  col_to_extract <- date_df %>% 
    filter(date >= as.Date("2011-01-01") &
             date <= as.Date("2023-12-01")) %>% 
    filter(month == i) %>% 
    pull(num) %>% 
    as.character()
  
  sst_mean <- raster_sst %>% as.data.frame() %>% 
    select(col_to_extract) %>% 
    rowMeans() %>% 
    scale() %>% 
    as_tibble() %>% 
    bind_cols(sst_mean)
}

colnames(sst_mean)  <- c("Dec", "Nov", "Oct", 
                         "Sep", "Aug", "Jul", 
                         "Jun", "May", "Apr", 
                         "Mar", "Feb", "Jan",
                         "X", "Y")

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

sst_mean %>% 
  pivot_longer(-c(X,Y), values_to = "value", names_to = "month") %>% 
  mutate(month = factor(month, levels = month_names)) %>% 
  ggplot(aes(x=X, y=Y, fill= value)) +
  geom_raster() +
  facet_wrap(~month, ncol = 3) +
  scale_fill_distiller(palette = "Spectral", name = "T°") +
  labs(title = "Scaled temperature by month",
       subtitle = "Average between 2011 & 2021",
       y = "Longitude",
       x = "Latitude") +
  theme_bw()


library(GGally)
sst_mean %>% 
  select(month_names) %>% 
  ggpairs()

