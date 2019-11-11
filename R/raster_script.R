# Downloading and working with raster data
library(raster)

clim = getData('worldclim', var='bio', res=5, path='tmp')


ext = extent(-74, -69, 40, 45)
#crop
c2 = crop(clim, ext)
plot(c2[[1]]) # Basic plotting

library(ggplot2)

c2_df = as.data.frame(c2, xy = TRUE)
head(c2_df)

## ggplot 

ggplot() +
  geom_raster(data = c2_df, aes(x = x, y = y, fill = bio1)) +
  coord_quickmap()

base = ggplot() +
  geom_raster(data = c2_df, aes(x = x, y = y, fill = bio1/10)) +
  coord_quickmap()

#modifying base

base + theme_bw()

# test colors

base + 
  theme_bw() + 
  scale_fill_gradientn(colours=c('navy', 'white', 'darkred'), na.value = "black")


library(viridis)
base + 
  theme_bw() + 
  scale_fill_gradientn(colours=viridis(10), na.value = "black")


# distribution of data in rasters

ggplot() +
  geom_histogram(data = c2_df, aes(x = bio1))

