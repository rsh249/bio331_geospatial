#spatial bias demo

library(raster)
library(ggplot2)
library(rasterExtras)
library(RSpatial)
library(spocc)
library(dplyr)

wc = getData('worldclim', var='bio', res = 5, path='tmp/')
ext = extent(-125, -55, 20, 60)
wc = crop(wc, ext)

wc_df = as.data.frame(wc, xy=TRUE)





taxon = 'Vaccinium angustifolium'

spdist <- occ(query=taxon, limit=6500)

sp_df = occ2df(spdist)

sp_df = sp_df %>% filter(longitude>=ext[1], 
                         longitude<=ext[2], 
                         latitude>=ext[3], 
                         latitude <=ext[4]) #filter points to study area


gk = gkde(wc[[1]], 
          sp_df[,c('longitude', 'latitude')], 
          parallel=T, 
          nclus = 12, 
          dist.method='Haversine', 
          maxram=20, 
          bw=50)

gk_df=as.data.frame(gk, xy=TRUE)
ggplot() +
  geom_raster(data = gk_df, aes(x = x, y = y, fill = layer)) +
  #geom_point(data=sp_df, aes(x=longitude, y=latitude), col='green', cex=0.1) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=c('darkred', 'grey', 'navy'),
                       na.value = "black")

dplot = ggplot() +
  geom_raster(data = gk_df, aes(x = x, y = y, fill = layer^0.1)) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=c('darkred', 'grey', 'navy'),
                       na.value = "black")

ggsave(dplot, 
       filename = paste("figures/", taxon, "_density.png", sep=""),
       height=7.25, width = 7.25, units='in',
       dpi = 300)

# distance to occ visualization

d2pt = dist2point(wc[[1]], 
                  sp_df[,c('longitude', 'latitude')], 
                  parallel=TRUE, nclus = 12, 
                  dist.method='Haversine', 
                  maxram = 20)

d2_df = as.data.frame(d2pt, xy=T)

ggplot() +
  geom_raster(data = d2_df, aes(x = x, y = y, fill = sqrt(layer))) +
  #geom_point(data=sp_df, aes(x=longitude, y=latitude), col='green', cex=0.2) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=c('darkred', 'grey', 'navy'),
                       na.value = "black")

#spatial thinning
occ2thin = poThin(
  df = sp_df[c('longitude', 'latitude')],
  spacing = 25,
  dimension = nrow(sp_df),
  lon = 'longitude',
  lat = 'latitude'
)

sp_df = sp_df[-occ2thin,]

#test thinned data
gk = gkde(wc[[1]], 
          sp_df[,c('longitude', 'latitude')], 
          parallel=T, 
          nclus = 12, 
          dist.method='Haversine', 
          maxram=20, 
          bw=50)

gk_df=as.data.frame(gk, xy=TRUE)
ggplot() +
  geom_raster(data = gk_df, aes(x = x, y = y, fill = layer)) +
  #geom_point(data=sp_df, aes(x=longitude, y=latitude), col='green', cex=0.1) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=c('darkred', 'grey', 'navy'),
                       na.value = "black")
