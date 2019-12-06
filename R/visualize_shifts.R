library(raster)
library(ggplot2)
taxon = "Crotalus horridus"
files = list.files('data', pattern = taxon, full.names=T)
print(files)
gri.files = files[grep('.gri', files)]
gri.files


mods = stack(gri.files)

plot(mods)

##Plotting shifts
shift26 = mods[[2]] - mods[[1]]

shift26_df = as.data.frame(shift26, xy=T)
(base26 = ggplot() +
    geom_raster(data = shift26_df, 
                aes(x = x, y = y, fill = layer)) + 
    coord_quickmap() +
    theme_bw() + 
    scale_fill_gradientn(colours=c('navy', 'white', 'darkred'), na.value = "black"))



shift85 = mods[[3]] - mods[[1]] # future 2.6 scenario minus current prediction

shift85_df = as.data.frame(shift85, xy=T)

(base85 = ggplot() +
    geom_raster(data = shift85_df, 
                aes(x = x, y = y, fill = layer),
                show.legend=F) + 
    coord_quickmap() +
    theme_bw() + 
    scale_fill_gradientn(colours=c('navy', 'white', 'darkred'), na.value = "black") +
    theme_void())

ext = extent(-74, -69, 40, 45)
shift85.NE = crop(shift85, ext)
shift85.NE_df = as.data.frame(shift85.NE, xy=T)

(base85.NE = ggplot() +
    geom_raster(data = shift85.NE_df, 
                aes(x = x, y = y, fill = layer)) + 
    coord_quickmap() +
    theme_bw() + 
    scale_fill_gradientn(colours=c('navy', 'white', 'darkred'), na.value = "black"))
