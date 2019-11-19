get = 'wget http://landsat-pds.s3.amazonaws.com/c1/L8/042/034/LC08_L1TP_042034_20170616_20170629_01_T1/LC08_L1TP_042034_20170616_20170629_01_T1_B4.TIF'

system(get)


ls = stack('LC08_L1TP_042034_20170616_20170629_01_T1_B4.TIF')


ls_df = as.data.frame(ls, xy = TRUE)
head(ls_df)

## ggplot 

ggplot() +
  geom_raster(data = ls_df, aes(x = x, y = y, fill = LC08_L1TP_042034_20170616_20170629_01_T1_B4)) +
  coord_quickmap()
