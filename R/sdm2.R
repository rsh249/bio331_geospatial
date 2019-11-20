# SDM Attempt #1

library(raster)
library(ggplot2)
library(rasterExtras)
library(RSpatial)
library(spocc)
library(dplyr)
library(ENMeval)

#accepting a commandline argument
argv <-commandArgs(TRUE)
if (length(argv) > 0) {
  taxon <- as.numeric(argv[1])
} else {
  print("NO TAXON SUBMITTED. TERMINATING JOB!")
  q('no') #exit without saving
}

# set up predictor data parameters
res=2.5
path='tmp/'
predvars = c(1,2,3,12,16) # subset of bioclim variables as predictors


wc = getData('worldclim', 
             var='bio', 
             res = res, 
             path = path)


wc_future85 = future = getData('CMIP5', 
                             var='bio', 
                             res=res, 
                             rcp = 85, 
                             model='AC', 
                             year=70, 
                             path = path)

wc_future26 = future = getData('CMIP5', 
                               var = 'bio', 
                               res = res, 
                               rcp = 26, 
                               model='AC', 
                               year =70, 
                               path = path)

#subset predictor variables
preds = wc[[predvars]]
preds26 = wc_future26[[predvars]]
preds85 = wc_future85[[predvars]]

# crop to study area
ext = extent(-125, -55, 20, 60) 
wc = crop(wc, ext)
wc_df = as.data.frame(wc, xy=TRUE) # data frame for plotting

#downloading
spdist <- occ(query=taxon, limit=6500) # check limit for your species
sp_df = occ2df(spdist)

#filtering
sp_df = sp_df %>% filter(longitude >= ext[1], 
                         longitude <= ext[2], 
                         latitude >= ext[3], 
                         latitude <= ext[4]) #dplyr filter points to study area

#thinning
occ2thin = poThin(
  df = sp_df[c('longitude', 'latitude')],
  spacing = 25, #minimum distance between points in thinned data
  dimension = nrow(sp_df),
  lon = 'longitude',
  lat = 'latitude'
)

sp_df = sp_df[-occ2thin,] #thin using index returned from occ2thin

#plot
dist_plot = ggplot() +
  geom_raster(data = wc_df, aes(x = x, y = y, fill = bio12)) +
  geom_point(data=sp_df, aes(x=longitude, y=latitude), col='green', cex=0.1) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=c('darkred', 'grey', 'navy', 'green'),
                       na.value = "black")
ggsave(dist_plot, 
       filename = paste("figures/", taxon, "_data.png", sep=""),
       height=7.25, width = 7.25, units='in',
       dpi = 300)


# ENMevaluate

eval = ENMevaluate(occ=sp_df[,c('longitude', 'latitude')], 
                   env = preds, 
                   method='randomkfold', 
                   kfolds=100, 
                   parallel=TRUE, 
                   numCores = 12, 
                   fc=c("L", "Q", "LQ", "LQH", "LQHPT"), 
                   RMvalues=seq(0.5, 4, 0.5), 
                   rasterPreds=T)


# predict

bestmod = which(eval@results$AICc==min(eval@results$AICc))

pr = predict(preds, eval@models[[bestmod]], type = 'cloglog')
pr_df = as.data.frame(pr, xy=T)

#heatmap
max_plot = ggplot() +
  geom_raster(data = pr_df, aes(x = x, y = y, fill = layer)) +
  geom_point(data=sp_df, aes(x=longitude, y=latitude), col='red', cex=0.05) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=viridis::viridis(99),
                       na.value = "black")
ggsave(max_plot, 
       filename = paste("figures/", taxon, "_maxent.png", sep=""),
       height=7.25, width = 7.25, units='in',
       dpi = 300)
writeRaster(pr, 
            filename = paste('data/', taxon, '_model'), 
            overwrite=TRUE)
# Binary Range Map

#extract model estimated suitability for occurrence localities
est.loc = extract(pr,  eval@occ.pts)
#extract model estimated suitability for background
est.bg = extract(pr, eval@bg.pts)
#evaluate predictive ability of model
ev = evaluate(est.loc, est.bg)
#detect possible thresholds 
thr = threshold(ev)
#plot using "equal sensitivity and specificity" criteria
pr_thr = pr>thr$sensitivity
pr_thr_df = as.data.frame(pr_thr, xy=TRUE)

thr_plot = ggplot() +
  geom_raster(data = pr_thr_df, aes(x = x, y = y, fill = layer)) +
  geom_point(data=sp_df, aes(x=longitude, y=latitude), col='red', cex=0.2) +
  scale_fill_manual(values = c('black', 'blue')) +
  coord_quickmap() +
  theme_bw() 

ggsave(thr_plot, 
       filename = paste("figures/", taxon, "_predict.png", sep=""),
       height=7.25, width = 7.25, units='in',
       dpi = 300)



# ### Extend to Future Climate ### #


pr_future26 = predict(preds26, 
                      eval@models[[bestmod]], 
                      type = 'cloglog')
pr_df26 = as.data.frame(pr_future26, xy=T)


#heatmap
max_plot26 = ggplot() +
  geom_raster(data = pr_df26, aes(x = x, y = y, fill = layer)) +
  geom_point(data=sp_df, aes(x=longitude, y=latitude), col='red', cex=0.05) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=viridis::viridis(99),
                       na.value = "black")
ggsave(max_plot26, 
       filename = paste("figures/", taxon, "_maxent26.png", sep=""),
       height=7.25, width = 7.25, units='in',
       dpi = 300)

pr_future26 = predict(preds26, 
                      eval@models[[bestmod]], 
                      type = 'cloglog')
pr_df26 = as.data.frame(pr_future26, xy=T)

#heatmap
max_plot85 = ggplot() +
  geom_raster(data = pr_df85, aes(x = x, y = y, fill = layer)) +
  geom_point(data=sp_df, aes(x=longitude, y=latitude), col='red', cex=0.05) +
  coord_quickmap() +
  theme_bw() + 
  scale_fill_gradientn(colours=viridis::viridis(99),
                       na.value = "black")
ggsave(max_plot85, 
       filename = paste("figures/", taxon, "_maxent85.png", sep=""),
       height=7.25, width = 7.25, units='in',
       dpi = 300)

writeRaster(pr_future26, 
            filename = paste('data/', taxon, '_model26'), 
            overwrite=TRUE)

writeRaster(pr_future85, 
            filename = paste('data/', taxon, '_model85'), 
            overwrite=TRUE)
