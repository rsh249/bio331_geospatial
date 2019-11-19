library(spocc)
library(dplyr)
library(RSpatial)
library(ggplot2)
library(raster)
library(ENMeval)


ext = extent(-130, -40, -55, 55)
varindex=c(1,2,3,12)

clim=raster::stack(list.files('/usr/share/data/wc1.4/wc2-5/', pattern='bil', full.names = T))
cclgmbi=raster::stack(list.files('/usr/share/data/paleoclim/cclgmbi', full.names = T))
mrlgmbi=raster::stack(list.files('/usr/share/data/paleoclim/mrlgmbi', full.names = T))
melgmbi=raster::stack(list.files('/usr/share/data/paleoclim/melgmbi', full.names = T))

clim = clim[[varindex]]
cclgmbi = cclgmbi[[varindex]]
mrlgmbi = mrlgmbi[[varindex]]
melgmbi = melgmbi[[varindex]]

clim = raster::crop(clim, ext)
cclgmbi = raster::crop(cclgmbi, ext)
mrlgmbi = raster::crop(mrlgmbi, ext)
melgmbi = raster::crop(melgmbi, ext)

names(cclgmbi) = names(clim)
names(melgmbi) = names(clim)
names(mrlgmbi) = names(clim)



diva_df = occ(query = "Larrea divaricata", 
              from = "gbif",
              limit = 1377,
              has_coords = TRUE)
divaoccdat = occ2df(diva_df)
divaoccdat = filter(divaoccdat, latitude < 15) %>% filter (latitude > -60) %>%
  filter(longitude > -90) %>% filter(longitude < -35 )
divaloc=divaoccdat[,c('longitude', 'latitude')]
occ2thin = poThin(
  df = divaloc,
  spacing = 25,
  dimension = nrow(divaloc),
  lon = 'longitude',
  lat = 'latitude'
)
divaloc = divaloc[-occ2thin,]




evalRandomSample = function(tbl1, predictors){
  require(dplyr)
  require(maxnet)
  require(ENMeval)
  plotSize = floor(0.7*(nrow(tbl1)))
  #plotPercentage = runif(1, min = 0.2, max = 0.95)
  #plotSize = floor( (plotPercentage)*(nrow(tbl1)))
  s1 = sample_n(tbl1, size=plotSize)
  
  s1_eval = ENMevaluate(occ=as.data.frame(s1), env = predictors, method='block', parallel=F, numCores=8,
                        fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
  #tempBest=s1_eval@models[[which(s1_eval@results$AICc == min(s1_eval@results$AICc))]]
  
  best=which(s1_eval@results$AICc == min(s1_eval@results$AICc))
  
  #train maxent on s1
  # For picking model parameters on the complete set
  best_param = s1_eval@results[best, 1]
  best_arr = strsplit(as.character(best_param), "_")
  
  rm = best_arr[[1]][length(best_arr[[1]])]
  
  fc1 = best_arr[[1]][1:(length(best_arr[[1]]) - 1)]
  
  maxmatr = rbind(s1_eval@occ.pts, s1_eval@bg.pts)
  pres = c(rep(1, nrow(s1_eval@occ.pts)), rep(0, nrow(s1_eval@bg.pts)))
  
  maxmatr = maxmatr[,c("LON", 'LAT')]

  
  maxextr = raster::extract(predictors, maxmatr)
  
  #maxmatr = cbind(maxmatr, pres)
  
  best_mod = maxnet(
    p = pres,
    data = as.data.frame(maxextr),
    maxnet.formula(
      p = pres,
      data = as.data.frame(maxextr),
      classes = stringr::str_to_lower(fc1)
    ),
    regmult = as.numeric(rm)
  )
  return(best_mod)
  
}

plotRandomSample =function(bestmod, predictors, varname){
  require(raster)
  s1.predict1 = predict(predictors, bestmod, type = "cloglog")
  names(s1.predict1) = paste(varname, " Layer ")
  return(s1.predict1)
}

erand = evalRandomSample(divaloc, clim)

p1 = plotRandomSample(erand, clim, 'test')

p1_df = as.data.frame(p1, xy = TRUE)
head(p1_df)

## ggplot 

ggplot() +
  geom_raster(data = p1_df, aes(x = x, y = y, fill=test..Layer)) +
  coord_quickmap()



#must declare divaloc
#set up paleoclimate layers


sub_paleopredict <- function(x){
  tempBest = evalRandomSample(divaloc, predictors=clim)
  diva_cclgmbi = plotRandomSample(tempBest, predictors = cclgmbi, "cclgmbi")
  diva_melgmbi = plotRandomSample(tempBest, predictors = melgmbi, "melgmbi")
  diva_mrlgmbi = plotRandomSample(tempBest, predictors = mrlgmbi, "mrlgmbi")
  ret = list(x, diva_cclgmbi, diva_melgmbi, diva_mrlgmbi)
  return(ret)
}

sub_paleopredict <- function(x){
  tempBest = evalRandomSample(divaloc, predictors=clim)
  sp_cclgmbi = plotRandomSample(tempBest, predictors = cclgmbi, "cclgmbi")
  sp_melgmbi = plotRandomSample(tempBest, predictors = melgmbi, "melgmbi")
  sp_mrlgmbi = plotRandomSample(tempBest, predictors = mrlgmbi, "mrlgmbi")
  sp_wc2 = plotRandomSample(tempBest, predictors = clim, "wc2")
  ret = list(x, stack(sp_cclgmbi, sp_melgmbi, sp_mrlgmbi, sp_wc2))
  return(ret)
}


test = sub_paleopredict(1)

t_df = as.data.frame(test[[2]], xy = TRUE)
head(t_df)

## ggplot 

ggplot() +
  geom_raster(data = t_df, aes(x = x, y = y, fill=cclgmbi..Layer)) +
  coord_quickmap()

ggplot() +
  geom_raster(data = t_df, aes(x = x, y = y, fill=melgmbi..Layer)) +
  coord_quickmap()

ggplot() +
  geom_raster(data = t_df, aes(x = x, y = y, fill=mrlgmbi..Layer)) +
  coord_quickmap()

ggplot() +
  geom_raster(data = t_df, aes(x = x, y = y, fill=wc2..Layer)) +
  coord_quickmap()
