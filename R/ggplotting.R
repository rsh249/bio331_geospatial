rasterToPlot = function(raster, color = "viridis"){
  names(raster)= 'Suitability'
  df = as.data.frame(raster, xy = TRUE)
  plot = ggplot() + 
    geom_raster(data = df,
                aes(x = x, y = y, fill = Suitability))
  #remove needless text
  plot = plot + theme(legend.position = 'none',  #remove legend
                      axis.title = element_blank(), #remove axis title
                      axis.ticks = element_blank(), #remove axis ticks
                      axis.text = element_blank()) #remove axis text
  #color scale
  plot = plot + scale_fill_viridis_c(option = color)
  return(plot)
}
g_legend <- function(a.gplot){ 
  plot = a.gplot + theme(legend.position='right') +theme(legend.key.height = unit(1, 'in')) + theme(legend.key.width = unit(1, 'in')) + 
    labs(fill = "Location Probabiltiy")
  tmp <- ggplot_gtable(ggplot_build(plot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]]
  
  return(legend)
} 


plot_diva_cclgmbi = rasterToPlot(diva_cclgmbi)
plot_diva_melgmbi = rasterToPlot(diva_melgmbi)
plot_diva_mrlgmbi = rasterToPlot(diva_mrlgmbi)
plot_diva_wc2 = rasterToPlot(diva_wc2)

plot_legend = g_legend(plot_diva_wc2)
grid.newpage()
grid.draw(plot_legend)

plot_trid_cclgmbi = rasterToPlot(trid_cclgmbi)
plot_trid_melgmbi = rasterToPlot(trid_melgmbi)
plot_trid_mrlgmbi = rasterToPlot(trid_mrlgmbi)
plot_trid_wc2 = rasterToPlot(trid_wc2)

plot_combined_cclgmbi = rasterToPlot(combined_cclgmbi)
plot_combined_melgmbi =rasterToPlot(combined_melgmbi)
plot_combined_mrlgmbi =rasterToPlot(combined_mrlgmbi)
plot_combined_wc2 =rasterToPlot(combined_wc2)

all_plots = ggarrange(plot_diva_wc2, plot_diva_cclgmbi, plot_diva_melgmbi, plot_diva_mrlgmbi,
                      plot_trid_wc2, plot_trid_cclgmbi, plot_trid_melgmbi, plot_trid_mrlgmbi,
                      plot_combined_wc2, plot_combined_cclgmbi, plot_combined_melgmbi, plot_combined_mrlgmbi,
                      plot_legend,
                      ncol = 5, nrow = 3, labels = c('A','B','C','D','E','F','G','H','J','K','L','M','N'))


lay <- rbind(c(21,21,21,21,21,13),
             c(14,1,2,3,4,13),
             c(15,5,6,7,8,13),
             c(16,9,10,11,12,13),
             c(NA,17,18,19,20,13))
diva_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "L. Divaricata"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
trid_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "L. Tridentata"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
combined_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "Combined"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
wc2_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "wc2"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
cclgmbi_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "cclgmbi"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
melgmbi_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "melgmbi"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
mrlgmbi_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "mrlgmbi"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
main_title = ggplot() + geom_text(aes(x =0 , y = 0, label = "Divaricata and Tridentata predictions in paleo and modern climates"), size = 10) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())


all_plots = grid.arrange(plot_diva_wc2, plot_diva_cclgmbi, plot_diva_melgmbi, plot_diva_mrlgmbi,
                         plot_trid_wc2, plot_trid_cclgmbi, plot_trid_melgmbi, plot_trid_mrlgmbi,
                         plot_combined_wc2, plot_combined_cclgmbi, plot_combined_melgmbi, plot_combined_mrlgmbi,
                         plot_legend,
                         diva_label, trid_label, combined_label, wc2_label, cclgmbi_label, melgmbi_label, mrlgmbi_label, main_title,
                         layout_matrix = lay)
all_plots

ggsave("all_plots.png", plot = all_plots, device = "png", width = 18, height = 9, unit = "in")

#viridis library


plot_diva_wc2 + scale_fill_viridis_c(option = "plasma")
ggsave("scale_plasma.png", device = "png", width = 9, height = 9, unit = "in")

plot_diva_wc2 + scale_fill_viridis_c(option = "inferno")
ggsave("scale_inferno.png", device = "png", width = 9, height = 9, unit = "in")

plot_diva_wc2 + scale_fill_viridis_c(option = "magma")
ggsave("scale_magma.png", device = "png", width = 9, height = 9, unit = "in")

plot_diva_wc2 + scale_fill_viridis_c(option = "viridis")
ggsave("scale_viridis.png", device = "png", width = 9, height = 9, unit = "in")


plot_diva_wc2 + scale_fill_viridis_c(option = "cividis")
ggsave("scale_cividis.png", device = "png", width = 9, height = 9, unit = "in")

plot_diva_wc2 + scale_fill_gradient2( low = muted("red"), mid = "white",
                                      high = muted("blue"), midpoint = 0, space = "Lab",
                                      na.value = "grey50", guide = "colourbar", aesthetics = "fill")
plot_diva_wc2 +  scale_fill_gradientn(colours = terrain.colors(10), na.value = "navy")
ggsave("scale_terrain3.png", device = "png", width = 9, height = 9, unit = "in")