#### ASSEMBLAGE OF FUNCTIONS FOR PLOTTING MAPS FOR HYDROMETEOROLOGICAL MODELS

datamapcal <- function(data_stack, model_bio,st_height){
  faa <- raster::predict(data_stack, model_bio)
  faa[is.infinite(faa)]<- NA
  faa[faa>100]<- 100
  faa[faa<0]<- 0
  st_height[st_height <= 0.5]<-NA
  faa <- mask(faa,st_height)
  return(faa)
}
library("maps")
library("ggplot2")
theme_set(theme_bw())
theme_update(panel.grid = element_blank())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "small", returnclass = "sf")
maps_plot <- function(datamap = datamap, title = title){
  datamap <- as(datamap, "SpatialPixelsDataFrame")
  datamap <- as.data.frame(datamap)
  colnames(datamap) <- c("value", "x", "y")
  ggplot(data=world) +
    geom_sf(colour="gray70", fill="gray70",size = 0.1)+
    geom_tile(data=datamap, aes(x=x, y=y, fill=value)) +
    labs(title = bquote(~R[.(title)]^2))+
    scale_fill_gradientn(name=expression(~R^2),
                         na.value = "transparent", colours=rev(brewer.pal(9, 'Spectral')),
                         breaks = c(0,25,50,75,100),
                         labels = c("0%","25%","50%","75%","100%")
                         )+
    theme(legend.position = "none",plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    xlab("Lon") + ylab("Lat") + 
    coord_sf(crs = 4326, xlim=c(-180, 180), ylim=c(-60,85),expand=FALSE)+
    NULL
}

maperrorpred <- function(zeta1, model_bio, datamap, st_height){
  st_height[st_height <= 0.5]<-NA
  df <- stats::predict(model_bio, newdata = zeta1 %>% dplyr::select(-x,-y),
                       na.rm=TRUE, nsim = 200, se.fit = TRUE, alpha = 0.05)
  df1 <-  zeta1 %>% dplyr::select(x,y) %>% cbind(df[['se.fit']] %>% data.frame())
  map <- rasterFromXYZ(df1,crs = crs(datamap))
  map <- extend(map,st_height,progress="text")
  map <- mask(map,st_height)
  map[map>100]<- 100
  return(map)
}


map_error_plot <- function(datamaperror, title){
  datamaperror <- as(datamaperror, "SpatialPixelsDataFrame")
  datamaperror <- as.data.frame(datamaperror)
  colnames(datamaperror) <- c("value", "x", "y")
  ggplot(data=world) +
    geom_sf(colour="gray70", fill="gray70",size = 0.1)+
    geom_tile(data=datamaperror, aes(x=x, y=y, fill=value)) +
    labs(title = bquote("SE "~R[.(title)]^2))+
    scale_fill_gradientn(name=expression(~R^2),
                         na.value = "transparent", colours=c('white',"darkred"),
                         breaks = c(0,10,20,30),
                         labels = c("0%","10%", "20%", "30%")
                         )+
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    xlab("Lon") + ylab("Lat") + 
    coord_sf(crs = 4326, xlim=c(-180, 180), ylim=c(-60,85),expand=FALSE)+
    NULL
}

gplot_data <- function(x, maxpixels = 500000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  dat <- dplyr::as_tibble(data.frame(coords, dat))
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

maps_plot_var <- function(datamap = datamap, title = title){
  datamap <- as(datamap, "SpatialPixelsDataFrame")
  datamap <- as.data.frame(datamap)
  colnames(datamap) <- c("value", "x", "y")
  ggplot(data=world) +
    geom_sf(colour="gray70", fill="gray70",size = 0.1)+
    geom_tile(data=datamap, aes(x=x, y=y, fill=value)) +
    labs(title = bquote(.(title)))+
    scale_fill_gradientn(na.value = "transparent", colours=brewer.pal(9, 'Spectral'))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.justification=c(0,1), 
          legend.position=c(0.03, 0.70),
          legend.background = element_blank(),
          legend.key.size = unit(0.28, "cm"),
          legend.key.width = unit(0.28,"cm"),
          legend.key = element_blank(),
          legend.title = element_blank())+
    xlab("Lon") + ylab("Lat") + 
    coord_sf(crs = 4326, xlim=c(-180, 180), ylim=c(-60,80),expand=FALSE)+
    NULL
}




library(ggbiome)

gd_get_biomes_spdf_var <- function (merge_biomes = FALSE) 
{
  if (!(is.logical(merge_biomes))) {
    stop("merge_biomes must be logical.")
  }
  if (is.na(merge_biomes)) {
    stop("merge_biomes must be either TRUE or FALSE.")
  }
  biomes_df <- data.frame(mat = c(29.339, 13.971, 15.371, 
                                  17.51, 24.131, 27.074, 28.915, 29.201, 29.339, 13.971, 
                                  -9.706, -7.572, 4.491, 17.51, 15.371, 13.971, 17.51, 
                                  4.491, -7.572, -9.706, -6.687, -0.949, 3.098, 7.147, 
                                  10.165, 13.918, 18.626, 18.176, 17.51, 18.626, 13.918, 
                                  10.165, 7.147, 3.098, -0.949, 1.039, 1.998, 2.444, 3.118, 
                                  4.446, 7.758, 12.614, 18.72, 18.637, 18.626, -0.949, 
                                  -6.687, -4.395, -4.098, -1.592, 0.914, 4.155, 3.118, 
                                  2.444, 1.998, 1.039, -0.949, 18.72, 12.614, 7.758, 4.446, 
                                  3.118, 4.155, 15.716, 20.136, 19.392, 18.72, 18.72, 
                                  19.392, 20.136, 22.278, 23.756, 24.199, 24.714, 25.667, 
                                  26.105, 27.414, 27.772, 25.709, 21.736, 18.72, 17.51, 
                                  18.176, 18.626, 18.637, 18.72, 21.736, 25.709, 27.772, 
                                  28.418, 28.915, 27.074, 24.131, 17.51, -6.687, -8.896, 
                                  -9.706, -13.382, -15.366, -15.217, -8.373, -4.098, -1.592, 
                                  -4.098, -4.395, -6.687), 
                          map = c(21.3, 23, 174.6, 535.1,702.9, 847.9, 992.4, 532.1, 21.3, 
                                  23, 7.3, 87.2, 314.6,535.1, 174.6, 23, 535.1, 314.6, 87.2, 
                                  7.3, 202.6, 391.7, 529.9, 783.1, 956.9, 1116.5, 1269.3, 
                                  794.3, 535.1, 1269.3,1116.5, 956.9, 783.1, 529.9, 391.7, 
                                  514.8, 673.4, 968.5, 1630.6, 1839.7, 2028, 2224, 2355.7, 
                                  1837.6, 1269.3, 391.7, 202.6, 922.9, 1074.1, 1405.9, 1744.9, 
                                  2012.3, 1630.6, 968.5, 673.4, 514.8, 391.7, 2355.7, 2224, 
                                  2028,1839.7, 1630.6, 2012.3, 2930.1, 3377.7, 2917, 2355.7, 
                                  2355.7, 2917, 3377.7, 3896.5, 4343.1, 4415.2, 4429.8, 4279, 
                                  4113.7, 3344.4, 2790.6, 2574, 2414.3, 2355.7, 535.1, 794.3, 
                                  1269.3, 1837.6, 2355.7, 2414.3, 2574,2790.6, 1920.3, 
                                  992.4, 847.9, 702.9, 535.1, 202.6, 50.8, 7.3, 34.8, 98.8, 
                                  170.8, 533, 1074.1, 1405.9, 1074.1,922.9, 202.6),
                          biome = c(rep("Subtropical desert", 9), 
                                    rep("Temperate grassland/desert", 7), 
                                    rep("Woodland/shrubland", 13), 
                                    rep("Temperate forest", 16), 
                                    rep("Boreal forest", 12), 
                                    rep("Temperate rain forest", 10), 
                                    rep("Tropical rain forest", 14), 
                                    rep("Tropical seasonal forest/savanna", 13), 
                                    rep("Tundra", 12)))
  if (merge_biomes) {
    biome <- as.character(biomes_df$biome)
    biome = recode(biome,
                   `Woodland/shrubland` = 'WOOD',
                   `Temperate forest` = "TEMP",
                   `Temperate rain forest` = "TEMP",
                   `Boreal forest` = "BOR",
                   `Tundra` = "BOR",
                   `Tropical seasonal forest/savanna` = "TROP",
                   `Tropical rain forest` = 'TROP',
                   `Temperate grassland/desert` = 'DRY',
                   `Subtropical desert` = 'DRY')
    biomes_df$biome <- as.factor(biome)
  }
  list_pol <- sapply(as.character(unique(biomes_df$biome)), 
                     function(id_biome, df) sp::Polygon(cbind(df$map[df$biome == id_biome], 
                                                              df$mat[df$biome == id_biome])), 
                     df = biomes_df, 
                     USE.NAMES = TRUE)
  sp_biomes <- sp::SpatialPolygons(lapply(1:length(list_pol), 
                                          function(i, x) {
                                            sp::Polygons(list(x[[i]]), names(x)[i])
                                          }, x = list_pol))
  spdf_biomes <- sp::SpatialPolygonsDataFrame(sp_biomes, data.frame(biome = names(list_pol)), 
                                              match.ID = "biome")
  return(spdf_biomes)
}

biome_plot <- function (merge_biomes = FALSE) 
{
  if (!(is.logical(merge_biomes))) {
    stop("merge_biomes must be logical")
  }
  if (is.na(merge_biomes)) {
    stop("merge_biomes must be either TRUE or FALSE")
  }
  suppressMessages(biomes_df <- ggplot2::fortify(gd_get_biomes_spdf_var(merge_biomes = merge_biomes)))
  if (merge_biomes) {
    pal <- viridis::viridis(9)[c(2, 9, 4, 7, 8, 6, 1, 3)]
  }
  else {
    pal <- viridis::viridis(9)[c(2, 5, 4, 9, 7, 8, 6, 1, 
                                 3)]
  }

  plot <- ggplot2::ggplot() + 
    ggplot2::geom_polygon(data = biomes_df, 
                          ggplot2::aes_(x = ~long, y = ~lat, group = ~id, fill = ~id)) + 
    # ggplot2::scale_fill_manual("Biomes", values = pal) + 
    scale_fill_manual("Biome",
                      values = c("#377EB8","#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00"))+
    ggplot2::xlab("Mean annual precipitation (mm)") + 
    ggplot2::ylab(expression(paste("Mean annual temperature ", (degree * C))))
  return(plot)
}

