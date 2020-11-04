
library(raster)
library(ggspatial)
library(rgdal)
library(ggrepel)
library(tidyverse)
library(cowplot)

rasterOptions(tmpdir = '~/temp' )


# 1.1. load map -----------------------------------------------------------

crowther <- raster('~/[databases]/Crowther/Crowther_Nature_Files_Revision_01_WGS84_GeoTiff/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif')
LAI <- raster("~/[databases]/LAI/LAI.tif")
st_height <- raster("~/[databases]/Canopy_height/canopy_height.tif")
MAT <- raster("~/[databases]/CHELSA/CHELSA_bio10_1.tif")
MAP <- raster("~/[databases]/CHELSA/CHELSA_bio10_12.tif")
BIO_4 <- raster("~/[databases]/CHELSA/CHELSA_bio10_4.tif")
BIO_15 <- raster("~/[databases]/CHELSA/CHELSA_bio10_15.tif")
clay <- raster("~/[databases]/SoilGrids/clay_ll.tif")
sand <- raster("~/[databases]/SoilGrids/sand_ll.tif")
# bedrock <- raster("~/[databases]/SoilGrids/BDRICM_M_250m_ll.tif")
bedrock <- raster("~/[databases]/SoilGrids/BDTICM_M_250m_ll.tif")
nitrogen <- raster("~/[databases]/SoilGrids/nitrogen_ll.tif")
ph <- raster("~/[databases]/SoilGrids/ph_ll.tif")
# nitrogen[nitrogen>55]<- NA
# nitrogen <- nitrogen%/%10 #Dominand class
# TOTN <- raster("~/[databases]/WISE/TOTN.tif")
# elevation <- raster("~/[databases]/elevation/wc2.1_30s_elev.tif")
# MO <- raster("~/[databases]/SoilGrids/OCDENS_M_sl4_250m_ll.tif")
# LAI_google <- raster("~/[databases]/LAI/LAI_google.tif")
# PET <- raster("~/[databases]/PET/PET_he_annual/pet_he_yr")

st_height <- crop(st_height,LAI,progress="text")
MAT <- crop(MAT,LAI,progress="text")
MAP <- crop(MAP,LAI,progress="text")
BIO_4 <- crop(BIO_4,LAI,progress="text")
BIO_15 <- crop(BIO_15,LAI,progress="text")
clay <- crop(clay,LAI,progress="text")
sand <- crop(sand,LAI,progress="text")
bedrock <- crop(bedrock,LAI,progress="text")
nitrogen <- crop(nitrogen,LAI,progress="text")
ph <- crop(ph,LAI,progress="text")
# MO <- crop(MO,LAI,progress="text")
# PET <- crop(PET,LAI,progress="text")
# PPET <- MAP/PET
# PPET <- crop(PPET,LAI,progress="text")
# LAI_google <- crop(LAI_google,LAI,progress="text")
# TOTN <- crop(TOTN,LAI,progress="text")
# elevation <- crop(elevation,LAI,progress="text")

# clay <- aggregate(clay,fact=48,progress="text")
# sand <- aggregate(sand,fact=48,progress="text")
# nitrogen <- aggregate(nitrogen,fact=48,progress="text")
# ph <- aggregate(ph,fact=48,progress="text")
bedrock <- raster::aggregate(bedrock,fact=48,progress="text")
MAP <- aggregate(MAP,fact=12,progress="text")
MAT <- aggregate(MAT,fact=12,progress="text")
BIO_4 <- aggregate(BIO_4,fact=12,progress="text")
BIO_15 <- aggregate(BIO_15,fact=12,progress="text")
st_height <- aggregate(st_height,fact=12,progress="text")
# MO <- aggregate(MO,fact=48,progress="text")
# PPET <- aggregate(PPET,fact=12,progress="text")
# LAI_google <- aggregate(LAI_google,fact=22.2639,progress="text")
# TOTN <- aggregate(TOTN,fact=12,progress="text")
# elevation <- aggregate(elevation,fact=12,progress="text")

# 
MAT <- resample(MAT,LAI,progress="text")
MAP <- resample(MAP,LAI,progress="text")
BIO_4 <- resample(BIO_4,LAI,progress="text")
BIO_15 <- resample(BIO_15,LAI,progress="text")
clay<- resample(clay,LAI,progress="text")
sand <- resample(sand,LAI,progress="text")
bedrock <- resample(bedrock,LAI,progress="text")
nitrogen <- raster::resample(nitrogen,LAI,progress="text")
ph <- raster::resample(ph,LAI,progress="text")
# MO <- resample(MO,LAI,progress="text")
# PPET <- resample(PPET,LAI,progress="text")
# LAI_google <- resample(LAI_google,LAI,progress="text")
# TOTN <- raster::resample(TOTN,LAI,progress="text")
# elevation <- raster::resample(elevation,LAI,progress="text")
clay <- clay/10
sand <- sand/10
nitrogen <- nitrogen/100
ph <- ph/10


# nitrogen <- raster::as.factor(nitrogen)

writeRaster(MAT,"~/Chapter_2/data/rasters/MAT_resample.tif")
writeRaster(MAP,"~/Chapter_2/data/rasters/MAP_resample.tif")
writeRaster(BIO_4,"~/Chapter_2/data/rasters/BIO_4_resample.tif")
writeRaster(BIO_15,"~/Chapter_2/data/rasters/BIO_15_resample.tif")
writeRaster(LAI,"~/Chapter_2/data/rasters/LAI_resample.tif")
writeRaster(st_height,"~/Chapter_2/data/rasters/st_height_resample.tif")
writeRaster(clay,"~/Chapter_2/data/rasters/clay_resample.tif",overwrite=TRUE)
writeRaster(sand,"~/Chapter_2/data/rasters/sand_resample.tif",overwrite=TRUE)
writeRaster(bedrock,"~/Chapter_2/data/rasters/bedrock_resample.tif",overwrite=TRUE)
writeRaster(MO,"~/Chapter_2/data/rasters/MO_resample.tif")
writeRaster(PPET,"~/Chapter_2/data/rasters/PPET_resample.tif")
writeRaster(LAI_google,"~/Chapter_2/data/rasters/LAI_google_resample.tif")
writeRaster(nitrogen,"~/Chapter_2/data/rasters/nitrogen_resample.tif",overwrite=TRUE)
writeRaster(ph,"~/Chapter_2/data/rasters/ph_resample.tif",overwrite=TRUE)
writeRaster(TOTN,"~/Chapter_2/data/rasters/TOTN_resample.tif")
writeRaster(elevation,"~/Chapter_2/data/rasters/elevation_resample.tif")
# 
# MAP <- raster("~/Chapter_2/data/rasters/MAP_resample.tif")
# MAT <- raster("~/Chapter_2/data/rasters/MAT_resample.tif")
# LAI <- raster("~/Chapter_2/data/rasters/LAI_resample.tif")
# st_height <- raster("~/Chapter_2/data/rasters/st_height_resample.tif")
# clay <- raster("~/Chapter_2/data/rasters/clay_resample.tif")
# sand <- raster("~/Chapter_2/data/rasters/sand_resample.tif")
# bedrock <- raster("~/Chapter_2/data/rasters/bedrock_resample.tif")
# MO <- raster("~/Chapter_2/data/rasters/MO_resample.tif")
# PPET <- raster("~/Chapter_2/data/rasters/PPET_resample.tif")
# # nitrogen <- raster("~/Chapter_2/data/rasters/nitrogen_resample.tif")
# TOTN <- raster("~/Chapter_2/data/rasters/TOTN_resample.tif")
# 
# data <- stack(st_height,MAT,MAP,clay,sand,bedrock,MO,PPET,TOTN)
# 
# names(data) <- c('st_height','MAT','MAP','clay','sand','bedrock','MO','PPET','TOTN')

# 1.2. Mask preparation -------------------------------------------

# crowther_mask <- crowther
# crowther_mask[crowther_mask < 400] <- NA



# 1.3. Reduce raster resolution and extension -------------------------------------------

# crowther_mask_ext <- crop(crowther_mask, clay)
# bios_ext <-  crop(bios, clay)
# crowther <- crop(crowther, clay)
# 
# clay_resam <- resample(clay, crowther_mask_ext)
# bios_ext_resam <- resample(bios_ext, crowther_mask_ext)
# 
# data_maps <- stack(clay_resam, bios_ext_resam)
# 
# data_maps <- aggregate(data_maps, fact = 0.1/res(data_maps)) # aggregate output
# crowther_mask_ext <- aggregate(crowther_mask_ext, fact = 0.1/res(crowther_mask_ext)) # aggregate output
# 
# crowther[crowther > 0] <- 0 #Base map
# crowther <- aggregate(crowther, fact = 0.1/res(crowther)) # aggregate output
# save(crowther, file = "data/daily_0.1/crowther_base_map.Rdata")
# writeRaster(crowther, filename = "data/daily_0.1/crowther_base_map.tiff", overwrite = TRUE)
# 
# # 1.4.Forest mask aplication ------------
# 
# data_maps_masked <- mask(data_maps, crowther_mask_ext)
# names(data_maps_masked)<-c('clay','BIO_6','MAP')
# save(data_maps_masked,file = "data/daily_0.1/masked_map_0.1.Rdata")
# writeRaster(data_maps_masked,filename = "data/daily_0.1/masked_map_0.1.grd",bandorder='BIL', overwrite = TRUE)
# 


# 2.0 GLOBAL PREDICT ---------------

load("data/daily_0.1/masked_map_0.1.Rdata")
load("data/daily_0.1/crowther_base_map.Rdata")
map01 <- raster::predict(data_maps_masked, mod1, re.form=~0, 
                         const = data.frame(pl_dbh = 20), inf.rm = TRUE)
map01[is.infinite(map01)]<- NA
map01[map01>1]<- 1




# 3.1 Predict plot -----------------

gplot_data <- function(x, maxpixels = 100000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

gplot_wrld_r <- gplot_data(crowther)

g <- ggplot() +
  geom_tile(data = dplyr::filter(gplot_wrld_r, !is.na(value)), 
            aes(x = x, y = y), fill = "grey70") +
  # geom_sf(data=crowther,mapping = aes(values),color = "grey50") +
  layer_spatial(map01) +
  scale_fill_gradientn(na.value = "transparent", colours=c('white',"darkblue"))+
  coord_sf()

g


# 3.2 Error plot -----------------

library('bootpredictlme4')
zeta <- as(data_maps_masked, "SpatialPixelsDataFrame") %>% data.frame()
zeta1 <- zeta[complete.cases(zeta),]
df <- stats::predict(mod1, newdata = zeta1 %>% dplyr::select(-x,-y) %>% mutate(pl_dbh = 20),
                     na.rm=TRUE, nsim = 200, re.form=~0,se.fit = TRUE, alpha = 0.05)

df1 <-  zeta1 %>% dplyr::select(x,y) %>% cbind(df[['se.fit']] %>% data.frame())
map <- rasterFromXYZ(df1,crs = crs(map01))


g <- ggplot() +
  geom_tile(data = dplyr::filter(gplot_wrld_r, !is.na(value)), 
            aes(x = x, y = y), fill = "grey70") +
  layer_spatial(map) +
  scale_fill_gradientn(na.value = "transparent", colours=c('white',"darkred"))+
  coord_sf()
