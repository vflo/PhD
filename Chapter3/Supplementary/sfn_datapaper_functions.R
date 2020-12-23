library(sapfluxnetr)
library(purrr)
library(tidyverse)
library(viridis)

# Negate in

`%ni%` = Negate(`%in%`)

# Function to get number of trees with measurements per day
# To consider a day with measurements uses a threshold of number
# of timesteps


get_ntrees_day<- function(sfn_data_obj,n_threshold=0){
  sfn_data_obj %>% 
    sfn_metrics(period="1 day",
                .funs = list(~sum(!is.na(.))),
                solar=FALSE,
                interval='general') %>% magrittr::extract2('sapf') %>% 
    dplyr::select(-TIMESTAMP_coll) %>% 
    mutate(n_trees =rowSums(.[-1]>n_threshold),
           si_code=get_si_code(sfn_data_obj)) %>%
    dplyr::select(si_code,TIMESTAMP,n_trees) 
  
}


# TODO:
# average per species or select trees?

sfn_finger<- function(sfn_data_obj,years=c(2011,2012)){
  
  sfn_data_obj %>% 
    get_sapf_data() %>% 
    gather(-TIMESTAMP,key='tree',value='sf')->sfdata
  
  sfn_data_obj %>%  
    get_plant_md() -> plmdata
  
  sfdata %>% 
    full_join(dplyr::select(plmdata,pl_code,pl_species),by=c('tree'='pl_code')) %>% 
    mutate(year=lubridate::year(TIMESTAMP),
           doy = lubridate::yday(TIMESTAMP),
           hour = lubridate::hour(TIMESTAMP)) %>%
    
    ggplot(.,aes(x=hour,y=doy,fill=sf))+
    geom_raster(interpolate=TRUE)+
    # geom_tile(color= "white",size=0.01) + 
    viridis::scale_fill_viridis(name="sf",option ="C")+
    facet_grid(year~pl_species) +
    theme_light() +
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 16),
          strip.background = element_rect(fill = 'white', colour = 'darkgray'),
          strip.text = element_text(size = 16, colour = 'black', face = 'bold'))
  
  
}


# sfn per species

sfn_finger_species<- function(sfn_data_obj,
                              years=1990:2020,
                              species= get_species_md(sfn_data_obj) %>% pull(sp_name) %>% unique()){
  
  sfn_data_obj %>% 
    get_sapf_data() %>% 
    gather(-TIMESTAMP,key='tree',value='sf')->sfdata
  
  sfn_data_obj %>%  
    get_plant_md() -> plmdata
  
  name_legend <- bquote(
    atop(
      Sap~flow~per~phantom(),
      atop(sapwood~area~phantom(), (cm^3~cm^-2~h^-1))
    )
  )
  
  # name_legend <- expression(paste(
  #   "Sap flow per \nsapwood area \n(cm"^{3}*"cm"^{-2}*"h"^{-1}*")"
  # ))
  
  sfdata %>% 
    full_join(dplyr::select(plmdata,pl_code,pl_species),by=c('tree'='pl_code')) %>% 
    group_by(pl_species) %>% 
    mutate(year=lubridate::year(TIMESTAMP),
           doy = lubridate::yday(TIMESTAMP),
           hour = lubridate::hour(TIMESTAMP)) %>% 
    dplyr::filter(year%in%years, pl_species%in%species) %>% 
    dplyr::mutate(
      pl_species = dplyr::if_else(
        pl_species == 'Mortoniodendron anisophyllum', 'Mortoniodendron sp.', pl_species
      )
    ) %>% 
    group_by(pl_species,year,doy,hour) %>% 
    mutate(sf_species=mean(sf,na.rm=TRUE)) %>%
    distinct(pl_species,doy,hour,.keep_all = TRUE) %>%
    
    ggplot(.,aes(x=hour,y=doy,fill=sf_species))+
    geom_raster(interpolate=TRUE)+
    # geom_tile(color= "white",size=0.01) + 
    viridis::scale_fill_viridis(
      # name="Sap flow per\nsapwood area\n[cm³ cm⁻² h⁻¹]",
      name = name_legend,
      option ="C",
      na.value = 'transparent'
    )+
    xlab('Time of day') + ylab('Day of year')+
    scale_x_continuous(labels = c('6h', '12h', '18h'), breaks = c(6,12,18)) +
    facet_grid(year~pl_species)+
    theme_light() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.background = element_rect(fill = 'white', colour = 'darkgray'),
          strip.text.x = element_text(size = 10, colour = 'black', face = 'italic'),
          strip.text.y = element_text(size = 10, colour = 'black'))
  
  
}


# biome plot functions
qc_get_biomes_spdf <- function(merge_deserts = FALSE, parent_logger = 'test') {
  
  # Using calling handlers to logging
  withCallingHandlers({
    
    # STEP 0
    # Argument checks
    # Is merge_deserts logical?
    if (!(is.logical(merge_deserts))) {
      stop('merge_deserts must be logical')
    }
    # Is merge_deserts NA?
    if (is.na(merge_deserts)) {
      stop('merge_deserts must be either TRUE or FALSE')
    }
    
    # STEP 1
    # Create the data frame
    biomes_df <- data.frame(
      mat = c(29.339, 13.971, 15.371, 17.510, 24.131, 27.074, 28.915, 29.201, 29.339, 13.971, -9.706, -7.572,  4.491, 17.510,
              15.371, 13.971, 17.510,  4.491, -7.572, -9.706, -6.687, -0.949,  3.098,  7.147, 10.165, 13.918, 18.626, 18.176,
              17.510, 18.626, 13.918, 10.165,  7.147,  3.098, -0.949,  1.039,  1.998,  2.444,  3.118,  4.446,  7.758, 12.614,
              18.720, 18.637, 18.626, -0.949, -6.687, -4.395, -4.098, -1.592,  0.914,  4.155,  3.118,  2.444,  1.998,  1.039,
              -0.949, 18.720, 12.614,  7.758,  4.446,  3.118,  4.155, 15.716, 20.136, 19.392, 18.720, 18.720, 19.392, 20.136,
              22.278, 23.756, 24.199, 24.714, 25.667, 26.105, 27.414, 27.772, 25.709, 21.736, 18.720, 17.510, 18.176, 18.626,
              18.637, 18.720, 21.736, 25.709, 27.772, 28.418, 28.915, 27.074, 24.131, 17.510, -6.687, -8.896, -9.706, -13.382,
              -15.366, -15.217, -8.373, -4.098, -1.592, -4.098, -4.395, -6.687),
      map = c(21.3,  23.0, 174.6, 535.1, 702.9, 847.9, 992.4, 532.1,  21.3,  23.0,  7.3,  87.2, 314.6, 535.1, 174.6,  23.0,
              535.1, 314.6,  87.2,   7.3, 202.6, 391.7, 529.9, 783.1, 956.9,1116.5,1269.3, 794.3, 535.1,1269.3,1116.5, 956.9,
              783.1, 529.9, 391.7, 514.8, 673.4, 968.5,1630.6,1839.7,2028.0,2224.0,2355.7,1837.6,1269.3, 391.7, 202.6, 922.9,
              1074.1,1405.9,1744.9,2012.3,1630.6, 968.5,673.4, 514.8, 391.7,2355.7,2224.0,2028.0,1839.7,1630.6,2012.3,2930.1,
              3377.7,2917.0,2355.7,2355.7,2917.0,3377.7,3896.5,4343.1,4415.2,4429.8,4279.0,4113.7,3344.4,2790.6,2574.0,2414.3,
              2355.7, 535.1, 794.3,1269.3,1837.6,2355.7,2414.3,2574.0,2790.6,1920.3, 992.4, 847.9, 702.9, 535.1, 202.6,  50.8,
              7.3,  34.8,  98.8, 170.8, 533.0,1074.1,1405.9,1074.1, 922.9, 202.6),
      biome = c(rep('Subtropical desert', 9), rep('Temperate grassland desert', 7), rep('Woodland/Shrubland', 13),
                rep('Temperate forest', 16), rep('Boreal forest', 12), rep('Temperate rain forest', 10),
                rep('Tropical rain forest', 14), rep('Tropical forest savanna', 13), rep('Tundra', 12))
    )
    
    # STEP 2
    # Merge deserts if specified
    if (merge_deserts){
      
      biome <- as.character(biomes_df$biome)
      
      biome[grepl('desert', biome, fixed = TRUE)] <- 'Desert'
      
      biomes_df$biome <- as.factor(biome)
      
    }
    
    # STEP 3
    # Create SpatialPolygonsDataFrame object
    list_pol <- sapply(as.character(unique(biomes_df$biome)),
                       function(id_biome,df)
                         sp::Polygon(cbind(df$map[df$biome == id_biome],
                                           df$mat[df$biome == id_biome])),
                       df=biomes_df, USE.NAMES = TRUE)
    
    sp_biomes <- sp::SpatialPolygons(
      lapply(1:length(list_pol),
             function(i, x) {sp::Polygons(list(x[[i]]),
                                          names(x)[i])},
             x = list_pol)
    )
    
    spdf_biomes <- sp::SpatialPolygonsDataFrame(sp_biomes, data.frame(biome = names(list_pol)),
                                                match.ID = 'biome')
    
    # STEP 4
    # Return SpatialPolygonsDataFrame object
    return(spdf_biomes)
    
    # END FUNCTION
  },
  
  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger, 'qc_get_biomes_spdf', sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger, 'qc_get_biomes_spdf', sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger, 'qc_get_biomes_spdf', sep = '.'))})
  
}

# vis biome base function
vis_biome <- function(merge_deserts = FALSE, parent_logger = 'test') {
  
  # Using calling handlers to logging
  withCallingHandlers({
    
    # STEP 0
    # Argument checks
    # Is merge_deserts logical?
    if (!(is.logical(merge_deserts))) {
      stop('merge_deserts must be logical')
    }
    # Is merge_deserts NA?
    if (is.na(merge_deserts)) {
      stop('merge_deserts must be either TRUE or FALSE')
    }
    
    # STEP 1
    # Get biomes SpatialPointsDataFrame object
    suppressMessages(
      biomes_df <- fortify(qc_get_biomes_spdf(merge_deserts = merge_deserts))
    )
    
    # STEP 2
    # Make and return the plot object
    # 2.1 Make color palette
    if (merge_deserts){
      
      pal <- viridis::viridis(9)[c(2,9,3,4,6,7,8,1)]
      
    } else {
      
      pal <- viridis::viridis(9)[c(2,3,5,4,9,6,7,8,1)]
      
    }
    
    # 2.2 Make the plot object
    plot <- ggplot() +
      geom_polygon(
        data = biomes_df,
        aes(x = long, y = lat, group = id,
            fill = id)
      ) +
      # scale_x_continuous(limits = c(0, 6000)) +
      scale_fill_manual('Biomes', values = pal) +
      xlab('Mean annual precipitation (mm)') +
      ylab('Mean annual temperature (ºC)')
    
    # 2.3 Return the plot object
    return(plot)
    
    # END FUNCTION
  },
  
  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger,
                                                        'vis_biome',
                                                        sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger,
                                                       'vis_biome',
                                                       sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger,
                                                        'vis_biome',
                                                        sep = '.'))})
  
}

# biome and points
vis_location_biome <- function(data, merge_deserts = FALSE,
                               parent_logger = 'test') {
  
  # Using calling handlers to logging
  withCallingHandlers({
    
    # STEP 0
    # Argument checks
    # Is data a data.frame?
    if (!is.data.frame(data)) {
      stop('Provided data object is not a data.frame.',
           ' Please verify if it is the correct object')
    }
    # Does data contains a longitude variable?
    if (is.null(data$si_long)) {
      stop('There is no longitude variable in this dataset. ',
           'Please verify if it is the correct data')
    }
    # Does data contains a latitude variable?
    if (is.null(data$si_lat)) {
      stop('There is no latitude variable in this dataset. ',
           'Please verify if it is the correct data')
    }
    # Is merge_deserts logical?
    if (!(is.logical(merge_deserts))) {
      stop('merge_deserts must be logical')
    }
    # Is merge_deserts NA?
    if (is.na(merge_deserts)) {
      stop('merge_deserts must be either TRUE or FALSE')
    }
    
    # STEP 1
    # Get MAT and MAP if not provided
    if (!all(c('si_mat', 'si_map') %in% names(data))){
      data <- qc_get_biome(data, merge_deserts = merge_deserts)
    }
    
    # STEP 2
    # Make the plot
    # 2.1 Get biome plot
    plot <- vis_biome(merge_deserts = merge_deserts)
    
    # 2.2 Make the plot object
    plot <- plot +
      geom_point(data = data, aes(
        x = si_map, y = si_mat
      ),
      color = 'black', shape = 21, fill = 'white', size = 2, stroke = 0.5) +
      theme_bw() +
      coord_cartesian(xlim = c (0, 4500), ylim = c(-16, 30), expand = FALSE)
    
    # 2.3 Return the plot object
    return(plot)
    
    # END FUNCTION
  },
  
  # handlers
  warning = function(w){logging::logwarn(w$message,
                                         logger = paste(parent_logger,
                                                        'vis_location_biome',
                                                        sep = '.'))},
  error = function(e){logging::logerror(e$message,
                                        logger = paste(parent_logger,
                                                       'vis_location_biome',
                                                       sep = '.'))},
  message = function(m){logging::loginfo(m$message,
                                         logger = paste(parent_logger,
                                                        'vis_location_biome',
                                                        sep = '.'))})
  
}

##########################################################################################
#
# Split violin plot code taken from here:
# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
#
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
##########################################################################################