# libraries
library(sapfluxnetr)
library(tidyverse)
library(furrr)
plan('multicore')
options('future.global.maxsize'=2*1024*1024^2)

getTwoSeasons <- function(input.date,lat){
  library(lubridate)
  numeric.date <- (100*month(input.date)+day(input.date))
  ## input Seasons upper limits in the form MMDD in the "break =" option:
  cuts <- base::cut(numeric.date, breaks = c(0,415,1015,1231))
  # rename the resulting groups (could've been done within cut(...levels=) if "Winter" wasn't double
  if(unique(lat) >= 0){
    levels(cuts) <- c(
      "Winter",
      "Summer",
      paste("Winter",'1',sep = '_')
    )
  }else{
    levels(cuts) <- c(
      "Summer",
      "Winter",
      paste("Summer",'1',sep = '_')
    )  
  }
  return(cuts)
}

#### Plant data preparation

#sapwood sites
folder_sapwood <- "~/sapfluxnet_db/0.1.4/RData/sapwood"

sfn_metadata_sapwood <- read_sfn_metadata(folder = folder_sapwood, .write_cache = TRUE)
sapwood_sites <- sfn_metadata_sapwood$site_md$si_code 

purrr::map(
  as.list(sapwood_sites),
  purrr::possibly(
    function(.x){
      read_sfn_data(.x,
                    folder=folder_sapwood) %>%
        sfn_metrics(period = '1 day', 
                    .funs = sapfluxnetr:::.fixed_metrics_funs(0.95, FALSE)[c(1,3)],
                    # .funs = list(~ mean(., na.rm = TRUE), ~ sd(., na.rm = TRUE), ~ n()),
                    solar = TRUE,
                    interval = 'daylight',
                    # interval = 'general',
                    int_start = 6,
                    int_end = 18
                    ) %>% 
        metrics_tidyfier( metadata = sfn_metadata_sapwood,
                          interval = 'daylight'
                          # interval = 'general'
                          ) ->res

      colnames(res) <- sub("_daylight", "", colnames(res))
      res %>%   
      filter(
          !is.na(ta_mean),
          !is.na(vpd_mean)
        ) %>% 
        mutate(
          date_ymd = lubridate::round_date(TIMESTAMP,unit = 'day'),
          leaf_stability_period = TRUE,
          # site_max_swc = max(swc_shallow_mean, na.rm = TRUE),
          sub_zero = case_when(ta_mean <= 0 ~ "YES",
                               ta_mean > 0 ~ 'NO'),
          season = getTwoSeasons(TIMESTAMP, lat = si_lat),
          season_add = case_when(season == "Winter" ~ 
                                   paste("Winter",year(TIMESTAMP),sep = "_"),
                                 season == "Summer"~ 
                                   paste('Summer',year(TIMESTAMP),sep = "_"),
                                 season == "Winter_1" ~ 
                                   paste("Winter",year(TIMESTAMP)+1,sep = "_"),
                                 season == "Summer_1"~ 
                                   paste('Summer',year(TIMESTAMP)+1,sep = "_"))
        ) -> res
      
      save(res, file=paste0('~/Chapter_2/data/sites/',
                           unique(res$si_code), '.RData'))
      
    },
    otherwise = NULL
  )
)


#plant sites
folder_plant <- "~/sapfluxnet_db/0.1.4/RData/plant"

sfn_metadata_plant <- read_sfn_metadata(folder = folder_plant, .write_cache = TRUE)
plant_sites <- sfn_metadata_plant$site_md$si_code 

plant_sites <- plant_sites[!plant_sites%in%sapwood_sites]

purrr::map(
  as.list(plant_sites),
  purrr::possibly(
    function(.x){
      read_sfn_data(.x,
                    folder=folder_plant) %>%
        sfn_metrics(period = '1 day', 
                    .funs = sapfluxnetr:::.fixed_metrics_funs(0.95, FALSE)[c(1,3)],
                    # .funs = list(~ mean(., na.rm = TRUE), ~ sd(., na.rm = TRUE), ~ n()),
                    solar = TRUE,
                    interval = 'daylight',
                    # interval = 'general',
                    int_start = 6,
                    int_end = 18
        ) %>% 
        metrics_tidyfier( metadata = sfn_metadata_plant,
                          interval = 'daylight'
                          # interval = 'general'
        ) ->res
      
      colnames(res) <- sub("_daylight", "", colnames(res))
      res %>%   
        filter(
          !is.na(ta_mean),
          !is.na(vpd_mean)
        ) %>% 
        mutate(
          date_ymd = lubridate::round_date(TIMESTAMP,unit = 'day'),
          leaf_stability_period = TRUE,
          # site_max_swc = max(swc_shallow_mean, na.rm = TRUE),
          sub_zero = case_when(ta_mean <= 0 ~ "YES",
                               ta_mean > 0 ~ 'NO'),
          season = getTwoSeasons(TIMESTAMP, lat = si_lat),
          season_add = case_when(season == "Winter" ~ 
                                   paste("Winter",year(TIMESTAMP),sep = "_"),
                                 season == "Summer"~ 
                                   paste('Summer',year(TIMESTAMP),sep = "_"),
                                 season == "Winter_1" ~ 
                                   paste("Winter",year(TIMESTAMP)+1,sep = "_"),
                                 season == "Summer_1"~ 
                                   paste('Summer',year(TIMESTAMP)+1,sep = "_"))
        ) -> res
      
      save(res, file=paste0('~/Chapter_2/data/sites/',
                            unique(res$si_code), '.RData'))
      
    },
    otherwise = NULL
  )
)



gc(reset = TRUE)
