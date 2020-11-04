# library(sjstats)
library(tidyverse)
# library(randomForest)
# library(nnet)
# require(caret)
library(taxonlookup)
library(broom)
library(magrittr)
library(rjson)
library(sp)
# library(GSIF)
# library(AICcmodavg)
# library(lmerTest)
library(lme4)
# library(medfate)
source("~/Chapter_3_mod/soil_functions.R")
# source('sfn_NN.R') # for the analysis
library(purrr)
library(future)
library(furrr)
plan('multisession')
options('future.global.maxsize'=2*1024*1024^2)
# options(scipen = 999)

num_NA <- function(x, ...) {
  if (length(x) == 0)
    return(NA)
  return(x)
}


## 0.1 Plant names to be modelized ---------------

# path <- "~/Chapter_3/data/site_species/"
path <- "~/Chapter_2/data/site_daylight/"
site_names <- list.files(path = path)

# source('~/Chapter_2/Individual_plant_analyses/treatment_list_filter.R')
source('~/Chapter_3_mod/treatment_filter_only_control.R')

load("~/Chapter_2/data/swc_ERA5_land.RData") #nearest cell extract ERA5 land 9x9km 1980 to present

ERA5_land_data <- swc_ERA5_land %>% 
  rename(ERA5_swc = swc_ERA5_land,
         TIMESTAMP = TIMESTAMP_daylight) %>% 
  mutate(TIMESTAMP = lubridate::as_date(TIMESTAMP))
# ERA5_land_data %>%
#   group_by(si_code) %>% 
#   mutate(lead_ERA5_swc = lead(ERA5_swc),
#          n_ERA5 = n(),
#          D_ind = log(lead_ERA5_swc+1/ERA5_swc+1),
#          D_sum = sum(D_ind,na.rm = TRUE),
#          D = 1/(n_ERA5-1)*D_sum) -> ERA5_land_data
# D_data <- ERA5_land_data %>% dplyr::select(si_code,D) %>% unique()

dbh_sapw_mod <- readRDS(file="~/Chapter_3_mod/dbh_sapw_area_model.rds")
# load(file = '~/Chapter_3/filter_sites_sapw_units_wrong.Rdata')
# soilgrids.r <- REST.SoilGrids(c('SNDPPT', 'SLTPPT', 'CLYPPT', 'ORCDRC', 'BLD', 'CEC', 'PHIHOX','CRFVOL'))
# soilgrids.r <- REST.SoilGrids(c('bdod'))
## import PET
PET <- read_csv('~/Chapter_2/data/PET.csv') %>% 
  dplyr::rename(si_lat = lati,si_long = long)
BIOS <- read_csv("~/Chapter_2/data/CHELSA_BIOS.csv") %>% mutate(si_long = lon, si_lat = lat)

load("~/Chapter_2/data/sand.RData")
load("~/Chapter_2/data/clay.RData")
load("~/Chapter_2/data/elevation.RData")

load("~/Chapter_2/data/sw_ERA5_18.RData")
load("~/Chapter_2/data/sw_ERA5_6.RData")


sw <- sw_ERA5_18 %>% 
  left_join(sw_ERA5_6) %>% 
  mutate(sw = sw_ERA5_18 - sw_ERA5_6,
         si_code = as.character(si_code)) %>% 
  rename("TIMESTAMP" = "TIMESTAMP_daylight")%>% 
  mutate(TIMESTAMP = lubridate::as_date(TIMESTAMP))

rm(sw_ERA5_18)
rm(sw_ERA5_6)

# site_names <- site_names[!site_names %in%
#                              list.files("data/models/complete_bin_gam_models/")] %>%
#   sort(decreasing =TRUE)
# 
# site_names <- site_names[-which(site_names == "ISR_YAT_YAT.RData")]

site_names <- sample(site_names)


## 0.2 MODELIZATION ----------------------

# furrr::future_map(site_names,.progress=TRUE,.f=function(.x){
furrr::future_map(site_names,.progress=TRUE,.f=function(.x){
# purrr::map(site_names,function(.x){
  ## 1.0 Data load -----------------------
  
  load(paste0(path,.x))
  
  print(.x)
  
  ## 1.1 Conditionals to allow model calculation -----
  
  x2 <- x2 %>% mutate(st_treatment = case_when(is.na(st_treatment) ~ "None",
                                               !is.na(st_treatment) ~ st_treatment))
  if(any(colnames(x2) == "vpd_mean")
  ){ 
    
    if(!any(colnames(x2) == "ppfd_in_mean")){x2 <- x2 %>% mutate(ppfd_in_mean = c(NA))}
    if(!any(colnames(x2) == "sw_in_mean")){x2 <- x2 %>% mutate(sw_in_mean = c(NA))}
    if(!any(colnames(x2) == "swc_shallow_mean")){x2 <- x2 %>% mutate(swc_shallow_mean = c(NA))}
    
    metadata <- x2 %>%
      mutate(pl_species = as.character(pl_species),
             TIMESTAMP = lubridate::ymd(TIMESTAMP))
    
    metadata %>%
      pull(pl_species) %>%
      unique() %>%
      lookup_table(missing_action = 'NA', by_species = TRUE) %>%
      rownames_to_column('pl_species') %>%
      left_join(metadata, ., by = 'pl_species') -> x2
    

    x2 <- x2 %>% cbind(elev = elevation %>% 
                         filter(si_code == unique(x2$si_code)) %>% 
                         pull(elevation))
    
    ## 2.1 Data preparation ---------------
    
    
    x2 %>% left_join(ERA5_land_data, by = c("TIMESTAMP", "si_code")) %>%
      mutate(swc = ERA5_swc,
             # swc = swc_shallow_mean %>% as.double(),
             # swc = if_else(!is.na(swc), swc_shallow_mean %>% as.double(),
             #               ERA5_swc),
             max_swc = max(swc, na.rm = TRUE),
             min_swc = min(swc, na.rm = TRUE),
             max_min_range = max_swc - min_swc,
             rew = (swc-min_swc)/(max_swc - min_swc)
             ) %>% 
      filter(!is.na(vpd_mean)) %>% 
      left_join(sw,by = c("TIMESTAMP", "si_code")) %>% 
      mutate(TIMESTAMP = lubridate::ymd(TIMESTAMP))%>% 
      mutate(
        ba = ((pl_dbh/2)^2)*pi,
        log_ba = log(ba),
        log_pl_sapw_area = predict.lm(object = dbh_sapw_mod, tibble(log_ba,group)),
        pl_sapw_area = case_when(is.na(pl_sapw_area) ~ exp(log_pl_sapw_area),
                                 !is.na(pl_sapw_area) ~ pl_sapw_area),
        sapflow_mean = case_when(pl_sap_units == "cm3 cm-2 h-1" ~ sapflow_mean,# cm3 cm-2 h-1
                                 pl_sap_units == "cm3 h-1" ~ sapflow_mean/pl_sapw_area),#from cm3 h-1 to cm3 cm2 h-1
        si_elev = coalesce(si_elev,elev),
        ppfd_ERA = LakeMetabolizer::sw.to.par.base(sw),
        ppfd_sw_in = LakeMetabolizer::sw.to.par.base(sw_in_mean),
        ppfd = case_when(
          # We priorice ppfd calculated with on-site sw, and secondly ERA5 sw 
          is.na(sw_in_mean)~ppfd_ERA, 
          !is.na(sw_in_mean)~as.double(ppfd_sw_in),
          !is.na(ppfd_in_mean)~as.double(ppfd_in_mean)
          # !is.na(ppfd_in_mean)& is.na(sw_in_mean)~as.double(ppfd_sw_in)
          # is.na(ppfd_in_mean)&!is.na(sw_in_mean)~as.double(ppfd_sw_in)
          )
      ) -> x2
    
    
    x2%>% 
      filter(st_treatment %in% treatments,
             # !is.na(swc),
             !is.na(vpd_mean),
             !is.na(sapflow_mean),
             !is.na(ta_mean),
             vpd_mean>=0.3,
             swc > 0
      )%>%
      mutate(diff = swc / lag(swc, default = first(swc), order_by = TIMESTAMP)) %>% #slope to obtain swc recharge
      filter(#lag(diff,n=1)<=1,
             #lead(diff,n=1)<=1,
             diff<=1
      ) %>% #filter no rain days + - 1 day
      dplyr::mutate(
        sapflow_mean = case_when(pl_sens_meth == "HD" & !pl_sens_calib ~ 
                                   sapflow_mean+0.40488*sapflow_mean,# calibration correction from flo et al. 2019
                                 pl_sens_meth == "HD" & is.na(pl_sens_calib) ~ 
                                   sapflow_mean+0.40488*sapflow_mean,# calibration correction from flo et al. 2019
                                 pl_sens_meth == "HD" & pl_sens_calib ~ sapflow_mean,
                                 pl_sens_meth != "HD" ~ sapflow_mean),
        cond_coeff = 1.158e2 + 4.236e-1 * ta_mean, #[kPa m3 kg-1]
        conductance = cond_coeff*(sapflow_mean/(60*60*10^3))/vpd_mean,#cm3 cm-2 h-1 to m3 vapor cm-2 s-1
        conductance = case_when(!is.na(si_elev) ~ conductance*44.6*(273/(273+ta_mean))*
                                  (101.325*exp(-0.00012*si_elev)/101.325),
                                is.na(si_elev) ~ conductance*44.6*(273/(273+ta_mean))),#from m3 cm-2 s-1 to mol*cm-2 s-1
        G_sw = conductance*10^4,#from mol cm-2 s-1 to mol m-2 s-1
        G_max = quantile(conductance, probs = 0.95,na.rm = TRUE),
        G_min = min(conductance,na.rm = TRUE),
        G_range = G_max - G_min,
        hour = lubridate::hour(TIMESTAMP),
        DOY = lubridate::round_date(TIMESTAMP,unit = 'day') %>% as.numeric(),
        si_lat = mean(si_lat, na.rm =TRUE),
        si_long = mean(si_long, na.rm = TRUE),
        st_sand_perc = mean(st_sand_perc, na.rm = TRUE),
        st_clay_perc = mean(st_clay_perc, na.rm = TRUE),
        st_silt_perc = mean(st_silt_perc, na.rm = TRUE),
        pl_height = coalesce(pl_height, st_height)
      ) %>%
      arrange(TIMESTAMP) -> faa

    faa %>%
      mutate(max_min_range = max_swc - min_swc,
             rew = (swc-min_swc)/(max_swc - min_swc)) -> faa
    
  }else{faa <- tibble(swc = NA)}
  
  if (!is.null(faa) & any(!is.na(faa$swc))){
    
    # 2.2.1 Soil data import --------------
    
    # pnts <- data.frame(lon=unique(faa$si_long),
    #                    lat=unique(faa$si_lat),
    #                    si_code = unique(faa$si_code))%>%
    #   mutate(lon = ifelse(si_code == "CHE_DAV_SEE",9.855,lon),
    #          lat = ifelse(si_code == "CHE_DAV_SEE",46.815,lat),
    #          lon = ifelse(si_code == "GBR_DEV_CON",-3.73,lon),
    #          lat = ifelse(si_code == "GBR_DEV_CON",56.02,lat),
    #          lon = ifelse(si_code == "GBR_DEV_DRO",-3.73,lon),
    #          lat = ifelse(si_code == "GBR_DEV_DRO",56.02,lat),
    #          lon = ifelse(si_code == "BRA_CAX_CON",-51.7,lon),
    #          lat = ifelse(si_code == "BRA_CAX_CON",-1.9,lat)) %>%
    #   dplyr::select(-si_code)
    # coordinates(pnts) <- ~lon+lat
    # proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")
    # # pnts
    # pr <- httr::GET(paste0('https://rest.soilgrids.org/soilgrids/v2.0/properties/query?lon=',
    #                        pnts@coords[1,1],
    #                        '&lat=',
    #                        pnts@coords[1,2],
    #                        '&property=cec&property=cfvo&property=clay&property=nitrogen',
    #                        '&property=ocd&property=ocs&property=phh2o&property=sand',
    #                        '&property=silt&property=soc&property=bdod',
    #                        '&depth=15-30cm&value=mean&value=uncertainty'))
    # prr <- httr::content(pr)$properties$layers
    # 
    # ov <- tibble(bdod = num_NA(prr[1][[1]]$depths[[1]]$values$mean/prr[1][[1]]$unit_measure$d_factor),
    #              cec = num_NA(prr[2][[1]]$depths[[1]]$values$mean/prr[2][[1]]$unit_measure$d_factor),
    #              cfvo = num_NA(prr[3][[1]]$depths[[1]]$values$mean/prr[3][[1]]$unit_measure$d_factor),
    #              clay = num_NA(prr[4][[1]]$depths[[1]]$values$mean/prr[4][[1]]$unit_measure$d_factor),
    #              nitrogen = num_NA(prr[5][[1]]$depths[[1]]$values$mean/prr[5][[1]]$unit_measure$d_factor),
    #              ocd = num_NA(prr[6][[1]]$depths[[1]]$values$mean/prr[6][[1]]$unit_measure$d_factor),
    #              phh2o = num_NA(prr[7][[1]]$depths[[1]]$values$mean/prr[7][[1]]$unit_measure$d_factor),
    #              sand = num_NA(prr[8][[1]]$depths[[1]]$values$mean/prr[8][[1]]$unit_measure$d_factor),
    #              silt = num_NA(prr[9][[1]]$depths[[1]]$values$mean/prr[9][[1]]$unit_measure$d_factor),
    #              soc = num_NA(prr[10][[1]]$depths[[1]]$values$mean/prr[10][[1]]$unit_measure$d_factor))
    # 
    # ov <- over(soilgrids.r, pnts)
    # ov <- ov %>% dplyr::select(BLD = BLDFIE.M.sl4,
    #                            CEC = CECSOL.M.sl4,
    #                            CLYPPT = CLYPPT.M.sl4,
    #                            ORCDRC = ORCDRC.M.sl4,
    #                            PHIHOX = PHIHOX.M.sl4,
    #                            SLTPPT = SLTPPT.M.sl4,
    #                            SNDPPT = SNDPPT.M.sl4,
    #                            CRFVOL = CRFVOL.M.sl4)
    # ov <- ov[1,]
    
    # faa <- faa %>% cbind(ov %>% data.frame(), row.names = NULL)
    faa <- faa %>% cbind(sand = sand %>% 
                           filter(si_code == unique(faa$si_code)) %>% 
                           pull(sand))
    faa <- faa %>% cbind(clay = clay %>% 
                           filter(si_code == unique(faa$si_code)) %>% 
                           pull(clay))
    
    ## 2.2.2 Soil WWP and FC calculation ---------
    
    swp <- NA
    swp.sat <- NA
    if(any(is.na(unique(faa$st_sand_perc)),is.na(unique(faa$st_clay_perc)))){
      swp <- clapp.and.hornberger(faa$swc,
                                  percent.clay=unique(faa$clay),
                                  percent.sand=unique(faa$sand))[["suction"]] #mm
      swp.sat <- clapp.and.hornberger(faa$swc,
                                      percent.clay=unique(faa$clay),
                                      percent.sand=unique(faa$sand))[["suction.sat"]] #mm
      
      swp <- -(swp/101972) #MPa
      swp.sat <- -(swp.sat/101972) #MPa
    }
    if(all(!is.na(unique(faa$st_sand_perc)),!is.na(unique(faa$st_clay_perc)))){
      swp <- clapp.and.hornberger(faa$swc,
                                  percent.clay=unique(faa$st_clay_perc),
                                  percent.sand=unique(faa$st_sand_perc))[["suction"]] #mm
      swp.sat <- clapp.and.hornberger(faa$swc,
                                      percent.clay=unique(faa$st_clay_perc),
                                      percent.sand=unique(faa$st_sand_perc))[["suction.sat"]] #mm
      swp <- -(swp/101972) #MPa
      swp.sat <- -(swp.sat/101972) #MPa
    }
    
    faa <- faa %>% cbind(data.frame(swp = swp,swp.sat = swp.sat), row.names = NULL)
    
  
    faa%>%
      filter(!is.na(G_sw),
             !is.na(swc),
             !is.na(vpd_mean),
             !is.na(ppfd),
             G_sw > 0,
             swc>0,
             ppfd>0,
             swc <= 0.3
      )-> faa
    
    if (!is.null(faa) & any(!is.na(faa$swc))){
    
    faa %>%
      left_join(PET %>% dplyr::select(-c(si_lat,si_long)),by=c('si_code'))%>%
      left_join(BIOS %>% dplyr::select(-c(si_lat,si_long)),by=c('si_code'))%>%
      group_by(pl_code) %>% 
      mutate(pl_height = coalesce(pl_height, st_height),
             n = n()) %>%
      ungroup() %>% 
      summarise(si_code = unique(si_code),
                pl_dbh = weighted.mean(pl_dbh,n),
                pl_height = weighted.mean(pl_height,n),
                si_long = unique(si_long),
                si_lat = unique(si_lat),
                si_elev = unique(si_elev),
                MAP = unique(MAP),
                MAT = unique(MAT)/10,
                PET = unique(PET),
                PPET = MAP/PET,
                si_biome = unique(si_biome),
                min_timestamp = min(TIMESTAMP,na.rm=TRUE),
                max_timestamp = max(TIMESTAMP,na.rm=TRUE),
                # clay = unique(clay),
                # bdod = unique(bdod),  
                # cec = unique(cec), 
                # cfvo = unique(cfvo),
                # nitrogen = unique(nitrogen),  
                # ocd = unique(ocd), 
                # phh2o = unique(phh2o), 
                # silt  = unique(silt),
                # sand  = unique(sand),
                # soc  = unique(soc),
                pl_sap_units = unique(pl_sap_units),
                swc_shallow = ifelse(all(is.na(swc_shallow_mean)),
                                     NA_integer_,
                                     mean(swc_shallow_mean,na.rm = TRUE)
                                     ),
                sw_in_mean = ifelse(all(is.na(sw_in_mean)),
                                   NA_integer_,
                                   mean(sw_in_mean,na.rm = TRUE)
                                   ),
                ppfd_in_mean = ifelse(all(is.na(ppfd_in_mean)),
                                   NA_integer_,
                                   mean(ppfd_in_mean,na.rm = TRUE)
                                   ),
                ppfd_sw_in = ifelse(all(is.na(ppfd_sw_in)),
                                      NA_integer_,
                                      mean(ppfd_sw_in,na.rm = TRUE)
                                    ),
                ppfd = ifelse(all(is.na(ppfd)),
                                    NA_integer_,
                                    mean(ppfd,na.rm = TRUE)
                              )
      )-> env_data
    
    save(faa, file = paste0('./data/models_raw/complete_G_log/',.x))
    
    ## 3.0 GAM calculation -----------------  
    
    faa %>% 
      group_by(pl_code) %>% 
      # mutate(
      #        quan_filter = 'NO',
      #        quan_filter = 
      #          case_when(
      #            G_sw<quantile(G_sw,0.99,na.rm = TRUE) ~ 'YES',
      #            G_sw>quantile(G_sw,0.01,na.rm = TRUE) ~ 'YES')) %>% 
      mutate(
        log_swc = log(swc),
        min_swc = min(log_swc,na.rm = TRUE),
        max_swc = max(log_swc,na.rm = TRUE)) %>% 
      # filter(quan_filter == "YES") %>%
      # filter(log_swc>=-3) %>%
      ungroup() ->faa2
    
    if(nrow(faa2)>1){
      
      faa2%>%
        mutate(#ppfd = ppfd/1000,
               # vpd_cut=cut(vpd_mean, 10),
               vpd_cut=cut(vpd_mean, seq(0.3, 9.9, by=0.2)),
               # swp_cut=cut(log_swc, c(-12,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,100)),
               swc_cut=cut(swc, 5),
               ppfd_cut=cut(ppfd, seq(0, 1750, by=250)),
               max_vpd = max(vpd_mean,na.rm=TRUE),
               min_vpd = min(vpd_mean,na.rm=TRUE),
               max_swc = max(swc,na.rm=TRUE),
               min_swc = min(swc,na.rm=TRUE),
               max_ppfd = max(ppfd,na.rm=TRUE),
               min_ppfd = min(ppfd,na.rm=TRUE),
               range_vpd = max_vpd - min_vpd,
               range_swc = max_swc - min_swc,
               range_ppfd = max_ppfd - min_ppfd,
               # vpd_cut_mean = sapply(str_extract_all(vpd_cut,"-?[0-9]+(\\.[0-9]+)?"), 
               #                       function(x) mean(as.numeric(x))),
               # swc_cut_mean = sapply(str_extract_all(swc_cut,"-?[0-9]+(\\.[0-9]+)?"), 
               #                       function(x) mean(as.numeric(x))),
               # ppfd_cut_mean = sapply(str_extract_all(ppfd_cut,"-?[0-9]+(\\.[0-9]+)?"), 
               #                       function(x) mean(as.numeric(x))),
               n_day_plant = n_distinct(vpd_mean)  
                 )%>%
        filter(n_day_plant>=10,
               range_swc>0.05|range_vpd>0.5)->faa2


      if (all(n_distinct(faa2$vpd_mean) >= 9,
          n_distinct(faa2$swc) >= 9,
          n_distinct(faa2$ppfd) >= 9,
          max(faa2$max_ppfd) >= 400)){
        set.seed(7)
        
        faa2%>%
          group_by(pl_code,vpd_cut,swc_cut,ppfd_cut) %>% 
          summarise(G_sw = mean(G_sw,na.rm = TRUE),
                    sapflow_mean = mean(sapflow_mean,na.rm = TRUE),
                    # G_sw = quantile(G_sw,0.95,na.rm = TRUE),
                    n_dist_G_cut = n_distinct(G_sw),
                    log_vpd_mean = log(mean(vpd_mean)),
                    log_swc = log(mean(swc)),
                    log_ppfd = log(mean(ppfd)/1000),
                    pl_species = unique(pl_species),
                    n_day_plant = unique(n_day_plant)) %>% 
          ungroup() %>% 
          mutate(min_day_plant = min(n_day_plant))->faa2
        
        
        test_1 <- unique(faa2$pl_species) %>% length()>1
        test_2 <- faa2 %>% 
          group_by(pl_species) %>% 
          summarise(n_dis = n_distinct(pl_code)) %>%
          mutate(mode_plants = median(n_dis),
            logi=any(n_dis > 1)) %>%
          dplyr::select(logi) %>% 
          unique() %>% 
          as.logical()
        

        if(test_1 & test_2){

          type = 1

            model_G_sw_log <- try(lmer(G_sw~log_vpd_mean+log_swc+log_ppfd+
                                         (-1+log_vpd_mean+log_swc+log_ppfd|pl_species)+
                                         (1|pl_species:pl_code),REML = TRUE,
                                       data = faa2))
            model_G_sw_log_vpd <- try(lmer(G_sw~log_vpd_mean+
                                             (log_vpd_mean|pl_species)+
                                             (1|pl_species:pl_code),REML = TRUE,data = faa2))
            model_G_sw_log_swc <- try(lmer(G_sw~log_swc+
                                             (log_swc|pl_species)+
                                             (1|pl_species:pl_code),REML = TRUE,data = faa2))
            model_G_sw_log_ppfd <- try(lmer(G_sw~log_ppfd+
                                              (log_ppfd|pl_species)+
                                              (1|pl_species:pl_code),REML = TRUE,data = faa2))
            relgrad <- with(model_G_sw_log@optinfo$derivs,solve(Hessian,gradient))

            if(max(abs(relgrad))>0.0001){
              test_1 <- FALSE
              test_2 <- TRUE}

      }

        
        if(test_1 & !test_2){
          
          type = 2
# 
#           model_G_sw_log_null <- lmer(G_sw~ 1 +
#                                         (1|pl_species),
#                                       data = faa2)
          
          model_G_sw_log <- try(lmer(G_sw~log_vpd_mean+log_swc+log_ppfd+
                                       (log_vpd_mean+log_swc+log_ppfd|pl_species),REML = TRUE,
                                     data = faa2))
          model_G_sw_log_vpd <- try(lmer(G_sw~log_vpd_mean+
                                       (log_vpd_mean|pl_species),REML = TRUE,
                                       data = faa2))
          model_G_sw_log_swc <- try(lmer(G_sw~log_swc+
                                       (log_swc|pl_species),REML = TRUE,
                                       data = faa2))
          model_G_sw_log_ppfd <- try(lmer(G_sw~log_ppfd+
                                       (log_ppfd|pl_species),REML = TRUE,
                                       data = faa2))

        }
        
        if(!test_1 & test_2){

           type = 3
# 
#            model_G_sw_log_null <- lmer(G_sw~ 1 +
#                                          (1|pl_code),
#                                        data = faa2)
           
             model_G_sw_log <- try(lmer(G_sw~log_vpd_mean+log_swc+log_ppfd+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             model_G_sw_log_vpd <- try(lmer(G_sw~log_vpd_mean+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             model_G_sw_log_swc <- try(lmer(G_sw~log_swc+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             model_G_sw_log_ppfd <- try(lmer(G_sw~log_ppfd+
                                          (1|pl_code),REML = TRUE,
                                        data = faa2))
             }  
  
         # }
        
        if(!test_1 & !test_2){

          type = 4

          model_G_sw_log_null <- try(lm(G_sw~1,data = faa2))
          
          model_G_sw_log <- try(lm(G_sw~log_vpd_mean+log_swc+log_ppfd,data = faa2))
          model_G_sw_log_vpd <- try(lm(G_sw~log_vpd_mean,data = faa2))
          model_G_sw_log_swc <- try(lm(G_sw~log_swc,data = faa2))
          model_G_sw_log_ppfd <- try(lm(G_sw~log_ppfd,data = faa2))

          
        } 
        
        lm_model <- lm(G_sw~log_vpd_mean+log_swc+log_ppfd, data = faa2)
        
        # partR2(model_G_sw_log,
        #        partvars = c("log_vpd_mean", "log_swc", "log_ppfd"),
        #        R2_type = "marginal", nboot = 100, CI = 0.95,
        #        data = faa2)
        # r2glmm::r2beta(model_G_sw_log,method="nsj")
        
        print(unique(x2$si_code))
        
        model <- list(env_data,type,faa2,
                      model_G_sw_log,
                      model_G_sw_log_vpd,
                      model_G_sw_log_swc,
                      model_G_sw_log_ppfd,
                      lm_model
                      # model_G_sw_log_null
                      )
        
        save(model,file=paste0(getwd(),"/data/models/complete_G_log/",unique(x2$si_code),".RData"))
      }else{NULL}
      
    }else{NULL}
    
  }else{NULL}
  
  }else{NULL}
  
  rm(x2)
  gc()
  
}
)
