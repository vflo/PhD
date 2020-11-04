# library(sjstats)
library(tidyverse)
library(taxonlookup)
# library(broom)
library(magrittr)
library(rjson)
library(sp)
library(GSIF)
library(AICcmodavg)
library(lmerTest)
# library(medfate)
source("~/Chapter_3_mod/soil_functions.R")
library(purrr)
library(future)
library(furrr)
plan('multisession')
options('future.global.maxsize'=2*1024*1024^2)
# options(scipen = 999)

## 0.1 Plant names to be modelized ---------------

# path <- "~/Chapter_3_mod/data/species_data_subdaily/"
path <- "~/Chapter_3_mod/data/species_data/"
# path <- "~/Chapter_2/data/site_daylight/"
site_names <- list.files(path = path)

# source('~/Chapter_2/Individual_plant_analyses/treatment_list_filter.R')
source('~/Chapter_3_mod/treatment_filter_only_control.R')

load("~/Chapter_2/data/swc_ERA5_land.RData") #nearest cell extract ERA5 land 9x9km 2001 to present

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
site_names <- sample(site_names)

## 0.2 MODELIZATION ----------------------

furrr::future_map(site_names,.progress=TRUE,.f=function(.x){
  # purrr::map(site_names,function(.x){
  ## 1.0 Data load -----------------------
  
  load(paste0(path,.x))
  
  print(.x)
  
  ## 1.1 Conditionals to allow model calculation -----
  x2 <- x2 %>% 
    # mutate(TIMESTAMP_real = TIMESTAMP, 
    #        TIMESTAMP = lubridate::as_date(TIMESTAMP)) %>% 
    split(.[['si_code']],drop = TRUE)
  # 
  purrr::map(x2,function(x2){
    
    print(unique(x2$si_code))
    
    x2 <- x2 %>% mutate(st_treatment = case_when(is.na(st_treatment) ~ "None",
                                                 !is.na(st_treatment) ~ st_treatment))
    if(any(colnames(x2) == "vpd_mean")
    ){ 
      
      if(!any(colnames(x2) == "sw_in_mean")){x2 <- x2 %>% mutate(sw_in_mean = c(NA))}
      if(!any(colnames(x2) == "swc_shallow_mean")){x2 <- x2 %>% mutate(swc_shallow_mean = c(NA))}
      
      metadata <- x2 %>%
        mutate(pl_species = as.character(pl_species),
               TIMESTAMP = lubridate::as_date(TIMESTAMP))
      
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
      
      ERA5_land_data1 <- ERA5_land_data %>% filter(si_code == x2$si_code %>% unique())
      
      x2 %>% right_join(ERA5_land_data1, by = c("TIMESTAMP", "si_code")) %>%
        mutate(swc = ERA5_swc,
               # swc = swc_shallow_mean %>% as.double(),
               # swc = if_else(!is.na(swc), swc_shallow_mean %>% as.double(),
               #               ERA5_swc),
               max_swc = max(swc, na.rm = TRUE),
               min_swc = min(swc, na.rm = TRUE),
               max_min_range = max_swc - min_swc,
               rew = (swc-min_swc)/(max_swc - min_swc)) %>% 
        filter(!is.na(vpd_mean)) %>% 
        left_join(sw, by = c("TIMESTAMP", "si_code")) %>% 
        mutate(TIMESTAMP = lubridate::ymd(TIMESTAMP))%>% 
        mutate(
          ba = ((pl_dbh/2)^2)*pi,
          log_ba = log(ba),
          log_pl_sapw_area = predict.lm(object = dbh_sapw_mod, tibble(log_ba,group)),
          pl_sapw_area = case_when(is.na(pl_sapw_area) ~ exp(log_pl_sapw_area),
                                   !is.na(pl_sapw_area) ~ pl_sapw_area),
          sapflow_mean = case_when(pl_sap_units == "cm3 cm-2 h-1" ~ sapflow_mean,# cm3 cm-2 h-1
                                   pl_sap_units == "cm3 h-1" ~ sapflow_mean/pl_sapw_area),#from cm3 h-1 to cm3 cm2 h-1
          si_elev = coalesce(si_elev,elev)
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
          si_long = mean(si_long, na.rm = TRUE)
        ) %>%
        arrange(TIMESTAMP) -> faa
      
      
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
        mutate(swc_q_50 = quantile(swc,0.5,na.rm = TRUE)) %>% 
        filter(!is.na(G_sw),
               !is.na(swc),
               !is.na(vpd_mean),
               G_sw > 0,
               rew>0,
               swc <= 0.3
               # swc < swc_q_50
        )-> faa
      
      
      # faa %>% 
      #   ungroup() %>% 
      #   summarise(si_code = unique(si_code),
      #             pl_dbh = weighted.mean(pl_dbh,n),
      #             pl_height = weighted.mean(pl_height,n),
      #             si_long = unique(si_long),
      #             si_lat = unique(si_lat),
      #             MAP = unique(MAP),
      #             MAT = unique(MAT)/10,
      #             PET = unique(PET),
      #             PPET = MAP/PET,
      #             D = unique(D),
      #             si_biome = unique(si_biome)
      #   )-> env_data
      
      ## 3.0 GAM calculation -----------------  
      
      faa2 <- faa
      
      # faa %>%
      #   group_by(pl_code) %>%
      # mutate(
      #        quan_filter = 'NO',
      #        quan_filter =
      #          case_when(
      #            G_sw<quantile(G_sw,0.99,na.rm = TRUE) ~ 'YES',
      #            G_sw>quantile(G_sw,0.01,na.rm = TRUE) ~ 'YES')) %>%
      # # mutate(
      # #        # log_swp = log(-swp),
      # #        min_swc = min(swc,na.rm = TRUE),
      # #        max_swc = max(swc,na.rm = TRUE)) %>%
      # filter(quan_filter == "YES") %>%
      # # filter(rew>0) %>%
      # ungroup() ->faa2
      
      
      
      if(nrow(faa2)>1){
        
        faa2%>%
          ungroup() %>% 
          mutate(vpd_cut=cut(vpd_mean, seq(0.3, 9.9, by=0.2)),
                 # swc_cut=cut(swc, seq(0, 1, by=0.05)),
                 # vpd_cut=cut(vpd_mean, 10),
                 # vpd_cut=cut(vpd_mean, seq(0, 10, by=0.3)),
                 # swp_cut=cut(log_swp, seq(unique(min_swp), unique(max_swp), by=1)),
                 # swp_cut=cut(log_swp, c(-12,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,100)),
                 # swp_cut=cut(log_swp, c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,100)),
                 # swc_cut=cut(swc, 5),
                 rew_cut=cut(rew, 5),
                 # vpd_cut_mean = sapply(str_extract_all(vpd_cut,"-?[0-9]+(\\.[0-9]+)?"), 
                 #                       function(x) mean(as.numeric(x))),
                 # rew_cut_mean = sapply(str_extract_all(rew_cut,"-?[0-9]+(\\.[0-9]+)?"),
                 #                       function(x) mean(as.numeric(x))),
                 # cor_vpd_swc = cor(vpd_mean,swc),
                 n_day_plant = n_distinct(vpd_mean)  
          )->faa2
        
        
        if (
          all(n_distinct(faa2$vpd_mean) >= 1,
              n_distinct(faa2$rew) >= 1)
        ){
          set.seed(7)
          
          faa2%>%
            group_by(pl_code,vpd_cut,rew_cut) %>%
            summarise(n_dist_G_cut = n_distinct(G_sw),
                      G_sw_q = quantile(G_sw,prob = 0.95,na.rm = TRUE),
                      G_sw = mean(G_sw,na.rm = TRUE),
                      sapflow_mean = mean(sapflow_mean,na.rm = TRUE),
                      vpd_med = median(vpd_mean,na.rm = TRUE),
                      swc_med = median(swc,na.rm = TRUE),
                      rew_med = median(rew,na.rm = TRUE),
                      vpd_mean = mean(vpd_mean,na.rm = TRUE),
                      swc = mean(swc,na.rm = TRUE),
                      rew = mean(rew,na.rm = TRUE),
                      # rew_cut_mean = mean(rew_cut_mean,na.rm = TRUE),
                      pl_species = unique(pl_species),
                      n_day_plant = unique(n_day_plant),
                      si_code= unique(si_code),
                      pl_height = unique(pl_height),
                      st_height = unique(st_height),
                      si_lat = unique(si_lat),
                      si_long = unique(si_long),
                      si_clay = unique(st_clay_perc),
                      si_sand = unique(st_sand_perc),
                      clay_soil = unique(clay),
                      sand_soil = unique(sand),
                      max_swc = unique(max_swc)
                      # cor_vpd_swc = unique(cor_vpd_swc)
            ) %>%
            ungroup() %>%
            mutate(min_day_plant = min(n_day_plant),
                   swc_min = min(swc),
                   swc_max = max(swc),
                   range_swc = swc_max-swc_min,
                   vpd_min = min(vpd_mean),
                   vpd_max = max(vpd_mean),
                   range_vpd = vpd_max-vpd_min)->faa2
          #         
          faa2 %>%
            left_join(PET %>% select(-c(si_lat,si_long)),by=c('si_code'))%>%
            left_join(BIOS %>% select(-c(si_lat,si_long)),by=c('si_code'))%>%
            group_by(pl_code) %>% 
            mutate(pl_height = coalesce(pl_height, st_height),
                   n = n()) -> faa2
          
        }else{faa2 <- tibble(G_sw = NA,
                             sapflow_mean = NA,
                             n_dist_G_cut = NA,
                             vpd_mean = NA,
                             rew =NA,
                             pl_species = x2$pl_species %>% unique(),
                             pl_code = x2$pl_code %>% unique(),
                             n_day_plant = NA,
                             min_day_plant = NA,
                             si_code= x2$si_code %>% unique())}
        
      }else{faa2 <- tibble(G_sw = NA,
                           sapflow_mean = NA,
                           n_dist_G_cut = NA,
                           vpd_mean = NA,
                           rew =NA,
                           pl_species = x2$pl_species %>% unique(),
                           pl_code = x2$pl_code %>% unique(),
                           n_day_plant = NA,
                           min_day_plant = NA,
                           si_code= x2$si_code %>% unique())}
      
    }else{faa2 <- tibble(G_sw = NA,
                         sapflow_mean = NA,
                         n_dist_G_cut = NA,
                         vpd_mean = NA,
                         rew =NA,
                         pl_species = x2$pl_species %>% unique(),
                         pl_code = x2$pl_code %>% unique(),
                         n_day_plant = NA,
                         min_day_plant = NA,
                         si_code= x2$si_code %>% unique())}
    
    return(faa2)
    
  }) %>% bind_rows()->df1
  return(df1)
}
) %>% bind_rows()-> df      

df2 <- df %>% filter(
  # !pl_code %in% c('RUS_CHE_Y4_Lca_Jt_4',
  #                 'RUS_CHE_Y4_Lca_Jt_5',
  #                 'RUS_CHE_Y4_Lca_Jt_8',
  #                 'RUS_CHE_Y4_Lca_Jt_9'),
  !si_code %in% c('ISR_YAT_YAT'),
  !is.na(vpd_mean),
  log(swc) >= -5 ,
  range_swc>0.05|range_vpd>0.5)#%>%
# dplyr::select(-rew) %>%
# mutate(G_sw=G_sw/1000,
#        G_sw_q=G_sw_q/1000)
# df2 <- df
# df <- df2

save(df2,file = "~/Chapter_3_mod/data/species_data_binned_rew.RData")    
load(file = "~/Chapter_3_mod/data/species_data_binned_rew.RData")  


library(optimx)
model <- lmer(G_sw~log(vpd_mean)+log(rew)+(log(vpd_mean)+log(rew)||pl_species)+(1|pl_species:pl_code),
              data=df2 )

relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

save(model,file='~/Chapter_3_mod/data/model_rew.RData')


model_log <- lmer(log(G_sw)~log(vpd_mean)+log(rew)+(log(vpd_mean)+log(rew)||pl_species)+(1|pl_species:pl_code),
              data=df2 )

relgrad <- with(model_log@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

save(model_log,file='~/Chapter_3_mod/data/model_log_rew.RData')

# 
# 
# 
# 
# 
# 
#         output_lmer <- function( model,faa, env_data, type){ 
#           #model
#           suma_mod <- summary(model)
#           r2 <- r2(model)
#           Vcov <- vcov(model, useScale = FALSE)
#           betas <- fixef(model)
#           se <- sqrt(diag(Vcov))
#           tval <- betas / se
#           pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
#           pred_Gref <- predictSE(model,newdata=tibble(log_vpd = 0,
#                                                       log_rew = 0))
#           se_Gref <- pred_Gref$se.fit
#           Gref <- pred_Gref$fit
#           range_rew <- max(faa$rew,na.rm = TRUE)-min(faa$rew, na.rm = TRUE)
#           mean_rew <- mean(faa$rew,na.rm = TRUE)
#           range_vpd <- max(faa$vpd_mean,na.rm = TRUE)-min(faa$vpd_mean, na.rm = TRUE)
#           mean_vpd <- mean(faa$vpd_mean,na.rm = TRUE)
#           
#           
#           df <- tibble(specie_code = gsub(".RData","",.x),
#                        pl_height = env_data$pl_height,
#                        pl_dbh = env_data$pl_dbh,
#                        MAT =  env_data$MAT,
#                        MAP =  env_data$MAP,
#                        PET =  env_data$PET,
#                        PPET =  env_data$PPET,
#                        D = env_data$D,
#                        si_lat =  env_data$si_lat,
#                        si_long = env_data$si_long,
#                        r2_cond = r2$R2_conditional,  
#                        r2 = r2$R2_marginal,
#                        Gref = Gref,
#                        e_Gref = exp(Gref),
#                        b0 = betas[1],
#                        b1 = betas[2],
#                        b2 = betas[3],
#                        se_Gref = se_Gref,
#                        e_se_Gref = exp(se_Gref),
#                        se_b0 = se[1],
#                        se_b1 = se[2],
#                        se_b2 = se[3],
#                        n = model@frame %>% nrow(),
#                        pval_b0 = pval[[1]],
#                        pval_b1 = pval[[2]],
#                        pval_b2 = pval[[3]],
#                        range_rew = range_rew,
#                        mean_rew = mean_rew,
#                        range_vpd = range_vpd,
#                        mean_vpd = mean_vpd,
#                        median_vpd = median(faa$vpd_mean,na.rm = TRUE),
#                        median_rew = median(faa$rew, na.rm = TRUE),
#                        max_swc = faa$max_swc %>% mean(na.rm = TRUE),
#                        min_swc = faa$min_swc %>% mean(na.rm = TRUE),
#                        max_max_swc = faa$max_swc %>% max(na.rm = TRUE),
#                        min_min_swc = faa$min_swc %>% min(na.rm = TRUE),
#                        model_type = type
#           )
#           
#           return(df)
#           
#         }
#         
#         
#         output_lm <- function(model, faa, env_data, type){ 
#           #model
#           suma_mod <- summary(model)
#           r2 <- r2(model)
#           Vcov <- vcov(model, useScale = FALSE)
#           betas <- coef(model)
#           se <- sqrt(diag(Vcov))
#           tval <- betas / se
#           pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
#           pred_Gref <- predict.lm(model,newdata=tibble(log_vpd = 0,
#                                                        log_rew = 0),
#                                   type='response',
#                                   se.fit = TRUE)
#           Gref <- pred_Gref$fit
#           se_Gref <- pred_Gref$se.fit
#           range_rew <- max(faa$rew,na.rm = TRUE)-min(faa$rew, na.rm = TRUE)
#           mean_rew <- mean(faa$rew,na.rm = TRUE)
#           range_vpd <- max(faa$vpd_mean,na.rm = TRUE)-min(faa$vpd_mean, na.rm = TRUE)
#           mean_vpd <- mean(faa$vpd_mean,na.rm = TRUE)
#           
#           df <- tibble(specie_code = gsub(".RData","",.x),
#                        pl_height = env_data$pl_height,
#                        pl_dbh = env_data$pl_dbh,
#                        MAT =  env_data$MAT,
#                        MAP =  env_data$MAP,
#                        PET =  env_data$PET,
#                        PPET =  env_data$PPET,
#                        D = env_data$D,
#                        si_lat =  env_data$si_lat,
#                        si_long = env_data$si_long,
#                        r2_cond = r2$R2_adjusted,
#                        r2 = r2$R2_adjusted,
#                        Gref = Gref,
#                        e_Gref = exp(Gref),
#                        b0 = betas[1],
#                        b1 = betas[2],
#                        b2 = betas[3],
#                        se_Gref = se_Gref,
#                        e_se_Gref = exp(se_Gref),
#                        se_b0 = se[1],
#                        se_b1 = se[2],
#                        se_b2 = se[3],
#                        n = model$model %>% nrow(),
#                        pval_b0 = pval[[1]],
#                        pval_b1 = pval[[2]],
#                        pval_b2 = pval[[3]],
#                        range_rew = range_rew,
#                        mean_rew = mean_rew,
#                        range_vpd = range_vpd,
#                        mean_vpd = mean_vpd,
#                        median_vpd = median(faa$vpd_mean,na.rm = TRUE),
#                        median_rew = median(faa$rew, na.rm = TRUE),
#                        max_swc = faa$max_swc %>% mean(na.rm = TRUE),
#                        min_swc = faa$min_swc %>% mean(na.rm = TRUE),
#                        max_max_swc = faa$max_swc %>% max(na.rm = TRUE),
#                        min_min_swc = faa$min_swc %>% min(na.rm = TRUE),
#                        model_type = type
#           )
#           
#           return(df)
#           
#         }
#         
#         
#         
#         
#         ## 3.2 without transformation -----------------
#         
#         test_1 <- unique(faa$pl_code) %>% length()>1
#         test_2 <- unique(faa$si_code) %>% length()>1
#         if(test_1 & test_2){
#           test_2_1 <- isSingular(lmer(log_G ~ log_vpd + log_rew + (1|si_code/pl_code), data = faa))
#           test_2_3 <- performance::check_convergence(lmer(log_G ~ log_vpd + log_rew + (1|si_code/pl_code), data = faa))
#           test_2 <- all(!test_2_1, test_2_3)
#         }#more than 1 site and plant
#         
#         
#         if(test_1 & !test_2){ 
#           test_1_1 <- isSingular(lmer(log_G ~ log_vpd + log_rew + (1|pl_code), data = faa))
#           test_1_3 <- performance::check_convergence(lmer(log_G ~ log_vpd + log_rew + (1|pl_code), data = faa))
#           test_1 <- all(!test_1_1, test_1_3)
#         }#one site more than 1 plant
#         
#         if(test_1 & test_2){
#           
#           model <- lmer(log_G ~ log_vpd + log_rew + ( 1|si_code/pl_code), data = faa)
#           
#           type = 1
#           
#           df <- output_lmer(model, faa, env_data, type)
#           
#         }
#         
#         if(test_1 & !test_2){
#           
#           model <- lmer(log_G ~ log_vpd + log_rew + (1|pl_code), data = faa)
#           
#           type = 2
#           
#           df <- output_lmer(model, faa, env_data, type)
#           
#         }
#         if(!test_1 & !test_2){
#           
#           model <- lm(log_G ~ log_vpd + log_rew, data = faa)
#           
#           type = 3
#           
#           df <- output_lm(model, faa, env_data, type)
#           
#         } 
#         
#         ## 3.3 List output creation --------------------------
#         n_trees <- faa %>% dplyr::select(pl_code) %>% unique() %>% nrow()
#         n_days_tree <- faa %>% group_by(pl_code) %>% summarise(n = n())
#         pl_species <- faa %>% summarise(pl_species = unique(pl_species))
#         
#         x_models <- list(
#           model = model,
#           model_output = df,
#           n_trees = n_trees,
#           n_days_tree = n_days_tree,
#           pl_species = pl_species)
#         
#         print(unique(.x))
#         save(x_models,
#              file=paste0(getwd(),
#                          "/data/models_1_0_daylight/",
#                          .x))
#         
# 
# 
