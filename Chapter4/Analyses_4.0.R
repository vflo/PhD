library(sapfluxnetr)
library(tidyverse)
library(taxonlookup)
library(lubridate)
source('PCbiplot.R')
library(caret)
library(lmerTest)
library(magrittr)
library(rlang)
library(GSIF)
library(sp)
library(mgcv)
source('maps_functions.R')
source('calc_relip_mm.R')
library(partR2)
library(extrafont)
# font_import()
loadfonts()
library(ggplot2)
theme_set(theme_bw())
library(viridis)
library(gtable)
library(stringr)
library(raster)
library(ggspatial)
library(rgdal)
library(ggrepel)
library(tidyverse)
library(cowplot)
library(future)
library(furrr)
plan('multisession')
options('future.global.maxsize'=2*1024*1024^2)

rasterOptions(tmpdir = '~/temp' )

custom_palette <- c("#2b99ca","#e6545d", "#ACE4AA", "#40596f", "#f49f27")

# Author: Víctor Flo
# Date: 14/02/2020
# Content: NN models analyses. 
# Addition of taxons of the species using the package "taxonlookup".
# Addition of PET, AI mat, map, biome and soil data per site.
# Joining of RF_data, taxon data,PET, AI, mat, map, biome and soil. Creation of diferent metrics (e.g. SMD's)


#### DATA ####

# COMPLETE G_log MODELS DATA
# path <- "data/models/complete_G_log/"
# site_names <- list.files(path = path)
# models <- purrr::map(site_names,function(x){
#   load(paste0(path,x))
#   print(x)
#   model_env <- model[[1]]
#   model_type <- model[[2]]
#     faa <- model[[4]] %>% summary()
#     n <- faa$residuals %>% length()
# 
#   r2_G_log <- ifelse(class(model[[4]]) == "try-error", as.numeric("NA"),
#                      model[[4]] %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])
# 
#   # rel_r2 <- model[[4]] %>% calc.relip.mm(type='lmg')
#   # rel_r2 <- model[[4]] fii%>% partR2(partvars = c("log_vpd_mean", "log_swc", "log_ppfd"),
#   #                                 R2_type = "conditional", nboot = NULL, CI = 0.95,
#   #                                 max_level = 1, data = model[[3]])
# 
#   rel_r2 <- r2glmm::r2beta(model[[4]], method="nsj")
#   rel_r2 <- rel_r2 %>% as_tibble() %>% dplyr::select(Effect,Rsq)
#   vpd_rel = rel_r2[which(rel_r2$Effect=='log_vpd_mean'),2]$Rsq
#   swc_rel = rel_r2[which(rel_r2$Effect=='log_swc'),2]$Rsq
#   ppfd_rel = rel_r2[which(rel_r2$Effect=='log_ppfd'),2]$Rsq
# 
#   r2_G_log_vpd <- ifelse(class(model[[5]]) == "try-error", as.numeric("NA"),
#                          model[[5]] %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])
# 
#   r2_G_log_swc <- ifelse(class(model[[6]]) == "try-error", as.numeric("NA"),
#                          model[[6]] %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])
# 
#   r2_G_log_ppfd <- ifelse(class(model[[7]]) == "try-error", as.numeric("NA"),
#                           model[[7]] %>% MuMIn::r.squaredGLMM() %>% .[1,"R2c"])
# 
#   df <- tibble(n_days_complete = n,
#                r2_G_log = r2_G_log,
#                r2_G_log_vpd = r2_G_log_vpd,
#                r2_G_log_swc = r2_G_log_swc,
#                r2_G_log_ppfd = r2_G_log_ppfd,
#                model_type_complete = model_type,
#                vpd_rel = vpd_rel,
#                swc_rel = swc_rel,
#                ppfd_rel = ppfd_rel) %>%
#                # si_elev = model_env$si_elev,
#                # pl_sap_units = model_env$pl_sap_units,
#                # swc_shallow = model_env$swc_shallow,
#                # sw_in_mean = model_env$sw_in_mean,
#                # ppfd_in_mean = model_env$ppfd_in_mean,
#                # ppfd = model_env$ppfd,
#                # ppfd_sw_in = model_env$ppfd_sw_in
#     #            )
#     cbind(model_env %>% dplyr::select(-MAT,-MAP)
#           )
#   return(df)
# }) %>% bind_rows()
# save(models,file = "data/complete_G_log_results.RData")
load("data/complete_G_log_results.RData")

# path <- "data/models/complete_bin_gam_models_2/"
# site_names <- list.files(path = path)
# complete_gam <- purrr::map(site_names,function(x){
#   load(paste0(path,x))
#   print(x)
#   model_type <- model[[6]]
#   faa <- model[[1]] %>% summary()
#   r2_GAM_comp <- faa[['r.sq']]
#   faa <- model[[2]] %>% summary()
#   r2_GAM_vpd <- faa[['r.sq']]
#   faa <- model[[3]] %>% summary()
#   r2_GAM_swc <- faa[['r.sq']]
#   faa <- model[[4]] %>% summary()
#   r2_GAM_ppfd <- faa[['r.sq']]
#     df <- tibble(r2_GAM_comp = r2_GAM_comp,
#                  r2_GAM_vpd = r2_GAM_vpd,
#                  r2_GAM_swc = r2_GAM_swc,
#                  r2_GAM_ppfd = r2_GAM_ppfd,
#                  model_type_complete = model_type,
#                  si_code = x)
#     return(df)
#   }) %>% bind_rows()
# complete_gam <- complete_gam %>% mutate(si_code = str_remove(si_code,".RData"))
# save(complete_gam,file = "data/complete_gam_results.RData")
load("data/complete_gam_results.RData")


## import sites metadata
# folder_sapwood <- "~/sapfluxnet_db/0.1.4/RData/sapwood"
# 
# sfn_metadata_sapwood <- read_sfn_metadata(folder = folder_sapwood, .write_cache = TRUE)
# sapwood_sites <- sfn_metadata_sapwood$site_md$si_code 
# 
# plant_md <- sfn_metadata_sapwood[["plant_md"]]
# site_md <- sfn_metadata_sapwood[["site_md"]]
# stand_md <- sfn_metadata_sapwood[["stand_md"]]
# species_md <- sfn_metadata_sapwood[["species_md"]] %>%
#   dplyr::mutate(pl_species = .data$sp_name)
# env_md <- sfn_metadata_sapwood[["env_md"]]
# 
# metadata_sapwood <- plant_md %>%
#   dplyr::left_join(site_md,by = "si_code") %>%
#   dplyr::left_join(stand_md, by = "si_code") %>%
#   dplyr::left_join(species_md, by = c("si_code", "pl_species")) %>%
#   dplyr::left_join(env_md, by = "si_code") %>%
#   dplyr::select(.data$si_code,.data$pl_code,dplyr::everything())
# 
# 
# folder_plant <- "~/sapfluxnet_db/0.1.4/RData/plant"
# 
# sfn_metadata_plant <- read_sfn_metadata(folder = folder_plant, .write_cache = TRUE)
# 
# plant_sites <- sfn_metadata_plant$site_md$si_code 
# 
# plant_sites <- plant_sites[!plant_sites%in%sapwood_sites]
# 
# plant_md <- sfn_metadata_plant[["plant_md"]]
# site_md <- sfn_metadata_plant[["site_md"]]
# stand_md <- sfn_metadata_plant[["stand_md"]]
# species_md <- sfn_metadata_plant[["species_md"]] %>%
#   dplyr::mutate(pl_species = .data$sp_name)
# env_md <- sfn_metadata_plant[["env_md"]]
# 
# metadata_plant <- plant_md %>%
#   dplyr::left_join(site_md,by = "si_code") %>%
#   dplyr::left_join(stand_md, by = "si_code") %>%
#   dplyr::left_join(species_md, by = c("si_code", "pl_species")) %>%
#   dplyr::left_join(env_md, by = "si_code") %>%
#   dplyr::select(.data$si_code,.data$pl_code,dplyr::everything()) %>% 
#   filter(si_code %in% plant_sites)

# metadata <- bind_rows(metadata_sapwood, metadata_plant)
# save(metadata, file="~/Chapter4/data/metadata.RData")
load("~/Chapter4/data/metadata.RData")


## import CHELSA BIOS

chelsa <- read_csv("data/clim_biomes_chelsa.csv") %>%
  dplyr::select(si_code, si_mat, si_map, si_biome_chelsa = si_biome)

biomes <- read.csv("data/biome.csv")


BIOS <- read_csv("~/Chapter4/data/CHELSA_BIOS.csv") %>%
  dplyr::select(-lon,-lat,MAT = MAT, MAP = MAP,BIO_4,BIO_15,si_code) %>% 
  mutate(BIO_4 = BIO_4/1000,
         MAT = MAT/10)


load("~/Chapter4/data/LAI_google.RData")
LAI_google <- LAI
colnames(LAI_google) <- c("si_code","lon","lat","LAI_google")
load("~/Chapter4/data/LAI.RData")
load("~/Chapter4/data/canopy_height.RData")

## import Data Soil
load("~/Chapter4/data/clay.RData")
load("~/Chapter4/data/sand.RData")
load("~/Chapter4/data/nitrogen.RData")
load("~/Chapter4/data/ph.RData")
load("~/Chapter4/data/bedrock.RData")

## Joining and mutate
data <- models %>% 
  right_join(BIOS, by=c('si_code'))%>%
  left_join(complete_gam %>% dplyr::select(-model_type_complete), by = "si_code") %>% 
  left_join(clay %>% dplyr::select(si_code,clay),by=c('si_code')) %>% 
  left_join(sand %>% dplyr::select(si_code,sand),by=c('si_code')) %>% 
  left_join(nitrogen %>% dplyr::select(si_code,nitrogen),by=c('si_code')) %>%
  left_join(ph %>% dplyr::select(si_code,ph),by=c('si_code')) %>% 
  left_join(bedrock %>% dplyr::select(si_code,bedrock),by=c('si_code')) %>% 
  left_join(metadata %>% 
              dplyr::select(st_clay_perc,st_silt_perc,st_sand_perc,
                            si_code,st_height,st_lai) %>% 
              unique()) %>% 
  left_join(chelsa,by = 'si_code') %>%
  left_join(biomes,by = 'si_code') %>%
  left_join(LAI_google %>% ungroup() %>% 
              dplyr::select(si_code,LAI_google),by = 'si_code') %>% 
  left_join(LAI %>% ungroup()
            %>% dplyr::select(si_code,sd_LAI,LAI),by = 'si_code') %>% 
  left_join(canopy_height) %>%
  # left_join(TOTN, by = "si_code") %>%
  mutate(clay = coalesce(st_clay_perc, clay),
         sand = coalesce(st_sand_perc, sand),
         st_height1 = st_height,
         st_height = coalesce(st_height,pl_height),
         st_height = coalesce(st_height,canopy_height),
         lai = coalesce(st_lai,LAI_google),
         LAI = coalesce(lai,LAI))
         # MAT_bios = MAT_bios/10,
         # MAP = coalesce(MAP,MAP_bios),
         # MAT = coalesce(MAT,MAT_bios))

data %>% split(.[['si_code']],drop = FALSE) %>%
  purrr::map(function(x){
    x %>%
    mutate(si_biome_cal = ggbiome::gd_get_biome(si_mat = .$MAT,
                                                si_map = .$MAP,
                                                merge_deserts = TRUE)$si_biome %>%
             as.character())->faa
  return(faa)
}) %>% bind_rows() -> data

data_plant <- data %>% 
  mutate(PETP = 1/PPET,
         P_PET = MAP-PET,
         si_biome = recode(biome,
                           `wood` = 'WOOD',
                           `temp` = "TEMP",
                           `bor` = "BOR",
                           `tro` = "TROP",
                           `dry` = 'DRY')
         )


gA <- biome_plot(merge_biomes = TRUE)+
  geom_point(data=data_plant,aes(x=MAP,y=MAT),fill="white",shape=21)+
  scale_fill_manual(name = 'BIOME',
                    values = custom_palette)
gg_biomes <- plot_grid(gA,ncol = 1)

pdf(file ="gg_biomes.pdf",width=6,height=5)
gg_biomes
dev.off()

# gB <- ggbiome::vis_biome()+
#   geom_point(data=data_plant,aes(x=MAP,y=MAT),fill="white",shape=21)
# 
# gg_biomes <- plot_grid(gA,gB,ncol = 2,labels = "auto",align="hv")
# 
# pdf(file ="gg_biomes.pdf",width=11,height=3)
# gg_biomes
# dev.off()


library(ggpmisc)
source('color_branded.R')
dbh_sapw_mod <- readRDS(file="~/Chapter_3_mod/dbh_sapw_area_model.rds")
my.formula <- y ~ x
dbh_sapw_mod$model %>% 
  ggplot(aes(x=log_ba,y=log_pl_sapw_area,group=group,fill=group,color=group))+
  geom_point(shape=21)+
  geom_smooth(method = "lm")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  scale_fill_manual(values = c("grey10","grey50"))+
  scale_color_manual(values = c("grey10","grey50"))+
  # scale_color_branded(primary="green",other="yellow")+
  xlab("ln(basal area)")+
  ylab("ln(sapwood area)")+
  theme(legend.position = c(0.8, 0.2),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.title = element_blank()) -> dbh_sw

pdf(file ="dbh_sw.pdf",width=6,height=4)
dbh_sw
dev.off()

path <- "data/models/complete_G_log/"
site_names <- list.files(path = path)

summary_data <-purrr::map(site_names,function(x){
  load(paste0(path,x))
  print(x)
  faa <- model[[3]] %>% 
    cbind(si_code = gsub(".RData","",x)) %>% 
    dplyr::select(si_code,pl_code,pl_species) %>% 
    group_by(si_code) %>% 
    summarise(`n trees` = pl_code %>% unique() %>% length(),
           `n species` = pl_species%>% unique() %>% length())
  return(faa)
  }) %>% bind_rows()

resume <- data_plant %>% 
  left_join(summary_data) %>% 
  mutate(si_lat = round(si_lat,2) %>% as.character(),
         si_long = round(si_long,2) %>% as.character()) %>% 
  dplyr::select(si_code,
                si_lat,
                si_long,
                si_biome,
                n_days_complete,
                `n species`,
                `n trees`) %>% 
  filter(!is.na(n_days_complete))

write_csv(resume,path="resume_table.csv")


resume2 <- data_plant %>% 
  dplyr::select(si_code,
                `$R[2]_{vpd}$` = r2_G_log_vpd,
                `$R[2]_{swc}$` = r2_G_log_swc,
                `$R[2]_{ppfd}$` =r2_G_log_ppfd,
                `Relimp VPD` = vpd_rel,
                `Relimp SWC` = swc_rel,
                `Relimp PPFD` = ppfd_rel,
                MAT = MAT,
                MATsd = BIO_4,
                MAP = MAP,
                MAPsd = BIO_15,
                Bedrock = bedrock,
                Clay = clay,
                Sand = sand,
                pH = ph,
                `Total N` = nitrogen,
                `Stand height` = st_height,
                LAI = LAI,
                st_clay_perc,
                st_sand_perc,
                st_height1 = st_height1,
                pl_height,
                st_lai)  %>% 
  filter(!is.na(`$R[2]_{vpd}$`))%>% 
  dplyr::mutate(`$R[2]_{vpd}$` = format(round(`$R[2]_{vpd}$`,2), nsmall = 2),
                `$R[2]_{swc}$` = format(round(`$R[2]_{swc}$`,2), nsmall = 2),
                `$R[2]_{ppfd}$` = format(round(`$R[2]_{ppfd}$`,2), nsmall = 2),
                `Relimp VPD` = format(round(`Relimp VPD`,2), nsmall = 2),
                `Relimp SWC` = format(round(`Relimp SWC`,2), nsmall = 2),
                `Relimp PPFD` = format(round(`Relimp PPFD`,2), nsmall = 2),
                MAT = format(round(MAT,2), nsmall = 2),
                MATsd = format(round(MATsd,2), nsmall = 2),
                MAP = format(round(MAP,2), nsmall = 2),
                MAPsd = format(round(MAPsd,2), nsmall = 2),
                Bedrock = format(round(Bedrock,2), nsmall = 2),
                Clay = paste0(format(round(Clay,2), nsmall = 2),
                              ifelse(is.na(st_clay_perc)," b", " a")),
                Sand = paste0(format(round(Sand,2), nsmall = 2),
                              ifelse(is.na(st_sand_perc)," b", " a")),
                pH = format(round(pH,2), nsmall = 2),
                `Total N` = format(round(`Total N`,2), nsmall = 2),
                height = case_when(!is.na(st_height1)~as.character(" a") ,
                                   is.na(st_height1)&!is.na(pl_height)~as.character(" c"),
                                   is.na(st_height1)&is.na(pl_height)~as.character(" b")),
                `Stand height` = paste0(format(round(`Stand height`,2), nsmall = 2),
                                        as.character(height)),
                LAI = paste0(format(round(LAI,2), nsmall = 2),
                             ifelse(is.na(st_lai)," b", " a"))) %>%
  dplyr::select(-height,-pl_height,-st_lai,-st_height1,-st_clay_perc,-st_sand_perc)

write_csv(resume2,path="resume_table2.csv")

# Hidrometeorological correlations ---------------

# library(broom)
# path <- "data/models/complete_G_log/"
# # path2 <- "data/models_raw/complete_bin_gam_models_raw/"
# site_names <- list.files(path = path)
# raw_data <-purrr::map(site_names,function(x){
#   load(paste0(path,x))
#   print(x)
#   faa <- model[[3]] %>% 
#     dplyr::select(G_sw,vpd_mean,log_swc,ppfd,pl_code) %>% 
#     cbind(si_code = gsub(".RData","",x)) %>% 
#     distinct() %>% 
#     mutate(n = n())
#   faa %>%
#     group_by(pl_code) %>% 
#     do(glance(lm(log(vpd_mean)~log(ppfd), data = .))) %>%
#     dplyr::select(r2=r.squared) %>% 
#     ungroup() %>% 
#     summarise(r2=mean(r2))->fii
#   faa <- faa %>% cbind(r2=fii$r2)
#   
#   
#   return(faa)}) %>% bind_rows()
# 
# raw_data %>%
#   ggplot(aes(log(vpd_mean),log_swc))+
#   geom_point()
# 
# raw_data %>% 
#   ggplot(aes(log(vpd_mean),log(ppfd)))+
#   geom_point()
# 
# raw_data %>% 
#   ggplot(aes(log(ppfd),log_swc))+
#   geom_point()
# 
# lmer(log(vpd_mean)~log_swc+(log_swc|si_code),data=raw_data) %>% 
#   MuMIn::r.squaredGLMM()
# 
# lmer(log(vpd_mean)~log(ppfd)+(log(ppfd)|si_code/pl_code),data=raw_data) %>% 
#   MuMIn::r.squaredGLMM()
# 
# lmer(log(ppfd)~log_swc+(log_swc|si_code),data=raw_data) %>% 
#   MuMIn::r.squaredGLMM()


#### BIOME PLOTS ####

library(gridExtra)


theme0 <- function(...) theme( legend.position = "none",
                               panel.background = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.spacing = unit(0,"null"),
                               axis.ticks = element_blank(),
                               axis.text.x = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.ticks.length = unit(0,"null"),
                               axis.text = element_text(margin = unit(0,"null")),
                               panel.border=element_rect(color=NA),...)


#### VPD SM RAD G_LOG

p_vpd_swc <-ggplot(data_plant %>% filter(!is.na(si_biome)),aes(y=r2_G_log_vpd,x=r2_G_log_swc)) +
  geom_point(aes(color = si_biome),alpha = 1,show.legend=FALSE) +
  # stat_density2d(aes(fill=si_biome,alpha = (..level..)^0.5),
  #                geom = "polygon",
  #                bins = 12, 
  #                show.legend = TRUE) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),linetype = 2)+
  scale_x_continuous(expand = c(0.05,0)) +
  scale_y_continuous(expand = c(0.05,0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"points"),
        legend.position = c(0.7, 0.9), 
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18),
        axis.text = element_text(size=16),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_alpha_continuous(range = c(0, 1),
                         guide = 'none'
  )+
  scale_fill_manual(guide = 'none',
                    values = custom_palette)+
  scale_color_manual(guide = 'none',
                     values = custom_palette)+
  labs(y = expression(~R[VPD]^2),
       x = NULL,
       fill = 'BIOME')+ 
  xlim(-0.1, 1.1) + ylim(-0.1,1.1)


p_par_swc <-ggplot(data_plant%>% filter(!is.na(si_biome)),aes(y = r2_G_log_ppfd, 
                                                              x = r2_G_log_swc)) +
  geom_point(aes(color = si_biome), alpha = 1,show.legend=FALSE) +
  # stat_density2d(aes(fill = si_biome, alpha = (..level..)^0.5),
  #                geom = "polygon",
  #                bins = 12, 
  #                show.legend = FALSE) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),linetype = 2)+
  scale_x_continuous(expand = c(0.05,0)) +
  scale_y_continuous(expand = c(0.05,0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0,0,0,0),"points"),
        legend.position = c(0.7, 0.9), 
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=16),
        legend.title = element_text(size=18),
        axis.text = element_text(size=16))+
  scale_alpha_continuous(range = c(0, 1),
                         guide = 'none')+
  scale_fill_manual(values = custom_palette)+
  scale_color_manual(guide = 'none',
                     values = custom_palette)+
  labs(y =expression(~R[PPFD]^2),
       x =expression(~R[SWC]^2),
       fill = 'BIOME')+ 
  xlim(-0.1,1.1) + ylim(-0.1,1.1)


p_par_vpd <-ggplot(data_plant%>% filter(!is.na(si_biome)),aes(x=r2_G_log_ppfd,y=r2_G_log_vpd)) +
  geom_point(aes(color = si_biome),alpha = 1,show.legend=FALSE) +
  # stat_density2d(aes(fill=si_biome,alpha = (..level..)^0.5),
  #                geom = "polygon",
  #                bins = 12, 
  #                show.legend=FALSE) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),linetype = 2)+
  scale_x_continuous(expand=c(0.05,0)) +
  scale_y_continuous(expand=c(0.05,0)) +
  theme_bw() +
  theme(plot.margin=unit(c(0,0,0,0),"points"),
        legend.position = c(0.7, 0.9), 
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18),
        axis.text = element_text(size=16),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_alpha_continuous(range = c(0, 1),
                         guide = 'none'
  )+
  scale_fill_manual(values=custom_palette)+
  scale_color_manual(guide = 'none',
                     values=custom_palette)+
  labs( y = NULL, x =expression(~R[PPFD]^2))+ 
  xlim(-0.1,1.1) + ylim(-0.1,1.1)

p_par <- ggplot(data_plant%>% filter(!is.na(si_biome)),aes(x=r2_G_log_ppfd,colour = si_biome)) + 
  geom_density(alpha=0.5) +
  scale_x_continuous(breaks=NULL,expand=c(0.02,0)) +
  scale_y_continuous(breaks=NULL,expand=c(0.02,0)) +
  theme_bw() +
  theme(plot.margin=unit(c(0,0,0,0),"points"),
        legend.position = c(0.7, 0.9), 
        legend.justification = c(0, 1),
        plot.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18),
        axis.text = element_text(size=16),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_color_manual(guide = 'none',
                     values=custom_palette)+
  theme0()+
  xlim(-0.1,1.1)

p_swc <- ggplot(data_plant%>% filter(!is.na(si_biome)),aes(x=r2_G_log_swc,colour = si_biome)) + 
  geom_density(alpha=0.5) + 
  scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
  scale_y_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
  theme_bw() +
  theme(plot.margin=unit(c(0,0,0,0),"points"),
        legend.position = c(0.7, 0.9), 
        legend.justification = c(0, 1),
        plot.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18),
        axis.text = element_text(size=16),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_color_manual(guide = 'none',
                     values=custom_palette)+
  theme0()+
  xlim(-0.1,1.1)

p_vpd <- ggplot(data_plant%>% filter(!is.na(si_biome)),aes(x=r2_G_log_vpd,colour = si_biome)) + 
  geom_density(alpha=0.5) + 
  coord_flip()+
  scale_x_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
  scale_y_continuous(labels = NULL,breaks=NULL,expand=c(0.02,0)) +
  theme_bw() +
  theme(plot.margin=unit(c(0,0,0,0),"points"),
        legend.position = c(0.7, 0.9), 
        legend.justification = c(0, 1),
        plot.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=16),
        legend.title = element_text(size=18),
        axis.text = element_text(size=16),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_color_manual(guide = 'none',
                     values=custom_palette)+
  theme0()+
  xlim(-0.1,1.1)


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend_plot <- ggplot(data_plant%>% filter(!is.na(si_biome)) %>%rename(BIOME=si_biome),aes(x=r2_G_log_swc,fill = BIOME)) +
  geom_density()+
  theme_bw()+
  scale_fill_manual(values=custom_palette)+
  # theme(legend.title = element_text('BIOME'))+
  NULL

mylegend<-g_legend(legend_plot)

library(patchwork)
library(egg)
library(cowplot)
library(grid)

p_vpd_swc <- set_panel_size(p_vpd_swc,
                            width  = unit(6, "cm"),
                            height = unit(6, "cm"))

p_par_vpd <- set_panel_size(p_par_vpd,
                            width  = unit(6, "cm"),
                            height = unit(6, "cm"))

p_par_swc <- set_panel_size(p_par_swc,
                            width  = unit(6, "cm"),
                            height = unit(6, "cm"))

p_par <- set_panel_size(p_par,
                        width  = unit(6, "cm"),
                        height = unit(2, "cm"))

p_swc <- set_panel_size(p_swc,
                        width  = unit(6, "cm"),
                        height = unit(2, "cm"))

p_vpd <- set_panel_size(p_vpd,
                        width  = unit(2, "cm"),
                        height = unit(6, "cm"))

lege <- set_panel_size(g=mylegend,
                       width  = unit(6, "cm"),
                       height = unit(6, "cm"))

tab <- gtable(unit(rep(1, 137), "mm"), unit(rep(1, 150), "mm"))
tab <- gtable_add_grob(tab,p_vpd_swc, t = 48, l = 25)
tab <- gtable_add_grob(tab,p_par_vpd, t = 54, l = 95)
tab <- gtable_add_grob(tab,p_par_swc, t = 116, l = 25)
tab <- gtable_add_grob(tab,p_swc, t = 5, l = 33)
tab <- gtable_add_grob(tab,p_vpd, t = 48, l = 137)
tab <- gtable_add_grob(tab,p_par, t = 5, l = 96)
tab <- gtable_add_grob(tab,lege, t = 115, l = 100)
grid.newpage()
dim(tab)
pdf(file ="r2_G_log_test.pdf",width=6.5,height=6.5)
grid.arrange(tab) %>% grid.draw()
dev.off()



### G_LOG COMPARISON of individual variables ---------------------------

mod1 <- lm((r2_G_log_vpd-r2_G_log_swc)~si_biome , data = data_plant, weights = n_days_complete)
summary(mod1)
multcomp::cld(mod1 %>% emmeans::emmeans('si_biome'), 
              Letters = LETTERS) -> vpd_swc
emmeans::test(mod1 %>% emmeans::emmeans('si_biome'),0) -> vpd_swc_p
vpd_swc_p <- vpd_swc_p[order(factor(vpd_swc_p$si_biome, levels=unique(vpd_swc$si_biome))),]

mod1 <- lm((r2_G_log_vpd-r2_G_log_ppfd)~si_biome , data = data_plant, weights = n_days_complete)
summary(mod1)
multcomp::cld(mod1 %>% emmeans::emmeans('si_biome'), 
              Letters = LETTERS) -> vpd_ppfd
emmeans::test(mod1 %>% emmeans::emmeans('si_biome'),0) -> vpd_ppfd_p
vpd_ppfd_p <- vpd_ppfd_p[order(factor(vpd_ppfd_p$si_biome, levels=unique(vpd_ppfd$si_biome))),]

mod1 <- lm((r2_G_log_swc-r2_G_log_ppfd)~si_biome , data = data_plant, weights = n_days_complete)
summary(mod1)
multcomp::cld(mod1 %>% emmeans::emmeans('si_biome'), 
              Letters = LETTERS) -> swc_ppfd
emmeans::test(mod1 %>% emmeans::emmeans('si_biome'),0) -> swc_ppfd_p
swc_ppfd_p <- swc_ppfd_p[order(factor(swc_ppfd_p$si_biome, levels=unique(swc_ppfd$si_biome))),]


mod1 <- lm((r2_G_log_swc)~si_biome , data = data_plant, weights = n_days_complete)
summary(mod1)
anova(mod1)
multcomp::cld(mod1 %>% emmeans::emmeans('si_biome'), 
              Letters = LETTERS) -> swc_alone
emmeans::test(mod1 %>% emmeans::emmeans('si_biome'),0)


mod1 <- lm((r2_G_log_vpd)~si_biome , data = data_plant, weights = n_days_complete)
summary(mod1)
anova(mod1)
multcomp::cld(mod1 %>% emmeans::emmeans('si_biome'), 
              Letters = LETTERS) -> vpd_alone
emmeans::test(mod1 %>% emmeans::emmeans('si_biome'),0)


mod1 <- lm((r2_G_log_ppfd) ~ si_biome , data = data_plant, weights = n_days_complete)
summary(mod1)
anova(mod1)
multcomp::cld(mod1 %>% emmeans::emmeans('si_biome'), 
              Letters = LETTERS) -> ppfd_alone
emmeans::test(mod1 %>% emmeans::emmeans('si_biome'),0)
# plot(mod1)


mod1 <- lm(r2_G_log~si_biome , data = data_plant, weights = n_days_complete)
summary(mod1)
anova(mod1)
multcomp::cld(mod1 %>% emmeans::emmeans('si_biome'), 
              Letters = LETTERS) -> complete_alone
emmeans::test(mod1 %>% emmeans::emmeans('si_biome'),0)


df_comp1 <- tibble(Biome = vpd_swc$si_biome %>% as.character(),
                   vpd_swc = paste0(round(vpd_swc$emmean,3),
                                    vpd_swc$.group,
                                    " ",
                                   symnum(vpd_swc_p$p.value,corr = FALSE,na = FALSE,
                                         cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                                         symbols = c("***", "**", "*", ".", " "))
                                   ))
df_comp2 <- tibble(Biome = vpd_ppfd$si_biome %>% as.character(),
                  vpd_ppfd = paste0(round(vpd_ppfd$emmean,3),
                                   symnum(vpd_ppfd_p$p.value,corr = FALSE,na = FALSE,
                                          cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                                          symbols = c("***", "**", "*", ".", " ")),
                                   " ",
                                   vpd_ppfd$.group))
df_comp3 <- tibble(Biome = swc_ppfd$si_biome %>% as.character(),
                  swc_ppfd = paste0(round(swc_ppfd$emmean,3),
                                   symnum(swc_ppfd_p$p.value,corr = FALSE,na = FALSE,
                                          cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                                          symbols = c("***", "**", "*", ".", " ")),
                                   " ",
                                   swc_ppfd$.group))
df_comp4 <- tibble(Biome = vpd_alone$si_biome %>% as.character(),
                  vpd_alone = paste0(round(vpd_alone$emmean,3),
                                   " ",
                                   vpd_alone$.group))
df_comp5 <- tibble(Biome = swc_alone$si_biome %>% as.character(),
                  swc_alone = paste0(round(swc_alone$emmean,3),
                                   " ",
                                   swc_alone$.group))
df_comp6 <- tibble(Biome = ppfd_alone$si_biome %>% as.character(),                  
                   ppfd_alone = paste0(round(ppfd_alone$emmean,3),
                                   " ",
                                   ppfd_alone$.group))
df_comp7 <- tibble(Biome = complete_alone$si_biome %>% as.character(),
                  complete_alone = paste0(round(complete_alone$emmean,3),
                                   " ",
                                   complete_alone$.group)
                  )

df_comp <- df_comp4 %>% 
  left_join(df_comp5) %>% 
  left_join(df_comp6) %>% 
  left_join(df_comp7) %>% 
  left_join(df_comp1) %>% 
  left_join(df_comp2) %>% 
  left_join(df_comp3)

write_csv(df_comp,path="comparison_biome.csv")

#### BIOCLIMATIC MODELS --------------------------------------------
#### R2 COMPLETE MODEL 
fee <- data_plant %>%
  filter(!is.na(r2_G_log)) %>%
  dplyr::select(r2_G_log_vpd,PPET,st_height,MAT,MAP,clay,sand,
                nitrogen,LAI,bedrock,ph,BIO_4,BIO_15)

df2 = cor(fee, use = "complete.obs")
hc = caret::findCorrelation(df2, cutoff=0.7,names=TRUE) # putt any value as a "cutoff"
hc = sort(hc)
hc
# reduced_Data = fee[,hc]
# # print (reduced_Data)
# pairs(data_plant %>% dplyr::select(PPET,st_height,MAT,MAP,clay,sand,MO,LAI,TOTN))
lm(r2_G_log*100~(MAP+st_height+MAT+clay+nitrogen + ph+LAI+BIO_4+BIO_15),
   data=data_plant, weights = n_days_complete)->fuu
fuu %>% summary()
step(fuu) -> mod1
mod1 %>% summary()
car::vif(mod = mod1 )
anova(mod1)
dominanceanalysis::dominanceAnalysis(mod1)



#### R2 VPD NN
lm(r2_G_log_vpd*100~(MAP+st_height+MAT+clay+nitrogen + ph+LAI+BIO_4+BIO_15),
   data=data_plant, weights = n_days_complete)->fuu
fuu %>% summary()
step(fuu) -> mod2
mod2 %>% summary()
car::vif(mod2)
anova(mod2)
dominanceanalysis::dominanceAnalysis(mod2)


#### R2 SWC NN
lm(r2_G_log_swc*100~(MAP+st_height+MAT+clay+nitrogen + ph+LAI+BIO_4+BIO_15),
   data=data_plant, weights = n_days_complete)->fuu
fuu %>% summary()
step(fuu) -> mod3
mod3 %>% summary()
car::vif(mod3)
anova(mod3)
dominanceanalysis::dominanceAnalysis(mod3)


#### R2 PPFD NN
lm(r2_G_log_ppfd*100~(MAP+st_height+MAT+clay+nitrogen + ph+LAI+BIO_4+BIO_15+abs(lat)),
   data=data_plant, weights = n_days_complete)->fuu
fuu %>% summary()
step(fuu) -> mod4
mod4 %>% summary()
car::vif(mod4)
anova(mod4)
dominanceanalysis::dominanceAnalysis(mod4)
# 
# 


# #### R2 VPD - SWC NN
lm(vpd_rel~(MAP+st_height+MAT+clay+nitrogen + ph+LAI+BIO_4+BIO_15),
   data=data_plant, weights = n_days_complete)->fuu
fuu %>% summary()
step(fuu) -> mod5
mod5 %>% summary()
# car::vif(mod5)
# anova(mod5)
# 
# #### R2 VPD - PPFD NN
lm(swc_rel~(MAP+st_height+MAT+clay+nitrogen + ph+LAI+BIO_4+BIO_15),
   data=data_plant, weights = n_days_complete)->fuu
fuu %>% summary()
step(fuu) -> mod6
mod6 %>% summary()
# car::vif(mod6)
# anova(mod6)
# 
# #### R2 SWC - PPFD NN
lm(ppfd_rel~(MAP+st_height+MAT+clay+nitrogen + ph+LAI+BIO_4+BIO_15),
   data=data_plant, weights = n_days_complete)->fuu
fuu %>% summary()
step(fuu) -> mod7
mod7 %>% summary()
# car::vif(mod7)
# anova(mod7)

m1 <- summary(mod1)$coefficients
m1 <- m1 %>% 
  as_tibble(rownames = "Variable") %>% 
  dplyr::select(Variable,Estimate,p_val=`Pr(>|t|)`) %>% 
  mutate(esti = paste(round(Estimate,3),
                      symnum(p_val,corr = FALSE,na = FALSE,
                             cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                             symbols = c("***", "**", "*", ".", " ")))) %>% 
  dplyr::select(Variable,esti) %>% gather(key = key, value = value, 2:ncol(.)) %>% 
  spread(key = names(.)[1], value = "value") %>% 
  mutate(Variable = "VPD+SWC+PPFD")


m2 <- summary(mod2)$coefficients
m2 <- m2 %>% 
  as_tibble(rownames = "Variable") %>% 
  dplyr::select(Variable,Estimate,p_val=`Pr(>|t|)`) %>% 
  mutate(esti = paste(round(Estimate,3),
                      symnum(p_val,corr = FALSE,na = FALSE,
                             cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                             symbols = c("***", "**", "*", ".", " ")))) %>% 
  dplyr::select(Variable,esti) %>% gather(key = key, value = value, 2:ncol(.)) %>% 
  spread(key = names(.)[1], value = "value") %>% 
  mutate(Variable = "VPD")

m3 <- summary(mod3)$coefficients
m3 <- m3 %>% 
  as_tibble(rownames = "Variable") %>% 
  dplyr::select(Variable,Estimate,p_val=`Pr(>|t|)`) %>% 
  mutate(esti = paste(round(Estimate,3),
                      symnum(p_val,corr = FALSE,na = FALSE,
                             cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                             symbols = c("***", "**", "*", ".", " ")))) %>% 
  dplyr::select(Variable,esti) %>% gather(key = key, value = value, 2:ncol(.)) %>% 
  spread(key = names(.)[1], value = "value")%>% 
  mutate(Variable = "SWC")

m4 <- summary(mod4)$coefficients
m4 <- m4 %>% 
  as_tibble(rownames = "Variable") %>% 
  dplyr::select(Variable,Estimate,p_val=`Pr(>|t|)`) %>% 
  mutate(esti = paste(round(Estimate,3),
                      symnum(p_val,corr = FALSE,na = FALSE,
                             cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1),
                             symbols = c("***", "**", "*", ".", " ")))) %>% 
  dplyr::select(Variable,esti) %>% gather(key = key, value = value, 2:ncol(.)) %>% 
  spread(key = names(.)[1], value = "value")%>% 
  mutate(Variable = "PPFD")

mod_r2 <- tibble(r2 = c(summary(mod1)$r.squared %>% round(3),
                        summary(mod2)$r.squared %>% round(3),
                        summary(mod3)$r.squared %>% round(3),
                        summary(mod4)$r.squared %>% round(3)))

m_table <- m1 %>% bind_rows(m2) %>% bind_rows(m3) %>% bind_rows(m4) %>% 
  dplyr::select(Variable,`(Intercept)`,MAT, MATsd = BIO_4,MAP, MAPsd = BIO_15,
                clay,nitrogen,ph,
                st_height,LAI) %>% cbind(mod_r2)

write_csv(m_table,path="models_table.csv")



##### SPAGUETTI PLOTS ------------------------------------------------------


##### G LOG SPAGUETTI PLOTS ------------------------------------------------------
detach("package:lmerTest", unload=TRUE)
# detach("package:raster", unload=TRUE)
library(bootpredictlme4)
options(bootnsim = 100)
library(optimx)

neworder <- c("DRY","WOOD","TEMP",'BOR','TROP')


## VPD

path <- "data/models/complete_G_log/"
# path2 <- "data/models_raw/complete_bin_gam_models_raw/"
site_names <- list.files(path = path)
pred_G_log_vpd <-purrr::map(site_names,function(x){
  load(paste0(path,x))
  print(x)
  faa <- model[[3]] %>% 
    dplyr::select(G_sw,log_vpd_mean,pl_code,pl_species) %>% 
    mutate(G_sw = round(G_sw,6),
           log_vpd_mean = round(log_vpd_mean,6),
           pl_code = as.factor(pl_code),
           pl_species = as.factor(pl_species))
  data <- data_plant %>% filter(si_code == gsub(".RData","",x)) 
  model <- model[[4]]
  df <- if(attr(model, "class")=="lm"){df <- model$model}else{df <- model@frame}
  df <- df %>% rename("G_sw" = `G_sw`,"log_vpd_mean"= `log_vpd_mean`,
                      "log_swc"= `log_swc`,'log_ppfd'= `log_ppfd`)
  df <- df %>% mutate(G_sw = (G_sw),log_vpd_mean = log_vpd_mean,
                      log_swc = log(0.5),log_ppfd = log(0.75),
                      G_sw = round(G_sw,6), log_vpd_mean = round(log_vpd_mean,6))
  predict <- predict(model,newdata = df)
  df <- df %>% cbind(G_pred = predict) %>% 
    cbind(data %>% dplyr::select(-ppfd)) %>% 
    left_join(faa) 
  return(df)
}) %>% 
  bind_rows()%>% 
  filter(G_pred > 0) %>% 
  group_by(si_biome) %>% 
  mutate(max_vpd = max(log_vpd_mean) %>% unique()) %>%
  ungroup()

mod_G_log_vpd <- lmer(G_pred~log_vpd_mean*si_biome+(1|si_code/pl_code),
                          weights=r2_G_log, data = pred_G_log_vpd %>% 
                            mutate(r2_G_log = case_when(r2_G_log<0~0,
                                                        r2_G_log>=0~r2_G_log),
                                   si_biome = as.factor(si_biome)),
                          REML = TRUE, 
                          control = lmerControl(optimizer ="Nelder_Mead"))


save(mod_G_log_vpd,file="~/Chapter4/data/spa_model_G_log_vpd.RData")

gg2_G_log_vpd <- pred_G_log_vpd %>% 
  split(.[['si_biome']],drop = FALSE) %>%
  purrr::map(function(x){
    print(x$si_biome %>% unique())
    vpd_sim <- seq(0.3,exp(x$max_vpd) %>% max(),0.05)
    vpd_df <- tibble(log_vpd_mean=log(vpd_sim),si_biome = x$si_biome %>% unique())
    res <- predict(object = mod_G_log_vpd, 
                   newdata = vpd_df,
                   re.form=NA,
                   se.fit=TRUE
    )
    res <- tibble(vpd = vpd_sim ,
                  # predi = res,
                  predi = res$fit,
                  se = res$se.fit,
                  si_biome = x$si_biome %>% unique())
    return(res)}) %>% bind_rows()


save(gg2_G_log_vpd,file="~/Chapter4/data/gg2_G_log_vpd.RData")
load(file="~/Chapter4/data/gg2_G_log_vpd.RData")

fuu <- arrange(mutate(pred_G_log_vpd,si_biome=factor(si_biome,levels=neworder)),si_biome)

fuu %>% 
  filter(!is.na(si_biome),G_pred<800) %>%
  # group_by(pl_code)%>%
  # arrange(vpd_mean,G_pred, .by_group = TRUE) %>%
  # ungroup() %>%
  ggplot()+
  geom_line(aes(x=exp(log_vpd_mean),y=G_pred,group = pl_code,alpha = 0.1*r2_G_log),
            color = 'grey', show.legend =FALSE)+
  geom_line(data = gg2_G_log_vpd,aes(x=vpd,y=predi,group = si_biome, color = si_biome),
            size=1.5,show.legend=FALSE)+
  ggplot2::geom_line(data = gg2_G_log_vpd, aes(x = vpd, y = (predi + se),group = si_biome),
                     size = 0.3,linetype=2) +
  ggplot2::geom_line(data = gg2_G_log_vpd, aes(x = vpd, y = (predi - se),group = si_biome), 
                     size = 0.3,linetype=2) +
  facet_wrap(.~si_biome,scale = "free",nrow = 1)+
  scale_color_manual(values=custom_palette)+
  labs(
    y = "",
    # y = expression(paste(G[s], " [mol ",m^{-2}, " ", s^{-1},"]" )), 
    x = expression(paste('VPD [kPa]')))+
  guides(color=guide_legend(title="Biome"))+
  scale_alpha_continuous(range = c(0.1, 0.6),
                         guide = 'none'
  )+
  ylim(0,800)+
  theme(strip.background =element_blank(),
        legend.position = c(0.85, 0.25),
        plot.margin = unit(c(0,0,0,0),"points")) -> spa1

## SWC

pred_G_log_swc <-purrr::map(site_names,function(x){
  load(paste0(path,x))
  print(x)
  faa <- model[[3]] %>% 
    dplyr::select(G_sw,log_swc,pl_code,pl_species) %>% 
    mutate(G_sw = round(G_sw,6),
           log_swc = round(log_swc,6),
           pl_code = as.factor(pl_code),
           pl_species = as.factor(pl_species))
  data <- data_plant %>% filter(si_code == gsub(".RData","",x)) 
  model <- model[[4]]
  df <- if(attr(model, "class")=="lm"){df <- model$model}else{df <- model@frame}
  df <- df %>% rename("G_sw" = `G_sw`,"log_vpd_mean"= `log_vpd_mean`,"log_swc"= `log_swc`,'log_ppfd'= `log_ppfd`)
  df <- df %>% mutate(G_sw = (G_sw),
                      log_swc = log_swc, log_vpd_mean = log(1), log_ppfd = log(0.75),
                      G_sw = round(G_sw,6), log_swc = round(log_swc,6))
  predict <- predict(model,newdata = df)
  df <- df %>% cbind(G_pred = predict) %>% 
    cbind(data %>% dplyr::select(-ppfd)) %>% 
    left_join(faa) 
  return(df)
}) %>% 
  bind_rows()%>% 
  filter(G_pred > 0) %>% 
  group_by(si_biome) %>% 
  mutate(max_swc = max(exp(log_swc)) %>% unique(),
         min_swc = min(exp(log_swc)) %>% unique()) %>%
  ungroup()

mod_G_log_swc <- lmer(G_pred~log_swc*si_biome+(1|si_code/pl_code),
                          weights=r2_G_log, data = pred_G_log_swc %>% 
                            mutate(r2_G_log = case_when(r2_G_log<0~0,
                                                            r2_G_log>=0~r2_G_log)),
                          REML = TRUE, 
                          control = lmerControl(optimizer ="Nelder_Mead"))


save(mod_G_log_swc,file="~/Chapter4/data/spa_model_G_log_swc.RData")


gg2_G_log_swc <- pred_G_log_swc %>% 
  split(.[['si_biome']],drop = FALSE) %>%
  purrr::map(function(x){
    swc_sim <- seq(x$min_swc %>% min(),x$max_swc %>% max(),0.01)
    swc_df <- tibble(log_swc=log(swc_sim),
                     # swc = swc_sim,
                     si_biome = x$si_biome %>% unique())
    res <- predict(object = mod_G_log_swc, 
                   newdata = swc_df,re.form=NA,
                   se.fit=TRUE)
    res <- tibble(swc = swc_sim ,
                  predi = res$fit,
                  # predi = res,
                  se = res$se.fit,
                  si_biome = x$si_biome %>% unique())
    return(res)}) %>% bind_rows()

save(gg2_G_log_swc,file="~/Chapter4/data/gg2_G_log_swc.RData")
load(file="~/Chapter4/data/gg2_G_log_swc.RData")

fuu <- arrange(mutate(pred_G_log_swc,si_biome=factor(si_biome,levels=neworder)),si_biome)

fuu %>% 
  filter(!is.na(si_biome)) %>%
  # group_by(pl_code)%>%
  # arrange(swc,G_pred, .by_group = TRUE) %>% 
  # ungroup() %>% 
  ggplot()+
  geom_line(aes(x=exp(log_swc),y=(G_pred),group = pl_code,alpha = 0.1*r2_G_log),
            color = 'grey', show.legend =FALSE)+
  geom_line(data = gg2_G_log_swc,
            aes(x=(swc),y=(predi),group = si_biome,
                color = si_biome),size=1.5,show.legend=FALSE)+
  ggplot2::geom_line(data = gg2_G_log_swc, aes(x = swc, y = (predi + se),group = si_biome), size = 0.3,linetype=2) +
  ggplot2::geom_line(data = gg2_G_log_swc, aes(x = swc, y = (predi - se),group = si_biome), size = 0.3,linetype=2) +
  facet_wrap(.~si_biome,scales = 'free',nrow = 1)+
  scale_color_manual(values=custom_palette)+
  labs(y = expression(paste(G[Asw], " [mol ",m^{-2}, " ", s^{-1},"]" )), 
       x = expression(paste('SWC [',m^{3},  " ",m^{-3},"]")))+
  guides(color=guide_legend(title="Biome"))+
  scale_alpha_continuous(range = c(0.1, 0.6),
                         guide = 'none'
  )+
  ylim(0,800)+
  theme(strip.background =element_blank(),
        strip.text.x = element_blank(),
        legend.position = c(0.85, 0.25),
        plot.margin = unit(c(0,0,0,0),"points")) -> spa2


## PPFD
pred_G_log_ppfd <- purrr::map(site_names,function(x){
  load(paste0(path,x))
  print(x)
  faa <- model[[3]] %>% 
    dplyr::select(G_sw,log_ppfd,pl_code,pl_species) %>% 
    mutate(G_sw = round(G_sw,6),
           log_ppfd = round(log_ppfd,6),
           pl_code = as.factor(pl_code),
           pl_species = as.factor(pl_species))
  data <- data_plant %>% filter(si_code == gsub(".RData","",x)) 
  model <- model[[4]]
  df <- if(attr(model, "class")=="lm"){df <- model$model}else{df <- model@frame}
  df <- df %>% rename("G_sw" = `G_sw`,"log_vpd_mean"= `log_vpd_mean`,"log_swc"= `log_swc`,'log_ppfd'= `log_ppfd`)
  df <- df %>% mutate(G_sw = (G_sw),
                      log_ppfd = log_ppfd, log_vpd_mean = log(1), log_swc = log(0.5),
                      G_sw = round(G_sw,6), log_ppfd = round(log_ppfd,6))
  predict <- predict(model,newdata = df)
  df <- df %>% cbind(G_pred = predict) %>% 
    cbind(data %>% dplyr::select(-ppfd)) %>% 
    left_join(faa) 
  return(df)
}) %>% 
  bind_rows()%>% 
  filter(G_pred > 0) %>% 
  group_by(si_biome) %>% 
  mutate(max_ppfd = max(exp(log_ppfd)) %>% unique(),
         min_ppfd = min(exp(log_ppfd)) %>% unique()) %>%
  ungroup()

mod_G_log_ppfd <- lmer(G_pred~log_ppfd*si_biome+(1|si_code/pl_code),
                           weights=r2_G_log, data = pred_G_log_ppfd %>% 
                             mutate(r2_G_log = case_when(r2_G_log<0~0,
                                                             r2_G_log>=0~r2_G_log)),
                           REML = TRUE, 
                           control = lmerControl(optimizer ="Nelder_Mead"))


save(mod_G_log_ppfd,file="~/Chapter4/data/spa_model_G_log_ppfd.RData")

gg2_G_log_ppfd <- pred_G_log_ppfd %>% 
  split(.[['si_biome']],drop = FALSE) %>%
  purrr::map(function(x){
    print(x$si_biome %>% unique())
    ppfd_sim <- seq(x$min_ppfd %>% min(),x$max_ppfd %>% max(),0.02)
    ppfd_df <- tibble(log_ppfd = log(ppfd_sim),si_biome = x$si_biome %>% unique())
    res <- predict(object = mod_G_log_ppfd, 
                   newdata = ppfd_df,re.form=NA,
                   se.fit=TRUE)
    res <- tibble(ppfd = ppfd_sim ,
                  predi = res$fit,
                  # predi = res,
                  se = res$se.fit,
                  si_biome = x$si_biome %>% unique())
    return(res)}) %>% bind_rows()

save(gg2_G_log_ppfd,file="~/Chapter4/data/gg2_G_log_ppfd.RData")
load(file="~/Chapter4/data/gg2_G_log_ppfd.RData")

fuu <- arrange(mutate(pred_G_log_ppfd,si_biome=factor(si_biome,levels=neworder)),si_biome)

fuu %>% 
  filter(!is.na(si_biome)) %>%
  group_by(pl_code)%>%
  # arrange(ppfd,G_pred, .by_group = TRUE) %>% 
  ungroup() %>% 
  ggplot()+
  geom_line(aes(x=exp(log_ppfd),y=(G_pred),group = pl_code,alpha = 0.1*r2_G_log),
            color = 'grey', show.legend =FALSE)+
  geom_line(data = gg2_G_log_ppfd,aes(x=ppfd,y=(predi),group = si_biome, color = si_biome),
            size=1.5,show.legend=FALSE)+
  ggplot2::geom_line(data = gg2_G_log_ppfd, aes(x = ppfd, y = predi + se,group = si_biome), 
                     size = 0.3,linetype=2) +
  ggplot2::geom_line(data = gg2_G_log_ppfd, aes(x = ppfd, y = predi - se,group = si_biome), 
                     size = 0.3,linetype=2) +
  facet_wrap(.~si_biome,scales = 'free',nrow = 1)+
  scale_color_manual(values=custom_palette)+
  labs(y="",
    # y = expression(paste(G[s], " [mol ",m^{-2}, " ", s^{-1},"]" )), 
    x = expression(paste('PPFD [',mu,"mol ",m^{-2},  " ",s^{-1},"]"))
       )+
  guides(color=guide_legend(title="Biome"))+
  scale_alpha_continuous(range = c(0.1, 0.6),
                         guide = 'none'
  )+
  ylim(0,800)+
  theme(strip.background =element_blank(),
        strip.text.x = element_blank(),
        legend.position = c(0.85, 0.25),
        plot.margin = unit(c(0,0,0,0),"points")) -> spa3


# spa_plots <- plot_grid(spa1,spa2,spa3,ncol = 1)
spa_plots <-cowplot::align_plots(spa1,spa2,spa3, align = "hv")
spa_plots2 <- plot_grid(ggdraw(spa_plots[[1]]),
                        ggdraw(spa_plots[[2]]),
                        ggdraw(spa_plots[[3]]),ncol = 1)

pdf(file ="spa_plots.pdf",width=12,height=6.5)
spa_plots2
dev.off()





#### MAPS -------------------------------------------------------------

library(RColorBrewer)

crowther <- raster('~/[databases]/Crowther/Crowther_Nature_Files_Revision_01_WGS84_GeoTiff/Crowther_Nature_Biome_Revision_01_WGS84_GeoTiff.tif')

MAP <- raster("~/Chapter4/data/rasters/MAP_resample.tif")
MAT <- raster("~/Chapter4/data/rasters/MAT_resample.tif")
MAT <- MAT/10
BIO_4 <- raster("~/Chapter4/data/rasters/BIO_4_resample.tif")
BIO_15 <- raster("~/Chapter4/data/rasters/BIO_15_resample.tif")
BIO_4 <- BIO_4/1000
LAI <- raster("~/Chapter4/data/rasters/LAI_google_resample.tif")
st_height <- raster("~/Chapter4/data/rasters/st_height_resample.tif")
clay <- raster("~/Chapter4/data/rasters/clay_resample.tif")
sand <- raster("~/Chapter4/data/rasters/sand_resample.tif")
bedrock <- raster("~/Chapter4/data/rasters/bedrock_resample.tif")
nitrogen <- raster("~/Chapter4/data/rasters/nitrogen_resample.tif")
ph <- raster("~/Chapter4/data/rasters/ph_resample.tif")

data_stack <- stack(st_height,MAT,MAP,clay,sand,
                    LAI,nitrogen,ph,bedrock,
                    BIO_4,BIO_15)

names(data_stack) <- c('st_height','MAT','MAP','clay','sand',
                       'LAI','nitrogen','ph','bedrock',
                       'BIO_4','BIO_15')

# 1.2. Mask preparation -------------------------------------------

# crowther_mask <- crowther
# crowther_mask[crowther_mask < 400] <- NA
# gplot_wrld_r <- gplot_data(crowther)

# 2.0 PREDICT Complete ---------------
# 
map01 <- datamapcal(data_stack, mod1,st_height)
map02 <- datamapcal(data_stack, mod2,st_height)
map03 <- datamapcal(data_stack, mod3,st_height)
map04 <- datamapcal(data_stack, mod4,st_height)

# 3.1 Predict plot -----------------

g1 <- maps_plot(map01, title = "VPD + SWC + PPFD")
g2 <- maps_plot(map02,title = "VPD ")
g3 <- maps_plot(map03,title = "SWC ")
g4 <- maps_plot(map04,title = "PPFD ")+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.03, 0.65),
        legend.background = element_blank(),
        legend.key.size = unit(0.45, "cm"),
        legend.key.width = unit(0.45,"cm"),
        legend.key = element_blank(),
        legend.title = element_blank())


# 3.2 Error plot -----------------

library('bootpredictlme4')
zeta <- as(data_stack, "SpatialPixelsDataFrame") %>% data.frame()
zeta1 <- zeta[complete.cases(zeta),]

mapse01 <- maperrorpred(zeta1, mod1, map01, st_height)
mapse02 <- maperrorpred(zeta1, mod2, map02, st_height)
mapse03 <- maperrorpred(zeta1, mod3, map03, st_height)
mapse04 <- maperrorpred(zeta1, mod4, map04, st_height)

se1 <- map_error_plot(datamaperror = mapse01, title = "VPD + SWC + PPFD")
se2 <- map_error_plot(datamaperror = mapse02, title = "VPD")
se3 <- map_error_plot(datamaperror = mapse03, title = "SWC")
se4 <- map_error_plot(datamaperror = mapse04, title = "PPFD")+
  theme(legend.justification=c(0,1), 
        legend.position=c(0.03, 0.65),
        legend.background = element_blank(),
        legend.key.size = unit(0.45, "cm"),
        legend.key.width = unit(0.45,"cm"),
        legend.key = element_blank(),
        legend.title = element_blank())

gg_maps <- plot_grid(g2,se2,g3,se3,g4,se4,ncol = 2)

pdf(file ="G_log_maps.pdf",width=13,height=9.5)
gg_maps
dev.off()


#### VARIABLES
data_stack$MAP[data_stack$MAP >= 3000]<-3000
MAP <- maps_plot_var(data_stack$MAP,title = "MAP")
MAPsd <- maps_plot_var(data_stack$BIO_15,title = "MAPsd")
MAT <- maps_plot_var(data_stack$MAT,title = "MAT")
MATsd <- maps_plot_var(data_stack$BIO_4,title = "MATsd")
# Bedrock <- maps_plot_var(data_stack$bedrock,title = "Bedrock")
Clay <- maps_plot_var(data_stack$clay,title = "Clay")
Sand <- maps_plot_var(data_stack$sand,title = "Sand")
pH <- maps_plot_var(data_stack$ph,title = "pH")
data_stack$nitrogen[data_stack$nitrogen >= 10]<-10
TotN <- maps_plot_var(data_stack$nitrogen,title = "Total N")
data_stack$st_height[data_stack$st_height <= 0.5]<-NA
St_height <- maps_plot_var(data_stack$st_height,title = "Stand height")
LAI <- maps_plot_var(data_stack$LAI,title = "LAI")
gg_maps_var <- plot_grid(MAT,MATsd,MAP,MAPsd,Clay,Sand,pH,TotN,St_height,LAI,ncol = 2)
pdf(file ="var_maps.pdf",width=8.5,height=10)
gg_maps_var
dev.off()

# RELATIVE IMPORTANCE MAP ------------------
library(ggtern)
library(tricolore)

map05 <- datamapcal(data_stack, mod5,st_height)
map06 <- datamapcal(data_stack, mod6,st_height)
map07 <- datamapcal(data_stack, mod7,st_height)

fuu2 <- as.data.frame(map05)$layer
fuu3 <- as.data.frame(map06)$layer
fuu4 <- as.data.frame(map07)$layer

fuuu <- cbind(vpd=fuu2,swc=fuu3,ppfd=fuu4) %>% as_tibble()

trico <- tricolore::Tricolore(fuuu,p1 = 'vpd', p2 = 'swc', p3 = 'ppfd',
                            breaks = Inf, contrast = 0.5,lightness = 1,
                            chroma = 1, spread =1, show_data = FALSE,
                            show_center = FALSE, label_as = 'pct',
                            crop = FALSE)

tri <- coordinates(map02) %>% cbind(colortri = trico$rgb) %>%as_tibble()
tri <- tri %>% mutate(x=as.numeric(x),
                      y=as.numeric(y))
# coordinates(tri) <- ~ x + y
# # coerce to SpatialPixelsDataFrame
# gridded(tri) <- TRUE
# # coerce to raster
# rasterDF <- raster(tri)
# rasterDF
# map08 <- map05
# map08@data@values <- NA
triplot <- ggplot() +
  # ggspatial::layer_spatial(map05,fill="transparent") +
  geom_sf(data=world,colour="gray70", fill="gray70",size = 0.1)+
  geom_tile(data = tri, mapping = aes(x = x, y = y,fill=colortri)) +
  geom_point(data = data_plant,aes(x=si_long,y=si_lat), shape = 1)+
  scale_fill_identity(na.value = "transparent")+
  xlab("Lon") + ylab("Lat") + 
  coord_sf(crs = 4326, xlim=c(-180, 180), ylim=c(-60,85),expand=FALSE)+
  NULL


triplot +
  annotation_custom(
    ggplotGrob(
      trico$key+
        geom_Tline(Tintercept=.33, colour='grey') + 
        geom_Lline(Lintercept=.33, colour='grey') + 
        geom_Rline(Rintercept=.33, colour='grey')+
        theme_showgrid_major()+
        theme_showarrows()+
        theme(plot.background = element_rect( color = NA),text = element_text(size=10))+
        labs(L = '% VPD', T = '% SWC', R = '% PPFD')),
    xmin = -180, xmax = -100, ymin = -60, ymax = 20)+
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background=element_rect(fill = "transparent",colour = NA),
        plot.margin=unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))-> triplot_map

triplot_map <- plot_grid(triplot_map)

pdf(file ="triplot_map.pdf",width=10,height=5)
triplot_map
dev.off()