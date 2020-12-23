library(sapfluxnetr)
library(tidyverse)
library(magrittr)
library(taxonlookup)


# 1. Read metadata ---------------------------------------------------------------
# From previously written cache file

sfn_metadata_plant <- read_sfn_metadata(folder = '~/sapfluxnet_db/0.1.5/RData/plant', .write_cache = FALSE)
sfn_metadata_sapwood <- read_sfn_metadata(folder = '~/sapfluxnet_db/0.1.5/RData/sapwood', .write_cache = FALSE)
sfn_metadata_leaf <- read_sfn_metadata(folder = '~/sapfluxnet_db/0.1.5/RData/leaf', .write_cache = FALSE)

# 2. Aggregate all datasets -----------------------------------------------

# Join all metadata regardless of having sap flow per sapwood or per plant

sfn_allsites<- sfn_metadata_plant[['site_md']] %>% 
  full_join(dplyr::select(sfn_metadata_sapwood[['site_md']],-si_remarks))

sfn_allstands<- sfn_metadata_plant[['stand_md']] %>% 
  full_join(sfn_metadata_sapwood[['stand_md']])

sfn_allplants<- sfn_metadata_plant[['plant_md']] %>% 
  full_join(sfn_metadata_sapwood[['plant_md']]) %>% 
  distinct(pl_code,.keep_all = TRUE)

sfn_env <- sfn_metadata_plant[['env_md']] %>% 
  full_join(sfn_metadata_sapwood[['env_md']])


sfn_metadata_plant[['species_md']] %>% 
  group_by(si_code,sp_name) %>% 
  mutate(n_sp=sum(sp_ntrees)) -> sfn_sitesp_plant

sfn_metadata_sapwood[['species_md']] %>% 
  group_by(si_code,sp_name) %>% 
  mutate(n_sp=sum(sp_ntrees)) -> sfn_sitesp_sw


# 3. Fix errors and taxonize ----------------------------------------------


# Species level - plant
sfn_sitesp_plant%>% 
  ungroup() %>% 
  mutate(sp_name = case_when(
    sp_name =='Larix sibirica Ledeb.'~'Larix sibirica',
    sp_name == 'Leptolaena' ~'Leptolaena sp.',
    sp_name == 'Eschweillera sp.' ~'Eschweilera sp.',
    sp_name == 'Vacapoua americana' ~ 'Vouacapoua americana',
    sp_name == 'Brachulaena ramiflora' ~ 'Brachylaena ramiflora',
    sp_name == 'Cryptocaria spp.' ~ 'Cryptocarya sp.',
    TRUE ~ sp_name)) ->  sfn_sitespecies_plant_fix

sfn_sitespecies_plant_fix %>% 
pull(sp_name) %>%
  unique() %>%
  taxonlookup::lookup_table(missing_action = 'NA', by_species = TRUE) %>%
  rownames_to_column('sp_name') %>%
  left_join(sfn_sitespecies_plant_fix, ., by = 'sp_name') -> sfn_sitespecies_plant_tax1

# Species level - plant
sfn_sitesp_sw%>% 
  ungroup() %>% 
  mutate(sp_name = case_when(
    sp_name =='Larix sibirica Ledeb.'~'Larix sibirica',
    sp_name == 'Leptolaena' ~'Leptolaena sp.',
    sp_name == 'Eschweillera sp.' ~'Eschweilera sp.',
    sp_name == 'Vacapoua americana' ~ 'Vouacapoua americana',
    sp_name == 'Brachulaena ramiflora' ~ 'Brachylaena ramiflora',
    sp_name == 'Cryptocaria spp.' ~ 'Cryptocarya sp.',
    TRUE ~ sp_name)) ->  sfn_sitespecies_sw_fix

sfn_sitespecies_sw_fix %>% 
  pull(sp_name) %>%
  unique() %>%
  taxonlookup::lookup_table(missing_action = 'NA', by_species = TRUE) %>%
  rownames_to_column('sp_name') %>%
  left_join(sfn_sitespecies_sw_fix, ., by = 'sp_name') -> sfn_sitespecies_sw_tax1 


n_all_species<- union(sfn_sitespecies_plant_tax1$sp_name,sfn_sitespecies_sw_tax1$sp_name)

sfn_sitespecies_tax1<- dplyr::union(sfn_sitespecies_plant_tax1,sfn_sitespecies_sw_tax1)

# sfn_sitespecies_plant_tax %>% 
#   group_by(sp_name) %>% 
#   mutate(n_sp=sum(sp_ntrees)) %>% 
#   distinct(sp_name,.keep_all = TRUE) %>% 
#   arrange(desc(n_sp)) %>% View()
# 

# I had to do this bc of changes in tibble 3.0.0
# https://twitter.com/rushworth_a/status/1246719797955567623

sfn_sitespecies_tax1 <- as.data.frame(sfn_sitespecies_tax1)

sfn_sitespecies_tax1[sfn_sitespecies_tax1$sp_name == 'Myrtaceae sp.',
                     c('sp_name','genus','order','family','group')] <-
  c('Myrtaceae fam.','Myrtaceae fam.','Myrtaceae','Myrtales','Angiosperms')

sfn_sitespecies_tax1 <- as_tibble(sfn_sitespecies_tax1)


# Plant level
sfn_allplants %>% 
  mutate(pl_species = case_when(
    pl_species =='Larix sibirica Ledeb.'~'Larix sibirica',
    pl_species == 'Leptolaena spp.' ~'Leptolaena sp.',
    pl_species == 'Eschweillera sp.' ~'Eschweilera sp.',
    pl_species == 'Vacapoua americana' ~ 'Vouacapoua americana',
    pl_species == 'Brachulaena ramiflora' ~ 'Brachylaena ramiflora',
    pl_species == 'Cryptocaria spp.' ~ 'Cryptocarya sp.',
    TRUE ~ pl_species)) ->  sfn_allplants_fix

sfn_allplants_fix %>% 
pull(pl_species) %>%
  unique() %>%
  taxonlookup::lookup_table(missing_action = 'NA', by_species = TRUE) %>%
  rownames_to_column('pl_species') %>%
  left_join(sfn_allplants_fix, ., by = 'pl_species') -> sfn_allplants_tax

sfn_allplants_tax[sfn_allplants_tax$pl_species == 'Myrtaceae sp.',
                    c('pl_species','genus','order','family','group')] <-
  tibble(pl_species=rep('Myrtaceae fam.',6),genus=rep('Myrtaceae fam.',6),
         family=rep('Myrtaceae',6),order=rep('Myrtales',6),group=rep('Angiosperms',6))

# The following code fixes data for sapwood depth and sapwood area 
# USA_CHE_ASP: sapwood depth was in mm, convert to cm
# USA_SMI*: sapwood area in m2, convert to cm2

# sfn_allplants_tax %<>% 
# mutate(
#   pl_sapw_depth = ifelse(si_code=='USA_CHE_ASP',pl_sapw_depth*0.1,pl_sapw_depth),
#   pl_sapw_area = ifelse(si_code=='USA_SMI_SCB' | si_code=='USA_SMI_SER',pl_sapw_area*1E4,pl_sapw_area)
# ) 

# Fix number of plants per species/site
# Use plant level data to correct errors in number of 
# plants per species 

sfn_allplants_tax %>% 
  group_by(si_code,pl_species) %>% tally() %>% 
  right_join(sfn_sitespecies_tax1,by=c('si_code','pl_species'='sp_name')) %>% 
  mutate(sp_ntrees=n) %>% 
  dplyr::select(-sp_ntrees,-n_sp) %>% 
  rename(sp_name=pl_species,sp_ntrees=n) %>% 
  ungroup() ->sfn_sitespecies_tax

# 
# sfn_sitespecies_tax %>% 
#   group_by(sp_name) %>% 
#   summarise(nt=sum(sp_ntrees)) %>% View()


  # 4. Number of trees (per dataset, species, etc.) -------------------------

# Number of distinct sites
# based on dataset coding
sfn_allstands %>% 
  mutate(country_code=sapply(str_split(si_code,"_"),"[[",1),
         site_code=sapply(str_split(si_code,"_"),"[[",2),
         sites_geo=paste(country_code,site_code,sep='_')) %>% 
  distinct(sites_geo) %>% tally() ->n_sites_codes

# based on coordinates
sfn_allsites %>% 
  distinct(si_lat,si_long,.keep_all = TRUE) %>% tally() ->n_sites_coords

# TODO:

sfn_allsites %>% 
  separate(si_code,sep='_',into=c('country','site','stand','treatment'),
           remove=FALSE) ->n_sites_extracted


# Datasets per biome

sfn_allsites %>% 
  group_by(si_biome) %>% tally()
  

# number of trees and species per dataset, with coordinates
sfn_sites_nspecies <- sfn_sitespecies_tax %>% 
  group_by(si_code) %>% 
  summarise(nspecies=length(sp_name),
            ntrees=sum(sp_ntrees)) %>% 
  left_join(sfn_allsites %>% dplyr::select(si_code,si_lat,si_long))

# trees per dataset
sfn_sites_nspecies %>% 
  arrange(desc(ntrees)) ->sfn_sites_nspecies_dataset

sfn_sites_nspecies %>% 
  arrange(desc(ntrees)) %>% pull(ntrees) %>% quantile(c(0.25,0.5,0.75,.95))

sfn_sites_nspecies %>% 
  filter(ntrees>=4) %>% tally()

# Number of datasets per species
sfn_nspecies_dataset <- sfn_sitespecies_tax %>% 
  group_by(sp_name) %>% 
  tally() %>% 
  rename(n_sites=n) %>% 
  arrange(desc(n_sites)) 


# Number of species, taxonomic detail
sfn_nspecies_taxdetail<- sfn_sitespecies_tax %>% 
  distinct(sp_name,.keep_all=TRUE) %>% 
  dplyr::select(group, order, family,sp_name) 

# sfn_nspecies_taxdetail %>% 
#   group_by(group) %>% 
#   tally()


# Number of trees per species
sfn_species_ntrees<- sfn_allplants_tax %>% 
  group_by(pl_species) %>% 
  mutate(n_trees=n()) %>% 
  distinct(pl_species,.keep_all = TRUE) %>% 
  mutate(species = case_when(
    pl_species=='Myrtaceae sp.'~ paste0(family,' fam.'),
    TRUE ~ pl_species)) %>% 
  ungroup() %>% 
  dplyr::select(group,species,n_trees) %>% 
  arrange(desc(n_trees)) 


sfn_groups_ntrees<- sfn_allplants_tax %>% 
  group_by(group) %>% 
  tally()

# Number of trees per genus

sfn_genus_ntrees<- sfn_allplants_tax %>% 
  group_by(genus) %>% 
  mutate(n_trees=n()) %>% 
  distinct(genus,.keep_all = TRUE) %>% 
  mutate(genus_f = case_when(
    genus == 'Unknown'~ paste0(family,' fam.'),
    TRUE ~ genus)) %>% 
  ungroup() %>% 
  dplyr::select(genus_f,n_trees) %>% 
  arrange(desc(n_trees)) 

# 3. Methodologies -------------------------------------------------
# Levels: plant, sapwood, leaf

sfn_sites_plsw <- sfn_metadata_plant[['site_md']] %>% 
  semi_join( dplyr::select(
    sfn_metadata_sapwood[['site_md']],
    -si_remarks)) %>% 
  mutate(type='plant,sapwood')


sfn_sites_pl <- sfn_metadata_plant[['site_md']] %>% 
  anti_join(dplyr::select(
    sfn_metadata_sapwood[['site_md']],
    -si_remarks)) %>% 
  mutate(type='plant')


sfn_sites_sw <- dplyr::select(sfn_metadata_sapwood[['site_md']],-si_remarks) %>% 
  anti_join(sfn_metadata_plant[['site_md']]) %>% 
  mutate(type='sapwood')

sfn_sites_leaf <- sfn_metadata_leaf[['site_md']] %>% 
  mutate(type='leaf,plant,sapwood')

sfn_sites_alllevels <- sfn_sites_plsw %>% 
  full_join(sfn_sites_pl) %>% 
  full_join(sfn_sites_sw)  


# Measurement type: plant, sapwood, leaf
sfn_sites_type <- sfn_sites_alllevels %>% 
  mutate(
    type=ifelse(si_code%in%sfn_sites_leaf$si_code,
                'leaf,plant,sapwood',type),
    typef=factor(type)
  ) 


# Measurement method

sfn_method_ntrees <- sfn_allplants_tax %>% 
  group_by(pl_sens_meth) %>% 
  summarise(n_trees = n()) %>% 
  mutate(perc_trees = n_trees/sum(n_trees)*100) %>% 
  dplyr::select(pl_sens_meth,n_trees,perc_trees) %>% 
  arrange(desc(perc_trees))


# For alluvial plot
sfn_plants_type<- sfn_sites_type %>% 
  dplyr::select(si_code,typef) %>% 
  full_join( sfn_allplants_tax) %>% 
  dplyr::select(si_code,pl_sens_meth,typef,group) %>% 
  group_by(pl_sens_meth,typef,group) %>% tally()

sfn_plants_type %>% 
  group_by(typef) %>% 
  summarise(n_trees=sum(n))

# Methods per species
sfn_allplants_tax %>% 
  group_by(pl_species,pl_sens_meth) %>% tally() ->methods_species

# Installation height
sfn_allplants_tax ->installation_height
  

# 4. Percentage basal area ------------------------------------------------

dataset_trees_sp <- sfn_sites_type %>%
  right_join(sfn_allstands) %>% 
  right_join(
    sfn_sitespecies_tax %>%
      group_by(si_code) %>%
      mutate(total_ntrees = sum(sp_ntrees,na.rm=TRUE),
             nspecies = n_distinct(sp_name),
             # NOTE: This sum will be zero (0) when all are NAs
             percab_measured = sum(sp_basal_area_perc,na.rm=TRUE)),
    by='si_code'
  ) %>% 
    mutate(
      typeplant=ifelse(str_detect(type,'plant'),'plant',NA))



dataset_trees_sp %>% 
  dplyr::filter(percab_measured>100) %>% 
  dplyr::distinct(si_code,.keep_all = TRUE) %>% 
  dplyr::select(si_code,contains('contr'))
  




# 5. Treatments -----------------------------------------------------

sfn_allstands %>% 
  distinct(st_treatment) 
  
sfn_allplants_tax %>% 
  left_join(dplyr::select(sfn_allstands,si_code,st_treatment),by='si_code') %>% 
  filter(is.na(st_treatment) & !is.na(pl_treatment)) %>% 
  dplyr::select(si_code,st_treatment,pl_treatment) 


# FLUXNET

sfn_allsites %>% 
  pull(si_flux_network) %>% summary

# Dendroglobal

sfn_allsites %>% 
  pull(si_dendro_network) %>% summary

# 6. Save --------------------------------------------------------------------

save.image('sfn_datapaper_data.RData')

