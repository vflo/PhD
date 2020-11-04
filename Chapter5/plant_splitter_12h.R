path <- "~/Chapter5/data/site_species_daylight/"
site_names <- list.files(path = path)


purrr::map(site_names,function(x){

  load(paste0(path,x))
  
  print(x)

  x2 <- x2 %>% as_tibble() %>% split(.[['pl_code']],drop = TRUE)
  
  purrr::map(x2,function(res){

      res[which(res$sub_zero == "YES"),] %>% unique()->foo
      if(nrow(foo)>0){
        df <- tibble(date_ymd = foo$date_ymd,
                     season_add = foo$season_add)
        pre_sub_zero <- df %>%
          group_by(season_add) %>%
          filter(grepl('Winter', season_add))
        if(nrow(pre_sub_zero) > 0){
          pre_sub_zero <- pre_sub_zero%>%
            summarise(min_date = min(date_ymd),
                      max_date = max(date_ymd))%>%
            split(seq(nrow(.))) %>%
            map(function(x){
              seq(x$'min_date'-15*24*60*60,
                  x$'max_date'+35*24*60*60,
                  by="hour")
            }) %>%
            map(data_frame) %>%
            bind_rows()
        }else{pre_sub_zero = NA}}else{pre_sub_zero = NA}
    
      if(any(!is.na(pre_sub_zero))){
        res %>%
          filter(!as.numeric(res$TIMESTAMP) %in%
                   as.numeric(unlist(pre_sub_zero[,1]))) -> x2
      }else{x2 <- res}

     }
    ) 
  
  x2 <- x2%>%
    bind_rows()

  return(x2)
  
  }
)%>% 
  bind_rows()%>%
  mutate(pl_species = fct_recode(pl_species,
                                 "Vouacapoua americana"="Vacapoua americana",
                                 'Larix sibirica '='Larix sibirica Ledeb.',
                                 "Betula pubescens" = "Betula pubescens subsp. czerepanovii"),
         pl_species = as.character(pl_species)) %>%
  split(.[['pl_species']],drop = TRUE) ->res1


purrr::map(res1,function(x2){
  
  save(x2, file=paste0('~/Chapter5/data/species_data/',
                       unique(x2$pl_species), '.RData'))
}
)
