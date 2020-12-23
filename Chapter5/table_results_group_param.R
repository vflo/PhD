table_results_group_param <- function(data, 
                                      parameters,
                                      weight){
  
  require(xtable)
  require(tidyverse)
  require(rlang)
  
df <- parameters %>% as.list()
  
  
  df %>% purrr::map(function(x){
    
      param <- x

      data1 <- data %>%
        mutate(parame = .data[[param]],
               n_weight = .data[[weight]]) %>% 
        filter(!is.na(parame))
      wtd.t.test(x=data1$parame[data_hydra2$group=="Angiosperms"],
                 y=data1$parame[data_hydra2$group=="Gymnosperms"],
                 weight=data1$n_weight[data_hydra2$group=="Angiosperms"],
                 weighty=data1$n_weight[data_hydra2$group=="Gymnosperms"],
                 samedata=FALSE)-> res_m
      # lm(parame ~ group, data = data1, weight=n_weight) -> res_m
      # res_m%>% 
      #   summary() -> res_mod
      #   anova(res_m) -> res_anova
      
        sig_cor <- symnum(res_m$coefficients[3],corr = FALSE,na = FALSE, 
                          cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                          symbols = c("***", "**", "*", ".", " "))
        
        df_res <- tibble(Parameter = param,
                         # n_Angiosperms = data %>% filter(group == "Angiosperms") %>% nrow(),
                         # n_Gymnosperms = data %>% filter(group == "Gymnosperms") %>% nrow(),
                         Angiosperms = round(res_m$additional[2],3),
                         Gymnosperms = round(res_m$additional[3],3),
                         AOV = sig_cor[1] %>% as.character()
                         )
        
      # sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
      #                   cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
      #                   symbols = c("***", "**", "*", ".", " "))
      
      # df_res <- tibble(Parameter = param,
      #                  # n_Angiosperms = data %>% filter(group == "Angiosperms") %>% nrow(),
      #                  # n_Gymnosperms = data %>% filter(group == "Gymnosperms") %>% nrow(),
      #                  Angiosperms = round(res_mod$coefficients[1,1],3),
      #                  Gymnosperms = round(res_mod$coefficients[2,1],3)+
      #                                   round(res_mod$coefficients[1,1],3),
      #                  AOV = sig_cor[2] %>% as.character(),
      #                  R2 = res_mod$adj.r.squared %>% round(3)
      # )

    
  }) -> df_res
  
  res <- df_res %>% bind_rows() 
  return(res)
}
