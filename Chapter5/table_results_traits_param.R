table_results_traits_param <- function(data, 
                                       parameters, 
                                       traits, 
                                       parameters_log, 
                                       traits_log,
                                       weight){
  
  require(xtable)
  require(tidyverse)
  require(rlang)
  
  df <- expand.grid(traits = traits,
                    parameters = parameters,
                    stringsAsFactors = FALSE) %>% 
    cbind(expand.grid(traits_log = traits_log,
                      parameters_log = parameters_log,
                      stringsAsFactors = FALSE))
  df <- lapply(as.list(1:dim(df)[1]), function(x) df[x[1],])
  
  
  df %>% purrr::map(function(x){
    
    if((x[['parameters_log']] & x[["traits_log"]])){
      
      param <- x$parameters
      trai <- x$traits
      
      data <- data %>%
        mutate(log_par = log(.data[[param]]),
               log_trait = log(.data[[trai]]),
               n_weight = .data[[weight]])
      lm(log_par ~ log_trait, data = data, weight=n_weight) -> res_m
      res_m%>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      df_res <- tibble(Parameter = paste0("ln(",param,")"),
                       Trait = paste0("ln(",trai,")"),
                       n_species = nobs(res_m),
                       Intercept = paste(round(res_mod$coefficients[1,1],3),sig_cor[1]),
                       Slope = paste(round(res_mod$coefficients[2,1],3),sig_cor[2]),
                       R2 = res_mod$adj.r.squared %>% round(3)
      )
    }
    
    if((!x[['parameters_log']] & x[["traits_log"]])){
      
      param <- x$parameters
      trai <- x$traits
      
      data <- data %>%
        mutate(parame = (.data[[param]]),
               log_trait = log(.data[[trai]]),
               n_weight = .data[[weight]])
      lm(parame ~ log_trait, data = data, weight=n_weight) -> res_m
      res_m%>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      df_res <- tibble(Parameter = param,
                       Trait = paste0("ln(",trai,")"),
                       n_species = nobs(res_m),
                       Intercept = paste(round(res_mod$coefficients[1,1],3),sig_cor[1]),
                       Slope = paste(round(res_mod$coefficients[2,1],3),sig_cor[2]),
                       R2 = res_mod$adj.r.squared %>% round(3)
      )      
    }
    
    if((x[['parameters_log']] & !x[["traits_log"]])){
      
      param <- x$parameters
      trai <- x$traits
      
      data <- data %>%
        mutate(log_par = log(.data[[param]]),
               trait = .data[[trai]],
               n_weight = .data[[weight]])
      lm(log_par ~ trait, data = data, weight=n_weight) -> res_m
      res_m %>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      df_res <- tibble(Parameter = paste0("ln(",param,")"),
                       Trait = trai,
                       n_species = nobs(res_m),
                       Intercept = paste(round(res_mod$coefficients[1,1],3),sig_cor[1]),
                       Slope = paste(round(res_mod$coefficients[2,1],3),sig_cor[2]),
                       R2 = res_mod$adj.r.squared %>% round(3)
      )
    }
    
    if((!x[['parameters_log']] & !x[["traits_log"]])){
      
      param <- x$parameters
      trai <- x$traits
      
      data <- data %>%
        mutate(parame = .data[[param]],
               trait = .data[[trai]],
               n_weight = .data[[weight]])
      lm(parame ~ trait, data = data, weight=n_weight) -> res_m
      res_m%>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      df_res <- tibble(Parameter = param,
                       Trait = trai,
                       n_species = nobs(res_m),
                       Intercept = paste(round(res_mod$coefficients[1,1],3),sig_cor[1]),
                       Slope = paste(round(res_mod$coefficients[2,1],3),sig_cor[2]),
                       R2 = res_mod$adj.r.squared %>% round(3)
      )
    }    
    
    return(df_res)
    
  }) -> df_res
  
  res <- df_res %>% bind_rows() 
  return(res)
}