table_results_traits_param_MAP_height <- function(data, 
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
      lm(log_par ~ log_trait + MAP + pl_height , data = data%>% 
           dplyr::select(log_par,log_trait,MAP,pl_height,n_weight)%>% drop_na(), weight=n_weight) -> res_m
      stats::step(res_m,trace=0)->res_m
      res_m%>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      res_mod$coefficients[,1] <-  paste(round(res_mod$coefficients[,1],3),sig_cor)
      names_var<- rownames(res_mod$coefficients)
      coef_vector <- c(any(names_var == "(Intercept)"),
                       any(names_var == "log_trait"),
                       any(names_var == "MAP"),
                       any(names_var == "pl_height"))
      
      coefficients <- res_mod$coefficients %>% t()
      
      df_res <- tibble(Parameter = paste0("ln(",param,")"),
                       Trait = paste0("ln(",trai,")"),
                       n_species = nobs(res_m),
                       Intercept = ifelse(coef_vector[1],coefficients[1,"(Intercept)"],"NI"),
                       trait = ifelse(coef_vector[2],coefficients[1,"log_trait"],"NI"),
                       MAP = ifelse(coef_vector[3],coefficients[1,"MAP"],"NI"),
                       pl_height = ifelse(coef_vector[4],coefficients[1,"pl_height"],"NI"),
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
      lm(parame ~ log_trait + MAP + pl_height, data = data %>% 
           dplyr::select(parame,log_trait,MAP,pl_height,n_weight)%>% drop_na(), weight=n_weight) -> res_m
      stats::step(res_m,trace=0)->res_m
      res_m%>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      res_mod$coefficients[,1] <-  paste(round(res_mod$coefficients[,1],3),sig_cor)
      names_var<- rownames(res_mod$coefficients)
      coef_vector <- c(any(names_var == "(Intercept)"),
                       any(names_var == "log_trait"),
                       any(names_var == "MAP"),
                       any(names_var == "pl_height"))
      
      coefficients <- res_mod$coefficients %>% t()
      
      df_res <- tibble(Parameter = param,
                       Trait = paste0("ln(",trai,")"),
                       n_species = nobs(res_m),
                       Intercept = ifelse(coef_vector[1],coefficients[1,"(Intercept)"],"NI"),
                       trait = ifelse(coef_vector[2],coefficients[1,"log_trait"],"NI"),
                       MAP = ifelse(coef_vector[3],coefficients[1,"MAP"],"NI"),
                       pl_height = ifelse(coef_vector[4],coefficients[1,"pl_height"],"NI"),
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
      lm(log_par ~ trait + MAP + pl_height, data = data%>% 
           dplyr::select(log_par,trait,MAP,pl_height,n_weight)%>% drop_na(), weight=n_weight) -> res_m
      stats::step(res_m,trace=0)->res_m
      res_m %>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      res_mod$coefficients[,1] <-  paste(round(res_mod$coefficients[,1],3),sig_cor)
      names_var<- rownames(res_mod$coefficients)
      coef_vector <- c(any(names_var == "(Intercept)"),
                       any(names_var == "trait"),
                       any(names_var == "MAP"),
                       any(names_var == "pl_height"))
      
      coefficients <- res_mod$coefficients %>% t()
      
      df_res <- tibble(Parameter = paste0("ln(",param,")"),
                       Trait = trai,
                       n_species = nobs(res_m),
                       Intercept = ifelse(coef_vector[1],coefficients[1,"(Intercept)"],"NI"),
                       trait = ifelse(coef_vector[2],coefficients[1,"trait"],"NI"),
                       MAP = ifelse(coef_vector[3],coefficients[1,"MAP"],"NI"),
                       pl_height = ifelse(coef_vector[4],coefficients[1,"pl_height"],"NI"),
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
      lm(parame ~ trait + MAP + pl_height, data = data%>% 
           dplyr::select(parame,trait,MAP,pl_height,n_weight)%>% drop_na(), weight=n_weight) -> res_m
      stats::step(res_m,trace=0)->res_m
      res_m%>% 
        summary() -> res_mod
      
      sig_cor <- symnum(res_mod$coefficients[,4],corr = FALSE,na = FALSE, 
                        cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), 
                        symbols = c("***", "**", "*", ".", " "))
      
      res_mod$coefficients[,1] <-  paste(round(res_mod$coefficients[,1],3),sig_cor)
      names_var<- rownames(res_mod$coefficients)
      coef_vector <- c(any(names_var == "(Intercept)"),
                       any(names_var == "trait"),
                       any(names_var == "MAP"),
                       any(names_var == "pl_height"))
      
      
      coefficients <- res_mod$coefficients %>% t()
      
      df_res <- tibble(Parameter = param,
                       Trait = trai,
                       n_species = nobs(res_m),
                       Intercept = ifelse(coef_vector[1],coefficients[1,"(Intercept)"],"NI"),
                       trait = ifelse(coef_vector[2],coefficients[1,"trait"],"NI"),
                       MAP = ifelse(coef_vector[3],coefficients[1,"MAP"],"NI"),
                       pl_height = ifelse(coef_vector[4],coefficients[1,"pl_height"],"NI"),
                       R2 = res_mod$adj.r.squared %>% round(3)
      )
    }    
    
    return(df_res)
    
  }) -> df_res
  
  res <- df_res %>% bind_rows() 
  return(res)
}