plot_traits <- function(data,trait,parameter,parameter_se,weight,log_trait = FALSE,log_par = FALSE,sig=TRUE){
  require(tidyverse)
  require(ggplot2)
  theme_set(theme_bw())
  theme_update(panel.grid = element_blank())
  require(viridis)
  require(rlang)
  
  
  if(log_trait & log_par){
    
    lm(data[[parameter]] %>% log()~data[[trait]] %>% log(),weights = data[['n']]) %>% anova()->faaa

  gg <- data %>% filter(!is.na(.data[[trait]]), !is.na(.data[[parameter]])) %>% 
    ggplot(aes(x=.data[[trait]] %>% log(), y=.data[[parameter]] %>% log(), weight = .data[[weight]]))+
    geom_point(aes(size = .data[[weight]]), alpha = 0.5, show.legend = FALSE)+
    geom_smooth(method= 'lm', formula = y ~ x)+
    NULL
  
  gg_data <- ggplot_build(gg)
  head(gg_data$data[[2]])
  gg1 <- gg_data$data[[2]]

  
  gg2 <-ggplot()+
    geom_errorbar(data = data %>% filter(!is.na(.data[[trait]]),
                                         !is.na(.data[[parameter]])),
                  aes(x = .data[[trait]] %>% log(),
                      ymin = (.data[[parameter]] - .data[[parameter_se]]) %>% log(),
                      ymax = (.data[[parameter]] + .data[[parameter_se]]) %>% log()),
                  color = "grey60")+
    geom_point(data = data %>% filter(!is.na(.data[[trait]]),
                                      !is.na(.data[[parameter]])),
               aes(x=.data[[trait]] %>% log(),
               y=.data[[parameter]] %>% log(),
               size = .data[[weight]],
               color = group), alpha = 0.5, show.legend = FALSE)+
    # scale_color_manual(values = c("#bb3200","#00AFBB"))+
    scale_color_npg()+
    scale_size(range = c(2,8))+
    NULL
  
    if(sig){
      
      if(length(faaa$`Pr(>F)`)>1 & faaa$`Pr(>F)`[1] <0.05){
      gg2 <- gg2+
              geom_line(data = gg1,aes(x=x,y=y),color="grey20",size = 1.5)+
              geom_line(data = gg1,aes(x=x,y=ymin),color="grey20",linetype=3)+
              geom_line(data = gg1,aes(x=x,y=ymax),color="grey20",linetype=3)
      }
    }
  }
  

  if(log_trait & !log_par){
    
    lm(data[[parameter]] ~data[[trait]] %>% log()*data[['group']],weights = data[['n']]) %>% anova()->faaa
    feee<-faaa
    if(faaa$`Pr(>F)`[3] >0.05){lm(data[[parameter]] ~data[[trait]] %>% log()+
                                   data[['group']],weights = data[['n']]) %>% anova()->feee}
    if(faaa$`Pr(>F)`[3] >0.05 &
       feee$`Pr(>F)`[2] >0.05){lm(data[[parameter]] ~data[[trait]] %>% log(),weights = data[['n']]) %>% anova()->feee}
    faaa<-feee
    
    gg <- data %>% filter(!is.na(.data[[trait]]), !is.na(.data[[parameter]])) %>% 
      ggplot(aes(x=.data[[trait]] %>% log(), y=.data[[parameter]], weight = .data[[weight]]))+
      geom_point(aes(size = .data[[weight]]), alpha = 0.5, show.legend = FALSE)+
      geom_smooth(method= 'lm', formula = y ~ x)+
      geom_smooth(aes(group=group),method= 'lm', formula = y ~ x)+
      NULL
    
    gg_data <- ggplot_build(gg)
    head(gg_data$data[[2]])
    gg1 <- gg_data$data[[2]]
    gg3 <- gg_data$data[[3]] %>% split(.[['group']],drop = TRUE)
    
    gg2 <-ggplot()+
      geom_errorbar(data = data %>% filter(!is.na(.data[[trait]]),
                                           !is.na(.data[[parameter]])),
                    aes(x = .data[[trait]] %>% log(),
                        ymin = (.data[[parameter]] - .data[[parameter_se]]),
                        ymax = (.data[[parameter]] + .data[[parameter_se]])),
                    color = "grey60")+
      geom_point(data = data %>% filter(!is.na(.data[[trait]]),
                                        !is.na(.data[[parameter]])),
                 aes(x=.data[[trait]] %>% log(),
                     y=.data[[parameter]],
                     size = .data[[weight]],
                     color = group), alpha = 0.5, show.legend = FALSE)+
      # scale_color_manual(values = c("#bb3200","#00AFBB"))+
      scale_color_npg()+
      scale_size(range = c(2,8))+
      NULL
    
    if(sig){
      if(length(faaa$`Pr(>F)`)>1 & faaa$`Pr(>F)`[1] <0.05){
        gg2 <- gg2+
          geom_line(data = gg1,aes(x=x,y=y),color="grey20",size = 1.5)+
          geom_line(data = gg1,aes(x=x,y=ymin),color="grey20",linetype=3)+
          geom_line(data = gg1,aes(x=x,y=ymax),color="grey20",linetype=3)
      }
    }
  }
  
  if(!log_trait & log_par){
    
    lm(data[[parameter]] %>% log()~data[[trait]],weights = data[['n']]) %>% anova()->faaa

    gg <- data %>% filter(!is.na(.data[[trait]]), !is.na(.data[[parameter]])) %>% 
      ggplot(aes(x=.data[[trait]], y=.data[[parameter]] %>% log(), weight = .data[[weight]]))+
      geom_point(aes(size = .data[[weight]]), alpha = 0.5, show.legend = FALSE)+
      geom_smooth(method= 'lm', formula = y ~ x)+
      NULL
    
    gg_data <- ggplot_build(gg)
    head(gg_data$data[[2]])
    gg1 <- gg_data$data[[2]]

    gg2 <-ggplot()+
      geom_errorbar(data = data %>% filter(!is.na(.data[[trait]]),
                                           !is.na(.data[[parameter]])),
                    aes(x = .data[[trait]],
                        ymin = (.data[[parameter]] - .data[[parameter_se]]) %>% log(),
                        ymax = (.data[[parameter]] + .data[[parameter_se]]) %>% log()),
                    color = "grey60")+
      geom_point(data = data %>% filter(!is.na(.data[[trait]]),
                                        !is.na(.data[[parameter]])),
                 aes(x=.data[[trait]],
                     y=.data[[parameter]] %>% log(),
                     size = .data[[weight]],
                     color = group), alpha = 0.5, show.legend = FALSE)+
      # scale_color_manual(values = c("#bb3200","#00AFBB"))+
      scale_color_npg()+
      scale_size(range = c(2,8))+
      NULL
    
    if(sig){
      
      if(length(faaa$`Pr(>F)`)>1 & faaa$`Pr(>F)`[1] <0.05){
        gg2 <- gg2+
          geom_line(data = gg1,aes(x=x,y=y),color="grey20",size = 1.5)+
          geom_line(data = gg1,aes(x=x,y=ymin),color="grey20",linetype=3)+
          geom_line(data = gg1,aes(x=x,y=ymax),color="grey20",linetype=3)
      }
    }
  }
  
  if(!log_trait & !log_par){
    
    lm(data[[parameter]]~data[[trait]],weights = data[['n']]) %>% anova()->faaa
    
    gg <- data %>% filter(!is.na(.data[[trait]]), !is.na(.data[[parameter]])) %>% 
      ggplot(aes(x=.data[[trait]], y=.data[[parameter]], weight = .data[[weight]]))+
      geom_point(aes(size = .data[[weight]]), alpha = 0.5, show.legend = FALSE)+
      geom_smooth(method= 'lm', formula = y ~ x)+
      geom_smooth(aes(group=group),method= 'lm', formula = y ~ x)+
      NULL
    
    gg_data <- ggplot_build(gg)
    head(gg_data$data[[2]])
    gg1 <- gg_data$data[[2]]

    gg2 <-ggplot()+
      geom_errorbar(data = data %>% filter(!is.na(.data[[trait]]),
                                           !is.na(.data[[parameter]])),
                    aes(x = .data[[trait]],
                        ymin = (.data[[parameter]] - .data[[parameter_se]]),
                        ymax = (.data[[parameter]] + .data[[parameter_se]])),
                    color = "grey60")+
      geom_point(data = data %>% filter(!is.na(.data[[trait]]),
                                        !is.na(.data[[parameter]])),
                 aes(x=.data[[trait]],
                     y=.data[[parameter]],
                     size = .data[[weight]],
                     color = group), alpha = 0.5, show.legend = FALSE)+
      # scale_color_manual(values = c("#bb3200","#00AFBB"))+
      scale_color_npg()+
      scale_size(range = c(2,8))+
      NULL
    
    if(sig){
      
      if(length(faaa$`Pr(>F)`)>1 & faaa$`Pr(>F)`[1] <0.05){
        gg2 <- gg2+
          geom_line(data = gg1,aes(x=x,y=y),color="grey20",size = 1.5)+
          geom_line(data = gg1,aes(x=x,y=ymin),color="grey20",linetype=3)+
          geom_line(data = gg1,aes(x=x,y=ymax),color="grey20",linetype=3)
      }
    }
  }
  
  return(gg2)
  
  
}