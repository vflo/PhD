require(ggplot2)
require(lmerTest)
require(lme4)
require(nlme)
require(car)
library(plyr)
require(dplyr)
require(lattice)
require(sjPlot)
library(predictmeans)
require(DHARMa)
library(sjmisc)
library(MuMIn)

# library(stj.lmer)
# sjp.setTheme(base = theme_bw(),axis.textsize = .75)


Data <- read.csv ("~/PHD/Chapter1/Methods meta-analysis V6.csv",sep = ";")

########################
### DATA PREPARATION ###
########################

aggregated.methods <- revalue(Data$method, c("compensation_heat_pulse" = 'Pulse', 
                                             "heat_field_deformation" = 'Field', 
                                             "heat_ratio" = 'Pulse', 
                                             "stem_heat_balance" = 'Balance',
                                             "T-max" = 'Pulse', 
                                             "thermal_dissipation" = 'Dissipation', 
                                             "transient_thermal_dissipation" = 'Dissipation'))
Data$method <- revalue(Data$method, c("compensation_heat_pulse" = 'CHP', 
                                      "heat_field_deformation" = 'HFD', 
                                      "heat_ratio" = 'HR', 
                                      "stem_heat_balance" = 'SHB',
                                      "T-max" = 'T-max', 
                                      "thermal_dissipation" = 'TD', 
                                      "transient_thermal_dissipation" = 'TTD'))

var <- ifelse(Data$Clearwater.correction=="yes", 1, 0)+ ifelse(Data$wound.correction=="yes", 1, 0)
var <- ifelse(var == "1", "corrected", "uncorrected")
new_methods <- paste(Data$method, var, sep = ".")

#Different agrupations of calibrations
calibrations_pairs <- paste(Data$study,Data$calibration,sep="_")
calibrations_calibrated <- paste(calibrations_pairs,Data$calibrated,sep="_")
calibrations_corrected <- paste(calibrations_pairs,Data$correction,sep='_')
calibrations <- paste(calibrations_calibrated,var,sep="_")

#Complete raw data creation
subdata <- cbind(Data,aggregated.methods,calibrations,calibrations_pairs,calibrations_calibrated,
                 calibrations_corrected,correction=var,new_methods)

subdata <- rename(subdata,TC=plant.material.2,MF=aggregated.methods,RC=radial.correction)
subdata$TC <- revalue(subdata$TC, c("whole plant without roots"="without roots"))

#porosity names
from <- c('A','C','E','F','monocot')
to <- c('Ring porous','Semi-diffuse porous','Diffuse porous','Tracheids','monocot')
gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- gsub(pattern[i], replacement[i], x, ...)
  x
}
subdata$porosity <- gsub2(from, to, subdata$porosity)

#### !!!!!Atention if we are using density (actual.density ...) or total flow (actual.totalflow ...)!!!!!
# subdata$actual.density <- subdata$actual.density
# subdata$measured.density <- subdata$measured.density
subdata$actual.density[subdata$actual.density <= 0] <- NA
subdata$measured.density[subdata$measured.density <= 0 ] <- NA
subdata <- subdata[!is.na(subdata$actual.density),]
subdata <- subdata[!is.na(subdata$measured.density),]
# subdata <- subdata[which(subdata$measuring.type == 'Sapflux density'),]

real <- (subdata$actual.density)
measured <- (subdata$measured.density)

subdata <- cbind(subdata,real,measured)

##### Meta-analyses data preparation #####
#log-log slope
Bar <- subdata %>% group_by(calibrations) %>%
  do(broom::tidy(lm(log(measured)~ log(real), data = .))) %>%
  select(calibrations, term, estimate) %>%
  tidyr::spread(term, estimate) %>%
  right_join(select(subdata,'study','specie','measuring.type','wood_density_in_de_article','density','correc.aplication.2', 'year', 'DBH','vessel_area', 'vessel_density',
                    'porosity','MF','RC','TC','time.resolution', 'calibrated', 'calibrations','calibrations_pairs',
                    'calibrations_calibrated','calibrations_corrected','method','correction','new_methods','CORRECT_APLICATION'), 
             by = 'calibrations') %>%
  unique()
Bar<-rename(Bar,log_slope=`log(real)`,log_intercept=`(Intercept)`)
#slope
Bar <- subdata %>% group_by(calibrations) %>%
  do(broom::tidy(lm(measured ~ real, data = .))) %>%
  select(calibrations, term, estimate) %>%
  tidyr::spread(term, estimate) %>%
  right_join(Bar,by = 'calibrations')
Bar<-rename(Bar,slope=real,intercept=`(Intercept)`)
#logratio
Bar <- subdata %>% group_by(calibrations) %>%
  mutate(logratio = log(measured/real)) %>%
  dplyr::summarize(log_ratio = mean(logratio))%>%
  right_join(Bar,by = 'calibrations')
#ratio
Bar <- subdata %>% group_by(calibrations) %>%
  mutate(Ratio = measured/real) %>%
  dplyr::summarize(ratio = mean(Ratio))%>%
  right_join(Bar,by = 'calibrations')
#z-cor
Bar <- subdata %>% group_by(calibrations) %>%
  dplyr::summarize(z = 0.5*log( (1+cor(measured,real,method="pearson")) / (1-cor(measured,real,method="pearson")) )) %>%
  right_join(Bar,by = 'calibrations')
#cor
Bar <- subdata %>% group_by(calibrations) %>%
  dplyr::summarize(cor = cor(measured,real,method="pearson")) %>%
  right_join(Bar,by = 'calibrations')
#range
Bar <- subdata %>% group_by(calibrations) %>%
  dplyr::summarize(range_mean = (min(real) + max(real))/2,min=min(real),max=max(real),median=median(real)) %>%
  right_join(Bar,by = 'calibrations')
# #RRMSE
# Bar <- subdata %>% group_by(calibrations) %>%
#   dplyr::mutate(sum = ((real-measured)/real)^2) %>%
#   dplyr::summarize(RRMSE=sqrt(sum(sum)/n()))%>%    
#   right_join(Bar,by = 'calibrations')

#RMSE
Bar <- subdata %>% group_by(calibrations) %>%
  do(broom::tidy(sjstats::rmse(lm(measured~ real, data = .)))) %>% summarize(rmse=x)%>%
  right_join(Bar,by = 'calibrations')

#NRMSE
Bar <- subdata %>% group_by(calibrations) %>%
  dplyr::mutate(dife1=((real-measured)^2))%>%
  dplyr::summarize(RMSE=((sqrt(sum(dife1)/n()))))%>%    
  right_join(Bar,by = 'calibrations')%>%
  mutate(NRMSE=RMSE/range_mean*100)

#RRMSE
Bar <- subdata %>% group_by(calibrations) %>%
  dplyr::mutate(dife2=(((real-measured)/real)^2))%>%
  dplyr::summarize(RRMSE=sqrt(sum(dife2)/n()))%>%    
  right_join(Bar,by = 'calibrations')

#Modelated error
Bar <- subdata %>% group_by(calibrations) %>%
  dplyr::mutate(rel_error=sqrt(((measured-real)/real)^2))%>%
  do(broom::tidy(lm( rel_error*100~ real, data = .))) %>%
  select(calibrations, term, estimate) %>%
  tidyr::spread(term, estimate) %>%
  dplyr::rename(slope_error=real,intercept_error=`(Intercept)`)%>%
  right_join(Bar,by = 'calibrations')

subdata <- subdata %>%
  dplyr::mutate(rel_error=sqrt(((measured-real)/real)^2))

# 
# 
# # Bar <- cbind(Bar,ratio = log(Bar$mean_calibration_measured/Bar$mean_calibration_real))
# 
# stderror.calib <- subdata %>% group_by(calibrations) %>%
#   do(broom::tidy(lm(measured ~ real, data = .))) %>%
#   select(calibrations, term, std.error) %>%
#   tidyr::spread(term, std.error) %>%
#   right_join(subdata[,c('study','specie','density','correc.aplication.2', 'year', 'DBH','vessel_area', 'vessel_density',
#                         'porosity', 'dens_cat','MF','RC','TC','time.resolution', 'calibrated', 'calibrations','calibrations_pairs',
#                         'calibrations_calibrated','calibrations_corrected','method','correction','new_methods','CORRECT_APLICATION')], 
#              by = 'calibrations') %>%
#   unique() %>%
#   data.frame()
# 
# slope_variance <- subdata %>% group_by(calibrations) %>%
#   do(broom::tidy((vcov(lm(measured ~ real, data = .))[2,2]))) %>%
#   right_join(subdata[,c('study','specie','density','correc.aplication.2', 'year', 'DBH','vessel_area', 'vessel_density',
#                         'porosity', 'dens_cat','MF','RC','TC','time.resolution', 'calibrated', 'calibrations','calibrations_pairs',
#                         'calibrations_calibrated','calibrations_corrected','method','correction','new_methods','CORRECT_APLICATION')], 
#              by = 'calibrations') %>%
#   unique() %>%
#   data.frame()
# 
# variability_ratio<- subdata %>% group_by(calibrations) %>%
#   mutate(logratio = log(measured/real)) %>%
#   summarize(ratio.var = ((sd(logratio))^2)) %>%
#   right_join(Bar[,c('study','specie','density','correc.aplication.2', 'year', 'DBH','vessel_area', 'vessel_density',
#                     'porosity', 'dens_cat','MF','RC','TC','time.resolution', 'calibrated', 'calibrations','calibrations_pairs',
#                     'calibrations_calibrated','calibrations_corrected','method','correction','new_methods','CORRECT_APLICATION')], 
#              by = 'calibrations') %>%
#   unique() %>%
#   data.frame()
# 
# n.calib <- subdata %>% group_by(calibrations) %>%
#   summarize(n = n()) %>%
#   right_join(Bar[,c('study','specie','density','correc.aplication.2', 'year', 'DBH','vessel_area', 'vessel_density',
#                     'porosity', 'dens_cat','MF','RC','TC','time.resolution', 'calibrated', 'calibrations','calibrations_pairs',
#                     'calibrations_calibrated','calibrations_corrected','method','correction','new_methods','CORRECT_APLICATION')], 
#              by = 'calibrations') %>%
#   unique() %>%
#   data.frame()
# # weight_mean <- data.frame(weight = 1/((stderror.calib.intercept$intercept.sd)^2))%>%
# #   cbind(MF = Bar$MF)%>%
# #   group_by(MF)%>%mutate(sum = sum(weight))
# # weight_mean <- weight_mean$weight/weight_mean$sum
# # 
# # 
# # weight_int = 1/((stderror.calib$X.Intercept.*sqrt(n.calib$n))^2)
# # weight_slope = 1/((stderror.calib$log_real*sqrt(n.calib$n))^2)