#### FAMILIES ####


source('data preparation 1-1.R')
subdata1 <- filter(Bar,calibrated=="no",CORRECT_APLICATION == "yes")

#resume table
 resume_table <- subdata1 %>% dplyr::select(study,method,specie,porosity,TC,DBH) %>% 
  group_by(study,method,specie,TC) %>% 
  dplyr::summarise(porosity = unique(porosity),
                   DBH = mean(DBH))
write.table(resume_table,file='resume_table.csv',dec='.',sep=';')

#Clearwater and not clearwater aplication on TD calibrations. Is there any difference?
TD_data <- subdata1 %>% filter(method == 'TD')
foo_slope <- lmer(slope  ~ Clearwater.correction + TC + (1|study) + (1|specie),data=TD_data)
foo_log_slope <- lmer(log_slope  ~ Clearwater.correction + TC + (1|study) + (1|specie),data=TD_data)
foo_log_ratio <- lmer(log_ratio~ Clearwater.correction + TC + (1|study) + (1|specie),data=TD_data)
foo_z <- lmer(z~ Clearwater.correction + TC + (1|study) + (1|specie),data=TD_data)

summary(foo_slope)
summary(foo_log_slope)
summary(foo_log_ratio)
summary(foo_z)

anova(foo_slope)
anova(foo_log_slope)
anova(foo_log_ratio)
anova(foo_z)


# models
foo_slope <- lmer(slope  ~ (MF + TC) + (1|study) + (1|specie),data=subdata1)
foo_log_slope <- lmer(log_slope  ~ (MF + TC) + (1|study) + (1|specie),data=subdata1)
foo_log_ratio <- lmer(log_ratio~ (MF+TC) + (1|study) + (1|specie),data=subdata1)
foo_z <- lmer(z~ (MF+TC) + (1|study) + (1|specie),data=subdata1)
# foo_intercept <- lmer(intercept  ~ MF  + (1|study) + (1|specie),data=subdata1)

r.squaredGLMM(foo_slope)
r.squaredGLMM(foo_log_slope)
r.squaredGLMM(foo_log_ratio)
r.squaredGLMM(foo_z)

# lsmeans
aux1<-emmeans::emmeans(foo_slope,pairwise~MF, adjust="tukey")
lsslopeMF <- multcomp::cld(aux1$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux2<-emmeans::emmeans(foo_slope,pairwise~TC, adjust="tukey")
lsslopeTC <- multcomp::cld(aux2$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

aux5<-emmeans::emmeans(foo_z,pairwise~MF, adjust="tukey")
lszMF <- multcomp::cld(aux5$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux6<-emmeans::emmeans(foo_z,pairwise~TC, adjust="tukey")
lszTC <- multcomp::cld(aux6$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

aux7<-emmeans::emmeans(foo_log_slope,pairwise~MF, adjust="tukey")
lslog_slopeMF <- multcomp::cld(aux7$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux8<-emmeans::emmeans(foo_log_slope,pairwise~TC, adjust="tukey")
lslog_slopeTC <- multcomp::cld(aux8$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

aux9<-emmeans::emmeans(foo_log_ratio,pairwise~MF, adjust="tukey")
lslog_ratioMF <- multcomp::cld(aux9$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux10<-emmeans::emmeans(foo_log_ratio,pairwise~TC, adjust="tukey")
lslog_ratioTC <- multcomp::cld(aux10$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

# aux11<-emmeans::emmeans(foo_intercept,pairwise~MF, adjust="tukey")
# lsinterceptMF <- multcomp::cld(aux11$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
# lsmeans::test(emmeans::emmeans(foo_intercept,~MF), adjust="tukey",null=0)

emmeans::test(emmeans::emmeans(foo_log_ratio,~MF), adjust="tukey",null=0)
emmeans::test(emmeans::emmeans(foo_slope,~MF), adjust="tukey",null=1)
emmeans::test(emmeans::emmeans(foo_log_slope,~MF), adjust="tukey",null=1)

#accuracy deviation %
(exp(summary(aux9)$emmeans[,'emmean'])-1)*100

#r from z-cor
tanh(summary(aux5)$emmeans[,'emmean'])

#n species
group_by(subdata1,MF)%>%summarize(n_distinct(specie))

# resume table
library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(foo_log_ratio,foo_slope,foo_log_slope,foo_z)

# plots
source('methodplot.R')
p1 <- MFplot(lslog_ratioMF,intercept=0,xlab = 'Method family',ylab = 'Ln-Ratio', ylimits= c(-0.8,0.4))
p2 <- TCplot(lslog_ratioTC,intercept=0,xlab = 'Calibration material',ylab = 'Ln-Ratio', ylimits= c(-0.8,0.4))

p3 <- MFplot(lsslopeMF,xlab = 'Method family',ylab = 'Slope', ylimits= c(0.4,1.3))
p4 <- TCplot(lsslopeTC,xlab = 'Calibration material',ylab = 'Slope', ylimits= c(0.4,1.3))

p5 <- MFplot(lslog_slopeMF,xlab = 'Method family',ylab = 'Slope (ln-ln)', ylimits= c(0.6,1.3))
p6 <- TCplot(lslog_slopeTC,xlab = 'Calibration material',ylab = 'Slope (ln-ln)', ylimits= c(0.6,1.3))

p7 <- MFplot(lszMF,intercept=-9999,xlab = 'Method family',ylab = 'Z-Cor', ylimits= c(1,3.5))
p8 <- TCplot(lszTC,intercept=-9999,xlab = 'Calibration material',ylab = 'Z-Cor', ylimits= c(1,3.5))

# p9 <- MFplot(lsinterceptMF,intercept=0,xlab = 'Method family',ylab = 'Intercept', ylimits= c(-4,6))

source("multiplot.R")
plotfam <- multiplot(p1,p3,p5,p7,p2,p4,p6,p8,byrow=TRUE,layout=matrix(c(1,2,3,4,5,6,7,8),nrow=4))



#### METHODS ####

source('data preparation 1-1.R')
# source('data preparation 1-1-density.R')
subdata1 <- filter(Bar,calibrated=="no",CORRECT_APLICATION == "yes")
subdata1$method <- factor(subdata1$method, levels (subdata1$method) [c(1,5,3,2,4,6,7)])
# subdata1$specie<-stringr::str_trunc(as.character(subdata1$specie),5,side='right',ellipsis = '')
# models
foo_slope <- lmer(slope  ~ (method + TC) + (1|study) + (1|specie),data=subdata1)
foo_log_slope <- lmer(log_slope  ~ (method + TC) + (1|study) + (1|specie),data=subdata1)
foo_log_ratio <- lmer(log_ratio~ (method+TC) + (1|study) + (1|specie),data=subdata1)
foo_z <- lmer(z~ (method+TC) + (1|study) + (1|specie),data=subdata1)
# foo_intercept <- lmer(intercept  ~ method  + (1|study) + (1|specie),data=subdata1)

r.squaredGLMM(foo_slope)
r.squaredGLMM(foo_log_slope)
r.squaredGLMM(foo_log_ratio)
r.squaredGLMM(foo_z)

# emmeans
aux1<-emmeans::emmeans(foo_slope,pairwise~method, adjust="tukey")
lsslopemethod <- multcomp::cld(aux1$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux2<-emmeans::emmeans(foo_slope,pairwise~TC, adjust="tukey")
lsslopeTC <- multcomp::cld(aux2$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

aux5<-emmeans::emmeans(foo_z,pairwise~method, adjust="tukey")
lszmethod <- multcomp::cld(aux5$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux6<-emmeans::emmeans(foo_z,pairwise~TC, adjust="tukey")
lszTC <- multcomp::cld(aux6$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

aux7<-emmeans::emmeans(foo_log_slope,pairwise~method, adjust="tukey")
lslog_slopemethod <- multcomp::cld(aux7$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux8<-emmeans::emmeans(foo_log_slope,pairwise~TC, adjust="tukey")
lslog_slopeTC <- multcomp::cld(aux8$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

aux9<-emmeans::emmeans(foo_log_ratio,pairwise~method, adjust="tukey")
lslog_ratiomethod <- multcomp::cld(aux9$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
aux10<-emmeans::emmeans(foo_log_ratio,pairwise~TC, adjust="tukey")
lslog_ratioTC <- multcomp::cld(aux10$emmeans,alpha=0.05,Letters=letters,adjust="tukey")

# aux11<-emmeans::emmeans(foo_intercept,pairwise~method, adjust="tukey")
# lsinterceptmethod <- multcomp::cld(aux11,alpha=0.05,Letters=letters,adjust="tukey")
# lsmeans::test(emmeans::emmeans(foo_intercept,~method), adjust="tukey",null=0)


emmeans::test(emmeans::emmeans(foo_log_ratio,~method), adjust="tukey",null=0)
emmeans::test(emmeans::emmeans(foo_slope,~method), adjust="tukey",null=1)
emmeans::test(emmeans::emmeans(foo_log_slope,~method), adjust="tukey",null=1)

#accuracy deviation %
(exp(summary(aux9)$emmeans[,'emmean'])-1)*100

#r from z-cor
tanh(summary(aux5)$emmeans[,'emmean'])

#n species
group_by(subdata1,method)%>%summarize(n_distinct(specie))

# resume table
tab_model(foo_log_ratio,foo_slope,foo_log_slope,foo_z)

# plots
source('methodplot.R')
p1 <- methodplot(lslog_ratiomethod,intercept=0,xlab = 'Method family',ylab = 'Ln-Ratio', ylimits= c(-1,0.6))
p2 <- TCplot(lslog_ratioTC,intercept=0,xlab = 'Calibration material',ylab = 'Ln-Ratio', ylimits= c(-1,0.6))

p3 <- methodplot(lsslopemethod,xlab = 'Method family',ylab = 'Slope', ylimits= c(0.2,1.3))
p4 <- TCplot(lsslopeTC,xlab = 'Calibration material',ylab = 'Slope', ylimits= c(0.2,1.3))

p5 <- methodplot(lslog_slopemethod,xlab = 'Method family',ylab = 'Slope (ln-ln)', ylimits= c(0.3,1.3))
p6 <- TCplot(lslog_slopeTC,xlab = 'Calibration material',ylab = 'Slope (ln-ln)', ylimits= c(0.3,1.3))

p7 <- methodplot(lszmethod,intercept=-9999,xlab = 'Method family',ylab = 'Z-Cor', ylimits= c(0.6,3.5))
p8 <- TCplot(lszTC,intercept=-9999,xlab = 'Calibration material',ylab = 'Z-Cor', ylimits= c(0.6,3.5))

# p9 <- methodplot(lsinterceptmethod,intercept=0,xlab = 'Method family',ylab = 'Intercept', ylimits= c(-5,17))

source("multiplot.R")
plotfam <- multiplot(p1,p3,p5,p7,p2,p4,p6,p8,byrow=TRUE,layout=matrix(c(1,2,3,4,5,6,7,8),nrow=4))




#### Resume studies correct aplication by methods ####
xxx <- subdata1 %>% group_by(study,CORRECT_APLICATION,method) %>% summarise() 
    # summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))
    
    
#### CORRECT APLICATION ####    
    source('data preparation 1-1.R')
    subdata_TD_CHP <- filter(Bar,calibrated=="no",method == "TD"|method =="CHP")
    
    foo_slope_TD_CHP <- lmer(slope  ~ (CORRECT_APLICATION*method + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdata_TD_CHP)
    foo_log_slope_TD_CHP <- lmer(log_slope  ~ (CORRECT_APLICATION*method + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdata_TD_CHP)
    foo_log_ratio_TD_CHP <- lmer(log_ratio~ (CORRECT_APLICATION*method + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdata_TD_CHP)
    foo_z_TD_CHP <- lmer(z~ (CORRECT_APLICATION*method + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdata_TD_CHP)
    
    aux1<-emmeans::emmeans(foo_slope_TD_CHP,pairwise~method*CORRECT_APLICATION, adjust="tukey")
    lsslopemethod_TD_CHP <- multcomp::cld(aux1$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
        # aux2<-emmeans::emmeans(foo_slope_TD_CHP,pairwise~TC, adjust="tukey")
        # lsslopeTC_TD_CHP <- multcomp::cld(aux2,alpha=0.05,Letters=letters,adjust="tukey")

    aux3<-emmeans::emmeans(foo_z_TD_CHP,pairwise~method*CORRECT_APLICATION, adjust="tukey")
    lszmethod_TD_CHP <- multcomp::cld(aux3$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
        # aux4<-emmeans::emmeans(foo_z_TD_CHP,pairwise~TC, adjust="tukey")
        # lszTC_TD_CHP <- multcomp::cld(aux4,alpha=0.05,Letters=letters,adjust="tukey")

    aux5<-emmeans::emmeans(foo_log_slope_TD_CHP,pairwise~method*CORRECT_APLICATION, adjust="tukey")
    lslog_slopemethod_TD_CHP <- multcomp::cld(aux5$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
        # aux6<-emmeans::emmeans(foo_log_slope_TD_CHP,pairwise~TC, adjust="tukey")
        # lslog_slopeTC_TD_CHP <- multcomp::cld(aux6,alpha=0.05,Letters=letters,adjust="tukey")

    aux7<-emmeans::emmeans(foo_log_ratio_TD_CHP,pairwise~method*CORRECT_APLICATION, adjust="tukey")
    lslog_ratiomethod_TD_CHP <- multcomp::cld(aux7$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
        # aux8<-emmeans::emmeans(foo_log_ratio_TD_CHP,pairwise~TC, adjust="tukey")
        # lslog_ratioTC_TD_CHP <- multcomp::cld(aux8,alpha=0.05,Letters=letters,adjust="tukey")
    
    # subdataTD <- filter(Bar,calibrated=="no",method == "TD")
    # subdataCHP <- filter(Bar,calibrated=="no",method == "CHP")
# # models TD
#     foo_slope_TD <- lmer(slope  ~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataTD)
#     foo_log_slope_TD <- lmer(log_slope  ~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataTD)
#     foo_log_ratio_TD <- lmer(log_ratio~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataTD)
#     foo_z_TD <- lmer(z~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataTD)
#     
# # models CHP
#     foo_slope_CHP <- lmer(slope  ~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataCHP)
#     foo_log_slope_CHP <- lmer(log_slope  ~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataCHP)
#     foo_log_ratio_CHP <- lmer(log_ratio~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataCHP)
#     foo_z_CHP <- lmer(z~ (CORRECT_APLICATION + TC) + (1|study/calibrations_pairs) + (1|specie),data=subdataCHP)
    
# # lsmeans TD
#     aux1<-emmeans::emmeans(foo_slope_TD,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lsslopemethod_TD <- multcomp::cld(aux1,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux2<-emmeans::emmeans(foo_slope_TD,pairwise~TC, adjust="tukey")
#     # lsslopeTC_TD <- multcomp::cld(aux2,alpha=0.05,Letters=letters,adjust="tukey")
#     
#     aux3<-emmeans::emmeans(foo_z_TD,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lszmethod_TD <- multcomp::cld(aux3,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux4<-emmeans::emmeans(foo_z_TD,pairwise~TC, adjust="tukey")
#     # lszTC_TD <- multcomp::cld(aux4,alpha=0.05,Letters=letters,adjust="tukey")
#     
#     aux5<-emmeans::emmeans(foo_log_slope_TD,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lslog_slopemethod_TD <- multcomp::cld(aux5,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux6<-emmeans::emmeans(foo_log_slope_TD,pairwise~TC, adjust="tukey")
#     # lslog_slopeTC_TD <- multcomp::cld(aux6,alpha=0.05,Letters=letters,adjust="tukey")
#     
#     aux7<-emmeans::emmeans(foo_log_ratio_TD,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lslog_ratiomethod_TD <- multcomp::cld(aux7,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux8<-emmeans::emmeans(foo_log_ratio_TD,pairwise~TC, adjust="tukey")
#     # lslog_ratioTC_TD <- multcomp::cld(aux8,alpha=0.05,Letters=letters,adjust="tukey")
#     
# # lsmeans CHP
#     aux9<-emmeans::emmeans(foo_slope_CHP,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lsslopemethod_CHP <- multcomp::cld(aux9,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux10<-emmeans::emmeans(foo_slope_CHP,pairwise~TC, adjust="tukey")
#     # lsslopeTC_CHP <- multcomp::cld(aux10,alpha=0.05,Letters=letters,adjust="tukey")
#     
#     aux11<-emmeans::emmeans(foo_z_CHP,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lszmethod_CHP <- multcomp::cld(aux11,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux12<-emmeans::emmeans(foo_z_CHP,pairwise~TC, adjust="tukey")
#     # lszTC_CHP <- multcomp::cld(aux12,alpha=0.05,Letters=letters,adjust="tukey")
#     
#     aux13<-emmeans::emmeans(foo_log_slope_CHP,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lslog_slopemethod_CHP <- multcomp::cld(aux13,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux14<-emmeans::emmeans(foo_log_slope_CHP,pairwise~TC, adjust="tukey")
#     # lslog_slopeTC_CHP <- multcomp::cld(aux14,alpha=0.05,Letters=letters,adjust="tukey")
#     
#     aux15<-emmeans::emmeans(foo_log_ratio_CHP,pairwise~CORRECT_APLICATION, adjust="tukey")
#     lslog_ratiomethod_CHP <- multcomp::cld(aux15,alpha=0.05,Letters=letters,adjust="tukey")
#     # aux16<-emmeans::emmeans(foo_log_ratio_CHP,pairwise~TC, adjust="tukey")
#     # lslog_ratioTC_CHP <- multcomp::cld(aux16,alpha=0.05,Letters=letters,adjust="tukey")
#     
#     # resume table
#     sjt.lmer(foo_log_ratio_TD,foo_slope_TD,foo_log_slope_TD,foo_z_TD,
#              p.kr = FALSE,
#              show.header = TRUE,
#              digits.est = 3,digits.std = 3, digits.p = 3, digits.ci = 3, digits.se = 3,
#              string.est = "Estimate",
#              string.ci = "Conf. Int.",
#              string.p = "p-value",
#              string.dv = "TD",
#              string.pred = "Coefficients",
#              depvar.labels =  c("Ln-Ratio","Slope","Slope (ln-ln)","Z-correlation"))
#     
#     sjt.lmer(foo_log_ratio_CHP,foo_slope_CHP,foo_log_slope_CHP,foo_z_CHP,
#              p.kr = FALSE,
#              show.header = TRUE,
#              digits.est = 3,digits.std = 3, digits.p = 3, digits.ci = 3, digits.se = 3,
#              string.est = "Estimate",
#              string.ci = "Conf. Int.",
#              string.p = "p-value",
#              string.dv = "CHP",
#              string.pred = "Coefficients",
#              depvar.labels =  c("Ln-Ratio","Slope","Slope (ln-ln)","Z-correlation"))

    
#### DENSITY ####
    source('data preparation 1-1.R')
    subdata_dens <- filter(Bar,calibrated=="no",CORRECT_APLICATION == "yes")
    
    foo_slope <- lmer(slope  ~ (method*density + TC) + (1|study) + (1|specie),data=subdata_dens)
    foo_log_slope <- lmer(log_slope  ~ (method*density  + TC) + (1|study) + (1|specie),data=subdata_dens)
    foo_log_ratio <- lmer(log_ratio~ (method*density + TC) + (1|study) + (1|specie),data=subdata_dens)
    foo_z <- lmer(z~ (method*density + TC) + (1|study) + (1|specie),data=subdata_dens)

    anova(foo_slope)
    anova(foo_log_slope)
    anova(foo_log_ratio)
    anova(foo_z)
    
    
    sjt.lmer(foo_log_ratio,foo_slope,foo_log_slope,foo_z,
             p.kr = FALSE,
             show.header = TRUE,
             digits.est = 3,digits.std = 3, digits.p = 3, digits.ci = 3, digits.se = 3,
             string.est = "Estimate",
             string.ci = "Conf. Int.",
             string.p = "p-value",
             string.dv = "Method models",
             string.pred = "Coefficients",
             depvar.labels =  c("Ln-Ratio","Slope","Slope (ln-ln)","Z-correlation"))
    
#lsmean log_ratio
    aux1_dens<-emmeans::emtrends(foo_log_ratio, ~method, var="density")
    lsratiomethod <- multcomp::cld(aux1_dens,alpha=0.05,Letters=letters,adjust="none")
    emmeans::test(emmeans::emtrends(foo_log_ratio,~method,var='density'),null=0)
#lsmean slope
    aux3_dens<- emmeans::emtrends(foo_slope, ~method, var="density")
    lsslopemethod <- multcomp::cld(aux3_dens,alpha=0.05,Letters=letters,adjust="tukey")
    emmeans::test(emmeans::emtrends(foo_slope,~method,var='density'),null=0)
#lsmean log_slope   
    aux2_dens<- emmeans::emtrends(foo_log_slope, ~method, var="density")
    lsslopemethod <- multcomp::cld(aux2_dens,alpha=0.05,Letters=letters,adjust="tukey")
    emmeans::test(emmeans::emtrends(foo_log_slope,~method,var='density'),null=0)
#lsmean z    
    aux4_dens<-(emmeans::emtrends(foo_z, ~method, var="density"))
    lszmethod <- multcomp::cld(aux4_dens,alpha=0.05,Letters=letters,adjust="tukey") 
    emmeans::test(emmeans::emtrends(foo_z,~method,var='density'),null=0)

# density plots    
    eff_df1 <- data.frame(effects::Effect(c('density','method'),foo_log_ratio))
    p1 <- ggplot()+
      geom_line(data=eff_df1,aes(density,fit,color=method))+
      geom_jitter(data=subdata_dens,aes(density,log_ratio,color=method))+
      theme_bw()+
      xlab ( expression(density)) +
      ylab ( expression(Ln-Ratio)) +
      labs(color="",fill="")+
      theme(legend.position="none")
    
    
    eff_df2 <- data.frame(effects::Effect(c('density','method'),foo_slope))
    p2 <- ggplot()+
      geom_line(data=eff_df2,aes(density,fit,color=method))+
      geom_jitter(data=subdata_dens,aes(density,slope,color=method))+
      theme_bw()+
      xlab ( expression(density)) +
      ylab ( expression(Slope)) +
      labs(color="",fill="")+
      theme(legend.position="none")+
      scale_y_continuous(limits = c(-0.1, 2.1))
    
    eff_df3 <- data.frame(effects::Effect(c('density','method'),foo_log_slope))
    p3 <- ggplot()+
      geom_line(data=eff_df3,aes(density,fit,color=method))+
      geom_jitter(data=subdata_dens,aes(density,log_slope,color=method))+
      theme_bw()+
      xlab ( "density") +
      ylab ( "Slope (ln-ln)" ) +
      labs(color="",fill="")+
      theme(legend.position="none")+
      scale_y_continuous(limits = c(-0.1, 2.1))


    eff_df4 <- data.frame(effects::Effect(c('density','method'),foo_z))
    p4 <- ggplot()+
      geom_line(data=eff_df4,aes(x=density,y=fit,color=method))+
      geom_jitter(data=subdata_dens,aes(density,z,color=method))+
      theme_bw()+
      xlab ( expression(density)) +
      ylab ( expression(Z-correlation)) +
      labs(color='',fill="")+
      theme(legend.position="none")

    
    source("multiplot.R")
    multiplot(p1,p2,p3,p4, cols=1)
   

#### POROSITY ####   
    source('data preparation 1-1.R')
    subdataTD <- filter(Bar,calibrated=="no",CORRECT_APLICATION == "yes",method == "TD",porosity!="Monocot")
    subdataCHP <- filter(Bar,calibrated=="no",CORRECT_APLICATION == "yes",method == "CHP",porosity!="Ring porous")
    # models TD
    foo_slope_TD <- lmer(slope  ~ (porosity + TC) + (1|study) + (1|specie),data=subdataTD)
    foo_log_slope_TD <- lmer(log_slope  ~ (porosity + TC) + (1|study) + (1|specie),data=subdataTD)
    foo_log_ratio_TD <- lmer(log_ratio~ (porosity + TC) + (1|study) + (1|specie),data=subdataTD)
    foo_z_TD <- lmer(z~ (porosity + TC) + (1|study) + (1|specie),data=subdataTD)
    
    # models CHP
    foo_slope_CHP <- lmer(slope  ~ (porosity + TC) + (1|study) + (1|specie),data=subdataCHP)
    foo_log_slope_CHP <- lmer(log_slope  ~ (porosity + TC) + (1|study) + (1|specie),data=subdataCHP)
    foo_log_ratio_CHP <- lmer(log_ratio~ (porosity + TC) + (1|study) + (1|specie),data=subdataCHP)
    foo_z_CHP <- lmer(z~ (porosity + TC) + (1|study) + (1|specie),data=subdataCHP)
    
    # lsmeans TD
    aux1<-emmeans::emmeans(foo_log_ratio_TD,pairwise~porosity, adjust="tukey")
    lslog_ratiomethod_TD <- multcomp::cld(aux1$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
    lsmeans::test(emmeans::emmeans(foo_log_ratio_TD,~porosity, adjust="tukey"),null=0)
    
    aux2<-emmeans::emmeans(foo_slope_TD,pairwise~porosity, adjust="tukey")
    lsslopemethod_TD <- multcomp::cld(aux2$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
    lsmeans::test(emmeans::emmeans(foo_slope_TD,~porosity, adjust="tukey"),null=1)
    
    aux3<-emmeans::emmeans(foo_log_slope_TD,pairwise~porosity, adjust="tukey")
    lslog_slopemethod_TD <- multcomp::cld(aux3$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
    lsmeans::test(emmeans::emmeans(foo_log_slope_TD,~porosity, adjust="tukey"),null=1)

    aux4<-emmeans::emmeans(foo_z_TD,pairwise~porosity, adjust="tukey")
    lszmethod_TD <- multcomp::cld(aux4,alpha=0.05,Letters=letters,adjust="tukey")
    
    #n species
    group_by(subdataTD,porosity,method)%>%summarize(n())


    # lsmeans CHP
    aux5<-emmeans::emmeans(foo_log_ratio_CHP,pairwise~porosity, adjust="tukey")
    lslog_ratiomethod_CHP <- multcomp::cld(aux5$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
    lsmeans::test(emmeans::emmeans(foo_log_ratio_CHP,~porosity, adjust="tukey"),null=0)
    
    aux6<-emmeans::emmeans(foo_slope_CHP,pairwise~porosity, adjust="tukey")
    lsslopemethod_CHP <- multcomp::cld(aux6$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
    lsmeans::test(emmeans::emmeans(foo_slope_CHP,~porosity, adjust="tukey"),null=1)
    
    aux7<-emmeans::emmeans(foo_log_slope_CHP,pairwise~porosity, adjust="tukey")
    lslog_slopemethod_CHP <- multcomp::cld(aux7$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
    lsmeans::test(emmeans::emmeans(foo_log_slope_CHP,~porosity, adjust="tukey"),null=1)
    
    aux8<-emmeans::emmeans(foo_z_CHP,pairwise~porosity, adjust="tukey")
    lszmethod_CHP <- multcomp::cld(aux8$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
    
    group_by(subdataCHP,porosity,method)%>%summarize(n())
    
    # resume table
    sjt.lmer(foo_log_ratio_TD,foo_slope_TD,foo_log_slope_TD,foo_z_TD,
             p.kr = FALSE,
             show.header = TRUE,
             digits.est = 3,digits.std = 3, digits.p = 3, digits.ci = 3, digits.se = 3,
             string.est = "Estimate",
             string.ci = "Conf. Int.",
             string.p = "p-value",
             string.dv = "TD",
             string.pred = "Coefficients",
             depvar.labels =  c("Ln-Ratio","Slope","Slope (ln-ln)","Z-correlation"))
    
    sjt.lmer(foo_log_ratio_CHP,foo_slope_CHP,foo_log_slope_CHP,foo_z_CHP,
             p.kr = FALSE,
             show.header = TRUE,
             digits.est = 3,digits.std = 3, digits.p = 3, digits.ci = 3, digits.se = 3,
             string.est = "Estimate",
             string.ci = "Conf. Int.",
             string.p = "p-value",
             string.dv = "CHP",
             string.pred = "Coefficients",
             depvar.labels =  c("Ln-Ratio","Slope","Slope (ln-ln)","Z-correlation"))
    
    


#### RMSE ####
source('data preparation 1-1-density.R')
source('data preparation 1-1-totalflow.R')
source('data preparation 1-1.R')
subdata <- filter(Bar,calibrated == "no",CORRECT_APLICATION == "yes")
subdata$method <- factor(subdata$method, levels (subdata$method) [c(1,5,3,2,4,6,7)])

# foo_rmse <- lmer(RMSE ~  range_mean * method    + (1|study) + (1|specie),data=subdata)
foo_rmse <- lmer(NRMSE ~  range_mean  * method    + (1|study) + (1|specie),data=subdata)
# foo_rmse <- lmer(NRMSE ~  (range_mean + I(range_mean^2)) * method    + (1|study) + (1|specie),data=subdata)

# #SFD
# foo_rmse <- lmer(rmse ~  range_mean * method    + (1|study) + (1|specie),data=subdata)
# #SF
# foo_rmse <- lmer(rmse ~  range_mean * method   + (1|study) + (1|specie),data=subdata)

summary(foo_rmse)
anova(foo_rmse)

plot(effects::effect('median*method',foo_rmse,default.levels=20),multiline=TRUE)

plot(effects::Effect(focal.predictors = c("range_mean","method"), mod = foo_rmse,
            xlevels=list(method=seq(0, 100, 1), median=0:100)))
plot(effects::Effect(focal.predictors = c("range_mean","method"), mod = foo_rmse,
                     xlevels=list(method=seq(0, 100, 1), range_mean=0:100)))

#SFD
lsmeans::test(emmeans::emmeans(foo_rmse,~range_mean * method, adjust="tukey",at = list(range_mean = 0)),null=0)
lsmeans::test(lsmeans::lstrends(foo_rmse,~method,var='range_mean'),null=1)
# lsmeans::test(lsmeans::lstrends(foo_rmse,~method,var='I(range_mean^2)'),null=1)
lsmeans::test(emmeans::emmeans(foo_rmse,~range_mean * method, adjust="tukey",at = list(range_mean = 25)),null=0)
lsmeans::test(emmeans::emmeans(foo_rmse,~range_mean*method , adjust="tukey",at = list(range_mean = 1300)),null=0)

# #SF
# lsmeans::test(emmeans::emmeans(foo_rmse,~range_mean*method , adjust="tukey",at = list(range_mean = 1300)),null=0)
# lsmeans::test(lsmeans::lstrends(foo_rmse,~method,var='range_mean'),null=1)
# 
# MuMIn::std.coef(foo_rmse,partial.sd = FALSE)



# subdata%>%group_by(study,calibrations,porosity,sensor.length,method)%>%summarise()%>%filter(porosity=='Ring porous' & method=='TD') -> g



#### conduit size ####
source('data preparation 1-1.R')
subdata_dens <- filter(Bar,calibrated=="no",CORRECT_APLICATION == "yes")
subdata_dens <- filter(mutate(subdata_dens,prop_area_vessel =vessel_area*vessel_density),prop_area_vessel < 1)

foo_slope <- lmer(slope  ~ (method*prop_area_vessel) + (1|study) + (1|specie),data=subdata_dens)
foo_log_slope <- lmer(log_slope  ~ (method*prop_area_vessel) + (1|study) + (1|specie),data=subdata_dens)
foo_log_ratio <- lmer(log_ratio~ (method*prop_area_vessel) + (1|study) + (1|specie),data=subdata_dens)
foo_z <- lmer(z~ (method*prop_area_vessel) + (1|study) + (1|specie),data=subdata_dens)

anova(foo_slope)
anova(foo_log_slope)
anova(foo_log_ratio)
anova(foo_z)


sjt.lmer(foo_log_ratio,foo_slope,foo_log_slope,foo_z,
         p.kr = FALSE,
         show.header = TRUE,
         digits.est = 3,digits.std = 3, digits.p = 3, digits.ci = 3, digits.se = 3,
         string.est = "Estimate",
         string.ci = "Conf. Int.",
         string.p = "p-value",
         string.dv = "Method models",
         string.pred = "Coefficients",
         depvar.labels =  c("Ln-Ratio","Slope","Slope (ln-ln)","Z-correlation"))

#lsmean log_ratio
aux1_dens<-lsmeans::lstrends(foo_log_ratio, ~method, var="prop_area_vessel")
lsratiomethod <- multcomp::cld(aux1_dens$emmeans,alpha=0.05,Letters=letters,adjust="none")
lsmeans::test(lsmeans::lstrends(foo_log_ratio,~method,var='prop_area_vessel'),null=0)
#lsmean slope
aux3_dens<- lsmeans::lstrends(foo_slope, ~method, var="prop_area_vessel")
lsslopemethod <- multcomp::cld(aux3_dens$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
lsmeans::test(lsmeans::lstrends(foo_slope,~method,var='prop_area_vessel'),null=0)
#lsmean log_slope   
aux2_dens<- lsmeans::lstrends(foo_log_slope, ~method, var="prop_area_vessel")
lsslopemethod <- multcomp::cld(aux2_dens$emmeans,alpha=0.05,Letters=letters,adjust="tukey")
lsmeans::test(lsmeans::lstrends(foo_log_slope,~method,var='prop_area_vessel'),null=0)
#lsmean z    
aux4_dens<-(lsmeans::lstrends(foo_z, ~method, var="prop_area_vessel"))
lszmethod <- multcomp::cld(aux4_dens$emmeans,alpha=0.05,Letters=letters,adjust="tukey") 
lsmeans::test(lsmeans::lstrends(foo_z,~method,var='prop_area_vessel'),null=0)

# density plots    
eff_df1 <- data.frame(effects::Effect(c('prop_area_vessel','method'),foo_log_ratio))
p1 <- ggplot()+
  geom_line(data=eff_df1,aes(prop_area_vessel,fit,color=method))+
  geom_jitter(data=subdata_dens,aes(prop_area_vessel,log_ratio,color=method))+
  theme_bw()+
  xlab ( expression(density)) +
  ylab ( expression(Ln-Ratio)) +
  labs(color="",fill="")+
  theme(legend.position="none")


eff_df2 <- data.frame(effects::Effect(c('prop_area_vessel','method'),foo_slope))
p2 <- ggplot()+
  geom_line(data=eff_df2,aes(prop_area_vessel,fit,color=method))+
  geom_jitter(data=subdata_dens,aes(prop_area_vessel,slope,color=method))+
  theme_bw()+
  xlab ( expression(density)) +
  ylab ( expression(Slope)) +
  labs(color="",fill="")+
  theme(legend.position="none")+
  scale_y_continuous(limits = c(-0.1, 2.1))

eff_df3 <- data.frame(effects::Effect(c('prop_area_vessel','method'),foo_log_slope))
p3 <- ggplot()+
  geom_line(data=eff_df3,aes(prop_area_vessel,fit,color=method))+
  geom_jitter(data=subdata_dens,aes(prop_area_vessel,log_slope,color=method))+
  theme_bw()+
  xlab ( "density") +
  ylab ( "Slope (ln-ln)" ) +
  labs(color="",fill="")+
  theme(legend.position="none")+
  scale_y_continuous(limits = c(-0.1, 2.1))


eff_df4 <- data.frame(effects::Effect(c('prop_area_vessel','method'),foo_z))
p4 <- ggplot()+
  geom_line(data=eff_df4,aes(x=prop_area_vessel,y=fit,color=method))+
  geom_jitter(data=subdata_dens,aes(prop_area_vessel,z,color=method))+
  theme_bw()+
  xlab ( expression(density)) +
  ylab ( expression(Z-correlation)) +
  labs(color='',fill="")+
  theme(legend.position="none")


source("multiplot.R")
multiplot(p1,p2,p3,p4, cols=1)


#### SFD vs SF ####
source('data preparation 1-1.R')
# source('data preparation 1-1-density.R')
subdata1 <- filter(Bar,calibrated=="no",CORRECT_APLICATION == "yes",method == "TD")
subdata1$measuring.type <- as.factor(subdata1$measuring.type)
foo_slope <- lmer(slope  ~ (measuring.type + TC) + (1|study) + (1|specie),data=subdata1)
foo_log_slope <- lmer(log_slope  ~ (measuring.type + TC) + (1|study) + (1|specie),data=subdata1)
foo_log_ratio <- lmer(log_ratio~ (measuring.type+TC) + (1|study) + (1|specie),data=subdata1)
foo_z <- lmer(z~ (measuring.type+TC) + (1|study) + (1|specie),data=subdata1)
summary(foo_log_ratio)
summary(foo_slope)
summary(foo_log_slope)
summary(foo_z)
anova(foo_log_ratio)
anova(foo_slope)
anova(foo_log_slope)
anova(foo_z)


subdata1 <- filter(Bar,calibrated=="no",method == "CHP")
