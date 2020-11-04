library(tidyverse)
library(sjPlot)
library(gghighlight)
library(cowplot)

# Load model object from Victor
load('model_filt.RData')

# Getting data from model
datamod <- model@frame

# renaming vars
names(datamod) <- c('Gsw','log_vpd','log_swc','pl_species','pl_code')

# models
# mod.victor <- lmerTest::lmer(Gsw ~ (log_vpd + log_swc | pl_species) + (1 | pl_species:pl_code),data=datamod)
sjPlot::plot_model(model,type='re')

summary(model)
coef(model)

View(coef(model)$pl_species)
View(coef(mod1)$`pl_species:pl_code`)


coef.victor<- rownames_to_column(coef(model)$pl_species,'pl_species')
names(coef.victor) <- c('pl_species','log_vpd_vic','log_swc_vic','int_vic')


# remove intercepts in vpd/swc responses, add a (independent?) species intercept

mod1 <- lmerTest::lmer(Gsw ~ (0+log_vpd|pl_species) + (0+log_swc|pl_species) + 
                         (1|pl_species) +(1| pl_species:pl_code),data=datamod,
                       REML=TRUE,
                       control = lmerControl(optimizer = "optimx",optCtrl = list(method = "nlminb")))



# rough check correlations vpd-swc per species
datamod %>% 
  group_by(pl_species) %>% 
  summarise(rcor = cor(log_vpd,log_swc)) -> datacor

coef.victor %>% 
  left_join(datacor) %>% 
  ggplot(aes(y=log_swc_vic,x=rcor,col=pl_species))+geom_point()+
  gghighlight(abs(rcor)>0.3 & log_swc_vic<0)


coef.victor %>% 
  summarise(neg = sum(log_swc_vic<0))


coef.mod1<- rownames_to_column(coef(mod1)$pl_species,'pl_species')
names(coef.mod1) <- c('pl_species','log_vpd_mod1','log_swc_mod1','int_mod1')


View(coef.mods)


coef.victor %>% left_join(coef.mod1)->coef.mods


plot_int<- ggplot(coef.mods,aes(x=int_vic,y=int_mod1,col=pl_species))+
  geom_point()+geom_abline()+
  gghighlight(int_vic<0)

plot_vpd <- ggplot(coef.mods,aes(x=-log_vpd_vic,y=-log_vpd_mod1,col=pl_species))+
  geom_point()+geom_abline()+
  gghighlight(log_vpd_vic<0)

plot_swc<- ggplot(coef.mods,aes(x=log_swc_vic,y=log_swc_mod1,col=pl_species))+
  geom_point()+geom_abline()+ 
  gghighlight(log_swc_vic<(-1*60))

cowplot::plot_grid(plot_int,plot_vpd,plot_swc,ncol=3)


plot_int_swc_vic<- ggplot(coef.mods,aes(x=int_vic,y=log_swc_vic,col=pl_species))+
  geom_point()+geom_abline()+
  gghighlight(int_vic<0)

plot_int_vpd_vic<- ggplot(coef.mods,aes(x=int_vic,y=-log_vpd_vic,col=pl_species))+
  geom_point()+geom_abline()+
  gghighlight(int_vic<0)


plot_int_swc_mod1<- ggplot(coef.mods,aes(x=int_mod1,y=log_swc_mod1,col=pl_species))+
  geom_point()+geom_abline()+
  gghighlight(int_vic<0)

plot_int_vpd_mod1<- ggplot(coef.mods,aes(x=int_mod1,y=-log_vpd_mod1,col=pl_species))+
  geom_point()+geom_abline()+
  gghighlight(int_vic<0)

cowplot::plot_grid(plot_int_swc_vic,plot_int_vpd_vic,plot_int_swc,plot_int_vpd,ncol=2,nrow=2)



