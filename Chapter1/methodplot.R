library(ggplot2)
MFplot <- function(data=NULL,intercept=1,xlab=NULL,ylab=NULL,ylimits=c(0,1)){ggplot(data)+
    geom_point(aes(x=MF,y=emmean,color=MF))+
    geom_errorbar(aes(x=MF, ymin=lower.CL, ymax=upper.CL,colour=MF), width=.1)+
    geom_abline(slope = 0,intercept = intercept,linetype=3)+
    theme_bw()+
    xlab (xlab) +
    ylab (ylab) +
    labs(color="",fill="")+
    theme(legend.position="none")+
    geom_text(aes(x=MF,y=upper.CL,label=.group),
              vjust = -0.5,
              color   = "black")+
    scale_y_continuous(limits = ylimits)}

methodplot <- function(data=NULL,intercept=1,xlab=NULL,ylab=NULL,ylimits=c(0,1)){ggplot(data)+
    geom_point(aes(x=method,y=emmean,color=method))+
    geom_errorbar(aes(x=method, ymin=lower.CL, ymax=upper.CL,colour=method), width=.1)+
    geom_abline(slope = 0,intercept = intercept,linetype=3)+
    theme_bw()+
    xlab (xlab) +
    ylab (ylab) +
    labs(color="",fill="")+
    theme(legend.position="none")+
    geom_text(aes(x=method,y=upper.CL,label=.group),
              vjust = -0.5,
              color   = "black")+
    scale_y_continuous(limits = ylimits)}

TCplot <- function(data=NULL,intercept=1,xlab=NULL,ylab=NULL,ylimits=c(0,1)){ggplot(data)+
    geom_point(aes(x=TC,y=emmean,color=TC))+
    geom_errorbar(aes(x=TC, ymin=lower.CL, ymax=upper.CL,colour=TC), width=.1)+
    geom_abline(slope = 0,intercept = intercept,linetype=3)+
    theme_bw()+
    xlab (xlab) +
    ylab (ylab) +
    labs(color="",fill="")+
    theme(legend.position="none")+
    geom_text(aes(x=TC,y=upper.CL,label=.group),
              vjust = -0.5,
              color   = "black")+
    scale_y_continuous(limits = ylimits)}