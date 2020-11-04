calc.relip.mm <- function(model,type='lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[,-1]
  Y <- getME(model,'y')
  s_resid <- sigma(model)
  s_effect <- getME(model,'theta')*s_resid
  s2 <- sum(s_resid^2,s_effect^2)
  V <- Diagonal(x = s2,n=nrow(X))
  YX <- cbind(Y,X)
  cov_XY <- solve( t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  importances <- calc.relimp(as.matrix(cov_XY),rela=T,type=type)
  return(importances)
}


r2MLM <- function(data,within_covs,between_covs,random_covs,
                  gamma_w,gamma_b,Tau,sigma2,has_intercept=T,clustermeancentered=T){
  if(has_intercept==T){
    if(length(gamma_b)>1) gamma <- c(1,gamma_w,gamma_b[2:length(gamma_b)])
    if(length(gamma_b)==1) gamma <- c(1,gamma_w)
    if(is.null(within_covs)==T) gamma_w <- 0
  }
  if(has_intercept==F){
    gamma <- c(gamma_w,gamma_b)
    if(is.null(within_covs)==T) gamma_w <- 0
    if(is.null(between_covs)==T) gamma_b <- 0
  }
  if(is.null(gamma)) gamma <- 0
  ##compute phi
  phi <- var(cbind(1,data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(has_intercept==F) phi <- var(cbind(data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(is.null(within_covs)==T & is.null(within_covs)==T & has_intercept==F) phi <- 0
  phi_w <- var(data[,within_covs],na.rm=T)
  if(is.null(within_covs)==T) phi_w <- 0
  phi_b <- var(cbind(1,data[,between_covs]),na.rm=T)
  if(is.null(between_covs)==T) phi_b <- 0
  ##compute psi and kappa
  var_randomcovs <- var(cbind(1,data[,c(random_covs)]),na.rm=T)
  if(length(Tau)>1) psi <- matrix(c(diag(Tau)),ncol=1)
  if(length(Tau)==1) psi <- Tau
  if(length(Tau)>1) kappa <- matrix(c(Tau[lower.tri(Tau)==TRUE]),ncol=1)
  if(length(Tau)==1) kappa <- 0
  v <- matrix(c(diag(var_randomcovs)),ncol=1)
  r <- matrix(c(var_randomcovs[lower.tri(var_randomcovs)==TRUE]),ncol=1)
  if(is.null(random_covs)==TRUE){
    v <- 0
    r <- 0
    m <- matrix(1,ncol=1)
  }
  if(length(random_covs)>0) m <- matrix(c(colMeans(cbind(1,data[,c(random_covs)]),na.rm=T)),ncol=1)
  ##total variance
  totalvar_notdecomp <- t(v)%*%psi + 2*(t(r)%*%kappa) + t(gamma)%*%phi%*%gamma + t(m)%*%Tau%*%m + sigma2
  totalwithinvar <- (t(gamma_w)%*%phi_w%*%gamma_w) + (t(v)%*%psi + 2*(t(r)%*%kappa)) + sigma2
  totalbetweenvar <- (t(gamma_b)%*%phi_b%*%gamma_b) + Tau[1]
  totalvar <- totalwithinvar + totalbetweenvar
  ##total decomp
  decomp_fixed_notdecomp <- (t(gamma)%*%phi%*%gamma) / totalvar_notdecomp
  decomp_varslopes_notdecomp <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalvar_notdecomp
  decomp_varmeans_notdecomp <- (t(m)%*%Tau%*%m) / totalvar_notdecomp
  decomp_sigma_notdecomp <- sigma2/totalvar_notdecomp
  decomp_fixed_within <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalvar
  decomp_fixed_between <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalvar
  decomp_fixed <- decomp_fixed_within + decomp_fixed_between
  decomp_varslopes <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalvar
  decomp_varmeans <- (t(m)%*%Tau%*%m) / totalvar
  decomp_sigma <- sigma2/totalvar
  ##within decomp
  decomp_fixed_within_w <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalwithinvar
  decomp_varslopes_w <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalwithinvar
  decomp_sigma_w <- sigma2/totalwithinvar
  ##between decomp
  decomp_fixed_between_b <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalbetweenvar
  decomp_varmeans_b <- Tau[1] / totalbetweenvar
  #NEW measures
  if (clustermeancentered==TRUE){
    R2_f <- decomp_fixed
    R2_f1 <- decomp_fixed_within
    R2_f2 <- decomp_fixed_between
    R2_fv <- decomp_fixed + decomp_varslopes
    R2_fvm <- decomp_fixed + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
    R2_f_w <- decomp_fixed_within_w
    R2_f_b <- decomp_fixed_between_b
    R2_fv_w <- decomp_fixed_within_w + decomp_varslopes_w
    R2_v_w <- decomp_varslopes_w
    R2_m_b <- decomp_varmeans_b
  }
  if (clustermeancentered==FALSE){
    R2_f <- decomp_fixed_notdecomp
    R2_fv <- decomp_fixed_notdecomp + decomp_varslopes_notdecomp
    R2_fvm <- decomp_fixed_notdecomp + decomp_varslopes_notdecomp + decomp_varmeans_notdecomp
    R2_v <- decomp_varslopes_notdecomp
    R2_m <- decomp_varmeans_notdecomp
  }
  if(clustermeancentered==TRUE){
    decomp_table <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                             decomp_fixed_within_w,"NA",decomp_varslopes_w,"NA",decomp_sigma_w,
                             "NA",decomp_fixed_between_b,"NA",decomp_varmeans_b,"NA"),ncol=3)
    rownames(decomp_table) <- c("fixed, within","fixed, between","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total","within","between")
    R2_table <- matrix(c(R2_f1,R2_f2,R2_v,R2_m,R2_f,R2_fv,R2_fvm,
                         R2_f_w,"NA",R2_v_w,"NA","NA",R2_fv_w,"NA",
                         "NA",R2_f_b,"NA",R2_m_b,"NA","NA","NA")
                       ,ncol=3)
    rownames(R2_table) <- c("f1","f2","v","m","f","fv","fvm")
    colnames(R2_table) <- c("total","within","between")
  }
  ##barchart
  if(clustermeancentered==TRUE){
    contributions_stacked <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                                      decomp_fixed_within_w,0,decomp_varslopes_w,0,decomp_sigma_w,
                                      0,decomp_fixed_between_b,0,decomp_varmeans_b,0),5,3)
    colnames(contributions_stacked) <- c("total","within","between")
    rownames(contributions_stacked) <- c("fixed slopes (within)",
                                         "fixed slopes (between)",
                                         "slope variation (within)",
                                         "intercept variation (between)",
                                         "residual (within)")
    barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
            ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
            density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),xlim=c(0,1),width=c(.3,.3))
    legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkred","steelblue","darkred","midnightblue","white"),
           cex=.7, pt.cex = 1,xpd=T,density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0))
  }
  if(clustermeancentered==FALSE){
    decomp_table <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes_notdecomp,decomp_varmeans_notdecomp,decomp_sigma_notdecomp),ncol=1)
    rownames(decomp_table) <- c("fixed","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total")
    R2_table <- matrix(c(R2_f,R2_v,R2_m,R2_fv,R2_fvm),ncol=1)
    rownames(R2_table) <- c("f","v","m","fv","fvm")
    colnames(R2_table) <- c("total")
    ##barchar
    contributions_stacked <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes_notdecomp,decomp_varmeans_notdecomp,decomp_sigma_notdecomp),4,1)
    colnames(contributions_stacked) <- c("total")
    rownames(contributions_stacked) <- c("fixed slopes",
                                         "slope variation",
                                         "intercept variation",
                                         "residual")
    barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
            ylim=c(0,1),col=c("darkblue","darkblue","darkblue","white"),ylab="proportion of variance",
            density=c(NA,30,40,NA),angle=c(0,0,135,0),xlim=c(0,1),width=c(.6))
    legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkblue","darkblue","darkblue","white"),
           cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,30,40,NA),angle=c(0,0,135,0))
  }
  Output <- list(noquote(decomp_table),noquote(R2_table))
  names(Output) <- c("Decompositions","R2s")
  return(Output)
}




r2MLMcomp <- function(data,within_covs_modA,between_covs_modA,random_covs_modA,
                      gamma_w_modA,gamma_b_modA,Tau_modA,sigma2_modA,
                      within_covs_modB,between_covs_modB,random_covs_modB,
                      gamma_w_modB,gamma_b_modB,Tau_modB,sigma2_modB){
  ##r2MLM function
  r2MLM <- function(data,within_covs,between_covs,random_covs,
                    gamma_w,gamma_b,Tau,sigma2,modelname){
    if(length(gamma_b)>1) gamma <- c(1,gamma_w,gamma_b[2:length(gamma_b)])
    if(length(gamma_b)==1) gamma <- c(1,gamma_w)
    if(is.null(within_covs)==T) gamma_w <- 0
    if(is.null(gamma)) gamma <- 0 
    ##compute phi
    phi <- var(cbind(1,data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
    phi_w <- var(data[,within_covs],na.rm=T)
    if(is.null(within_covs)==T) phi_w <- 0
    phi_b <- var(cbind(1,data[,between_covs]),na.rm=T)
    if(is.null(between_covs)==T) phi_b <- 0    
    ##compute psi and kappa
    var_randomcovs <- var(cbind(1,data[,c(random_covs)]),na.rm=T)
    if(length(Tau)>1) psi <- matrix(c(diag(Tau)),ncol=1)
    if(length(Tau)==1) psi <- Tau
    if(length(Tau)>1) kappa <- matrix(c(Tau[lower.tri(Tau)==TRUE]),ncol=1)
    if(length(Tau)==1) kappa <- 0
    v <- matrix(c(diag(var_randomcovs)),ncol=1)
    r <- matrix(c(var_randomcovs[lower.tri(var_randomcovs)==TRUE]),ncol=1)
    if(is.null(random_covs)==TRUE){
      v <- 0
      r <- 0
      m <- matrix(1,ncol=1)
    } 
    if(length(random_covs)>0) m <- matrix(c(colMeans(cbind(1,data[,c(random_covs)]),na.rm=T)),ncol=1)
    ##total variance
    totalvar_notdecomp <- t(v)%*%psi + 2*(t(r)%*%kappa) + t(gamma)%*%phi%*%gamma + t(m)%*%Tau%*%m + sigma2
    totalwithinvar <- (t(gamma_w)%*%phi_w%*%gamma_w) + (t(v)%*%psi + 2*(t(r)%*%kappa)) + sigma2
    totalbetweenvar <- (t(gamma_b)%*%phi_b%*%gamma_b) + Tau[1]
    totalvar <- totalwithinvar + totalbetweenvar
    ##total decomp 
    decomp_fixed_notdecomp <- (t(gamma)%*%phi%*%gamma) / totalvar
    decomp_fixed_within <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalvar
    decomp_fixed_between <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalvar
    decomp_fixed <- decomp_fixed_within + decomp_fixed_between
    decomp_varslopes <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalvar
    decomp_varmeans <- (t(m)%*%Tau%*%m) / totalvar
    decomp_sigma <- sigma2/totalvar
    ##within decomp
    decomp_fixed_within_w <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalwithinvar
    decomp_varslopes_w <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalwithinvar
    decomp_sigma_w <- sigma2/totalwithinvar
    ##between decomp
    decomp_fixed_between_b <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalbetweenvar
    decomp_varmeans_b <- Tau[1] / totalbetweenvar
    ##measures
    R2_f <- decomp_fixed
    R2_f1 <- decomp_fixed_within 
    R2_f2 <- decomp_fixed_between 
    R2_fv <- decomp_fixed + decomp_varslopes
    R2_fvm <- decomp_fixed + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
    R2_f_w <- decomp_fixed_within_w
    R2_f_b <- decomp_fixed_between_b
    R2_fv_w <- decomp_fixed_within_w + decomp_varslopes_w
    R2_v_w <- decomp_varslopes_w
    R2_m_b <- decomp_varmeans_b
    decomp_table <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                             decomp_fixed_within_w,"NA",decomp_varslopes_w,"NA",decomp_sigma_w,
                             "NA",decomp_fixed_between_b,"NA",decomp_varmeans_b,"NA"),ncol=3)
    rownames(decomp_table) <- c("fixed, within","fixed, between","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total","within","between")
    R2_table <- matrix(c(R2_f1,R2_f2,R2_v,R2_m,R2_f,R2_fv,R2_fvm,
                         R2_f_w,"NA",R2_v_w,"NA","NA",R2_fv_w,"NA",
                         "NA",R2_f_b,"NA",R2_m_b,"NA","NA","NA")
                       ,ncol=3)
    rownames(R2_table) <- c("f1","f2","v","m","f","fv","fvm")
    colnames(R2_table) <- c("total","within","between")
    ##barchart
    contributions_stacked <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                                      decomp_fixed_within_w,0,decomp_varslopes_w,0,decomp_sigma_w,
                                      0,decomp_fixed_between_b,0,decomp_varmeans_b,0),5,3)
    colnames(contributions_stacked) <- c("total","within","between")
    rownames(contributions_stacked) <- c("fixed slopes (within)",
                                         "fixed slopes (between)",
                                         "slope variation (within)",
                                         "intercept variation (between)",
                                         "residual (within)")
    barplot(contributions_stacked, main=paste0("Decomposition of Scaled Variance, Model ",modelname), horiz=FALSE,
            ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
            density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),xlim=c(0,1),width=c(.3,.3))
    legend(.33,-.1,legend=rownames(contributions_stacked),fill=c("darkred","steelblue","darkred","midnightblue","white"),
           cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0))
    Output <- list(noquote(decomp_table),noquote(R2_table))
    names(Output) <- c("Decompositions","R2s")
    return(Output)
  }
  ##compute decomp for Model A and B
  results_modA <- r2MLM(data,within_covs_modA,between_covs_modA,random_covs_modA,
                        gamma_w_modA,gamma_b_modA,Tau_modA,sigma2_modA,"A")
  decomp_modA <- results_modA$Decompositions
  results_modB <- r2MLM(data,within_covs_modB,between_covs_modB,random_covs_modB,
                        gamma_w_modB,gamma_b_modB,Tau_modB,sigma2_modB,"B")
  decomp_modB <- results_modB$Decompositions
  ##comparison measures
  delta_f1_t <- as.numeric(decomp_modA[1,1]) - as.numeric(decomp_modB[1,1]) 
  delta_f2_t <- as.numeric(decomp_modA[2,1]) - as.numeric(decomp_modB[2,1]) 
  delta_v_t <- as.numeric(decomp_modA[3,1]) - as.numeric(decomp_modB[3,1]) 
  delta_m_t <-as.numeric(decomp_modA[4,1]) - as.numeric(decomp_modB[4,1]) 
  delta_f1_w <- as.numeric(decomp_modA[1,2]) - as.numeric(decomp_modB[1,2]) 
  delta_v_w <- as.numeric(decomp_modA[3,2]) - as.numeric(decomp_modB[3,2]) 
  delta_f2_b <- as.numeric(decomp_modA[2,3]) - as.numeric(decomp_modB[2,3]) 
  delta_m_b <- as.numeric(decomp_modA[4,3]) - as.numeric(decomp_modB[4,3]) 
  delta_f_t <- delta_f1_t + delta_f2_t
  delta_fv_t <- delta_f1_t + delta_f2_t + delta_v_t
  delta_fvm_t <- delta_f1_t + delta_f2_t + delta_v_t + delta_m_t 
  delta_f1v_w <- delta_f1_w + delta_v_w
  ##comparison bar charts
  contributions_stacked_total <-  matrix(c(as.numeric(decomp_modA[1,1]),as.numeric(decomp_modA[2,1]),as.numeric(decomp_modA[3,1]),as.numeric(decomp_modA[4,1]),as.numeric(decomp_modA[5,1]),              as.numeric(decomp_modB[1,1]),as.numeric(decomp_modB[2,1]),as.numeric(decomp_modB[3,1]),as.numeric(decomp_modB[4,1]),as.numeric(decomp_modB[5,1])),5,2)
  colnames(contributions_stacked_total) <- c("Model A","Model B")
  barplot(contributions_stacked_total, main="Decomposition of Scaled Total Variance", horiz=FALSE,
          ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
          density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),width=c(.3,.3))
  legend(0.26,-.1,legend=c("fixed slopes (within)",
                           "fixed slopes (between)",
                           "slope variation (within)",
                           "intercept variation (between)",
                           "residual (within)"),fill=c("darkred","steelblue","darkred","midnightblue","white"),
         cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0))
  contributions_stacked_within <-  matrix(c(as.numeric(decomp_modA[1,2]),0,as.numeric(decomp_modA[3,2]),0,as.numeric(decomp_modA[5,2]),
                                            as.numeric(decomp_modB[1,2]),0,as.numeric(decomp_modB[3,2]),0,as.numeric(decomp_modB[5,2])),5,2)
  colnames(contributions_stacked_within) <- c("Model A","Model B")
  barplot(contributions_stacked_within, main="Decomposition of Scaled Within-Cluster Variance", horiz=FALSE,
          ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
          density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),width=c(.3,.3))
  legend(0.28,-.1,legend=c("fixed slopes (within)",
                           "slope variation (within)",
                           "residual (within)"),fill=c("darkred","darkred","white"),
         cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,30,NA),angle=c(0,0,0))
  contributions_stacked_between <-  matrix(c(0,as.numeric(decomp_modA[2,3]),0,as.numeric(decomp_modA[4,3]),0,
                                             0,as.numeric(decomp_modB[2,3]),0,as.numeric(decomp_modB[4,3]),0),5,2)
  colnames(contributions_stacked_between) <- c("Model A","Model B")
  barplot(contributions_stacked_between, main="Decomposition of Scaled Between-Cluster Variance", horiz=FALSE,
          ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
          density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),width=c(.3,.3))
  legend(0.26,-.1,legend=c("fixed slopes (between)",
                           "intercept variation (between)"),fill=c("steelblue","midnightblue"),
         cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,40),angle=c(45,135))
  ##table of R2 deltas
  R2_modA <- results_modA$R2s
  R2_modB <- results_modB$R2s
  R2_delta <- suppressWarnings(as.numeric(R2_modB) - as.numeric(R2_modA))
  R2_delta <- matrix(R2_delta,7,3)
  colnames(R2_delta) <- colnames(R2_modA)
  rownames(R2_delta) <- rownames(R2_modA)
  Output <- list(R2_modA,R2_modB,R2_delta)
  names(Output) <- c("Model A R2s","Model B R2s","R2 differences, Model B - Model A")
  return(Output)
}

