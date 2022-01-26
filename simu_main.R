#### library

load.lib<-c("survival", "survminer","survRM2","dplyr", "tidyr","ranger", "glmnet", "twang", "mvtnorm")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE,repos = "http://cran.us.r-project.org")
sapply(load.lib, require, character=TRUE)
# source functions

ddir <- "/PSWeight.SGA/design"
adir <- "/PSWeight.SGA/analysis"

#design
source(paste0(ddir,"print_Sumstat_sub.R"))
source(paste0(ddir,"Plot_SumStat_sub.R"))
source(paste0(ddir,'PSmethod_sub.R'))
source(paste0(ddir,"SumStat_sub.R"))

#analysis
source(paste0(adir,"print_PSWeight_sub.R"))
source(paste0(adir,"summary_PSWeight_sub.R"))
source(paste0(adir,"print_PSweightsum_sub.R"))
source(paste0(adir,"PSweight_sub.R"))

# helper functions
source("datagen_520.R")
source("helper_overall.R")
source("RMST_sub_overall.R")



all_params <- expand.grid(ns = c(1000, 3000),
                          gammas= c(0.5,1),
                          kappas = c(-1/2, -1/4),
                          taus = c(0, -log(2)),
                          PHs=c(TRUE, FALSE)
)

nsim <-1000
nmod <- 5
load(paste0(outdir,"_true_val_censor_100_2021-09-13_.Rdata"))
all_results <- NULL
  #####################################################################
for(j in 1:nrow(all_params)){
  n =  all_params[j, "ns"]; kappa=all_params[j, "kappas"];
  gamma =  all_params[j, "gammas"];  tau =  all_params[j, "taus"];
  PH =all_params[j, "PHs"]
  print(c(n,kappa,gamma, tau, PH))
  

  set.seed(123)
  ## 1) Generate data
  data1 = sim_data(n=1000000, tau,gamma, kappa, PH=PH, censoring= TRUE)
  data1 = as.data.frame(data1)
  data1[,"ID"] =seq(1,nrow(data1))
   
  
  # subgroup variable
   subgroup= paste0("G",1:2)
  
  
  main.form <- as.formula(paste("Z~",paste(paste("X",1:10,sep="", collapse="+"),paste("G",1:2,sep="", collapse="+"),sep="+")))
    # truncated time
  truncate <- 365

  EST_MHR_OW <- EST_MHR_IPW <- SE_MHR_OW <- SE_MHR_IPW <- COVER_MHR_OW <- COVER_MHR_IPW<- 
  EST_RATE <- EST_RATO <- SE_RATE <- SE_RATO <- COVER_RATE <- COVER_RATO <- array(NA,dim=c(5,5, nsim))
  
  # list to store the survival probabilities
  main_ow.prob<- main_ipw.prob <-gbm_ow.prob<- gbm_ipw.prob <- rf_ow.prob<- rf_ipw.prob <-lasso_ow.prob<- lasso_ipw.prob <-list() 
  rmst_ow.main<- rmst_ipw.main <- rmst_ow.gbm <- rmst_ipw.gbm<- rmst_ow.rf<-rmst_ipw.rf <-rmst_ow.lasso<-rmst_ipw.lasso <-rmst_ow.true<-rmst_ipw.true<-list()
  
  # read in the true parameter values
  idx <- ceiling(j/2)

  shr_ipw2.cox <- TRUEALL[[idx]]$shr_ipw2.cox
  shr_ow2.cox <- TRUEALL[[idx]]$shr_ow2.cox
  p_RATE2<- TRUEALL[[idx]]$p_RATE2
  p_RATO2 <-  TRUEALL[[idx]]$p_RATO2
  
  GASD_OW<- GASD_IPW<- NULL
  
  # simulation repetitions
  for (i in 1: nsim){
    sample_data <- data1[sample(1:dim(data1)[1], size= n ),c(-1:-4)]
    G = sample_data[,subgroup]
    X = as.matrix(select(sample_data, X1:X10) )
    trueX = as.matrix(select(sample_data, X1:X8) )
    Z = sample_data[,"Z"]
    full.form <- as.formula(paste("Z~X", paste(paste(paste("G[,", 1:2, sep=""), "]*X", sep=""), collapse=" + "), sep="+"))
    true.form <-  as.formula(paste("Z~trueX", paste(paste(paste("G[,", 1:2, sep=""), "]*trueX", sep=""), collapse=" + "), sep="+"))
    
    # estimate ps
    # true
    true.eh <- PSmethod_sub(ps.formula=true.form, subgroup=NULL, method="glm", data=sample_data)
    
    # logistic- main
    main.eh <- PSmethod_sub(ps.formula=main.form, subgroup=NULL, method="glm", data=sample_data)
      
    # gbm
    gbm.eh <- PSmethod_sub(ps.formula=main.form, subgroup=NULL, method="gbm",n.trees=20000,interaction.depth = 4,shrinkage=0.0005,
                            data=sample_data)
    
    # rfs
    rf.eh <- PSmethod_sub(ps.formula=main.form,data=sample_data, subgroup=NULL, method="RFs", num.trees = 1000, replace = T, 
                           seed = 1234 )
    adj_val = 0.00001
    if(range(rf.eh)[1]==0){rf.eh[which(rf.eh==0)]=adj_val}
    if(range(rf.eh)[2]==1){rf.eh[which(rf.eh==1)]=1-adj_val}
    
    lasso.eh <- PSmethod_sub(ps.formula=full.form, subgroup=NULL, method="pLASSO", data=sample_data)
    
    
    # check balance: initiate storage matrix
    gasd_ow <- array(NA, dim=c(10, 5, nmod))
    gasd_ipw <- array(NA, dim=c(10,5, nmod))
    
    covM<-sample_data[, paste("X",1:10,sep="")]
    
    gasd_ow[, , 1] <- SumStat_sub(zname = "Z", ps.estimate =main.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("overlap"))$ASD
    gasd_ow[, , 2] <- SumStat_sub(zname = "Z", ps.estimate =rf.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("overlap"))$ASD
    gasd_ow[, , 3] <- SumStat_sub(zname = "Z", ps.estimate =gbm.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("overlap"))$ASD
    gasd_ow[, , 4] <- SumStat_sub(zname = "Z", ps.estimate =lasso.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("overlap"))$ASD
    gasd_ow[, , 5] <- SumStat_sub(zname = "Z", ps.estimate =true.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("overlap"))$ASD
  
    gasd_ipw[, ,1] <- SumStat_sub(zname = "Z", ps.estimate =main.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("IPW"))$ASD
    gasd_ipw[, ,2] <- SumStat_sub(zname = "Z", ps.estimate =rf.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("IPW"))$ASD
    gasd_ipw[, ,3] <- SumStat_sub(zname = "Z", ps.estimate =gbm.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("IPW"))$ASD
    gasd_ipw[, ,4] <- SumStat_sub(zname = "Z", ps.estimate =lasso.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("IPW"))$ASD
    gasd_ipw[, ,5] <- SumStat_sub(zname = "Z", ps.estimate =true.eh ,trtgrp=1, subgroup = subgroup, data=sample_data, covM=covM, weight=c("IPW"))$ASD
    
    
    GASD_OW[[i]] <- gasd_ow
    GASD_IPW[[i]] <- gasd_ipw
      
    sample_formula =Surv(Y, Delta ) ~ Z
    
    # 1- estimated S-HR - IPW main
    shr_ipw <- smhr(formula=sample_formula, ps=main.eh, mweight="IPW",data= sample_data, subgroup)
    EST_MHR_IPW[,1,i] <- shr_ipw[,1]
    SE_MHR_IPW[,1,i] <- shr_ipw[,2]
    COVER_MHR_IPW[,1,i] <- (shr_ipw2.cox<EST_MHR_IPW[,1,i]+qnorm(0.975)*SE_MHR_IPW[,1,i])&(shr_ipw2.cox>EST_MHR_IPW[,1,i]-qnorm(0.975)*SE_MHR_IPW[,1,i])
    
    
    # 1- estimated S-HR - OW main
    shr_ow <- smhr(formula=sample_formula, ps=main.eh, mweight="overlap", data=sample_data, subgroup)
    EST_MHR_OW[,1,i] <- shr_ow[,1]
    SE_MHR_OW[,1,i] <- shr_ow[,2]
    COVER_MHR_OW[,1,i] <- (shr_ow2.cox<EST_MHR_OW[,1,i]+qnorm(0.975)*SE_MHR_OW[,1,i])&(shr_ow2.cox>EST_MHR_OW[,1,i]-qnorm(0.975)*SE_MHR_OW[,1,i])
    
    
    # 2- estimated S-HR - IPW RFs
    shr_ipw <- smhr(formula=sample_formula, ps=rf.eh,mweight="IPW",data= sample_data, subgroup)
    EST_MHR_IPW[,2,i] <- shr_ipw[,1]
    SE_MHR_IPW[,2,i] <- shr_ipw[,2]
    COVER_MHR_IPW[,2,i] <- (shr_ipw2.cox<EST_MHR_IPW[,2,i]+qnorm(0.975)*SE_MHR_IPW[,2,i])&(shr_ipw2.cox>EST_MHR_IPW[,2,i]-qnorm(0.975)*SE_MHR_IPW[,2,i])
    
    
    # 2- estimated S-HR - OW RFs
    shr_ow <- smhr(formula=sample_formula,ps=rf.eh,mweight="overlap", data=sample_data, subgroup)
    EST_MHR_OW[,2,i] <- shr_ow[,1]
    SE_MHR_OW[,2,i] <- shr_ow[,2]
    COVER_MHR_OW[,2,i] <- (shr_ow2.cox<EST_MHR_OW[,2,i]+qnorm(0.975)*SE_MHR_OW[,2,i])&(shr_ow2.cox>EST_MHR_OW[,2,i]-qnorm(0.975)*SE_MHR_OW[,2,i])
    
    # 3- estimated S-HR - IPW GBM
    shr_ipw <- smhr(formula=sample_formula, ps=gbm.eh,mweight="IPW",data= sample_data, subgroup)
    EST_MHR_IPW[,3,i] <- shr_ipw[,1]
    SE_MHR_IPW[,3,i] <- shr_ipw[,2]
    COVER_MHR_IPW[,3,i] <- (shr_ipw2.cox<EST_MHR_IPW[,3,i]+qnorm(0.975)*SE_MHR_IPW[,3,i])&(shr_ipw2.cox>EST_MHR_IPW[,3,i]-qnorm(0.975)*SE_MHR_IPW[,3,i])
    
    
    # 3- estimated S-HR - OW GBM
    shr_ow <- smhr(formula=sample_formula,ps=gbm.eh,mweight="overlap", data=sample_data, subgroup)
    EST_MHR_OW[,3,i] <- shr_ow[,1]
    SE_MHR_OW[,3,i] <- shr_ow[,2]
    COVER_MHR_OW[,3,i] <- (shr_ow2.cox<EST_MHR_OW[,3,i]+qnorm(0.975)*SE_MHR_OW[,3,i])&(shr_ow2.cox>EST_MHR_OW[,3,i]-qnorm(0.975)*SE_MHR_OW[,3,i])
    
    # 4- estimated S-HR - IPW pLASSO
    shr_ipw <- smhr(formula=sample_formula, ps=lasso.eh,mweight="IPW",data= sample_data, subgroup)
    EST_MHR_IPW[,4,i] <- shr_ipw[,1]
    SE_MHR_IPW[,4,i] <- shr_ipw[,2]
    COVER_MHR_IPW[,4,i] <- (shr_ipw2.cox<EST_MHR_IPW[,4,i]+qnorm(0.975)*SE_MHR_IPW[,4,i])&(shr_ipw2.cox>EST_MHR_IPW[,4,i]-qnorm(0.975)*SE_MHR_IPW[,4,i])
    
    
    # 4- estimated S-HR - OW pLASSO
    shr_ow <- smhr(formula=sample_formula,ps=lasso.eh,mweight="overlap", data=sample_data, subgroup)
    EST_MHR_OW[,4,i] <- shr_ow[,1]
    SE_MHR_OW[,4,i] <- shr_ow[,2]
    COVER_MHR_OW[,4,i] <- (shr_ow2.cox<EST_MHR_OW[,4,i]+qnorm(0.975)*SE_MHR_OW[,4,i])&(shr_ow2.cox>EST_MHR_OW[,4,i]-qnorm(0.975)*SE_MHR_OW[,4,i])
    
    # 5- estimated S-HR - IPW true
    shr_ipw <- smhr(formula=sample_formula, ps=true.eh, mweight="IPW",data= sample_data, subgroup)
    EST_MHR_IPW[,5,i] <- shr_ipw[,1]
    SE_MHR_IPW[,5,i] <- shr_ipw[,2]
    COVER_MHR_IPW[,5,i] <- (shr_ipw2.cox<EST_MHR_IPW[,5,i]+qnorm(0.975)*SE_MHR_IPW[,5,i])&(shr_ipw2.cox>EST_MHR_IPW[,5,i]-qnorm(0.975)*SE_MHR_IPW[,5,i])
    
    # 5- estimated S-HR - OW true
    shr_ow <- smhr(formula=sample_formula, ps=true.eh, mweight="overlap", data=sample_data, subgroup)
    EST_MHR_OW[,5,i] <- shr_ow[,1]
    SE_MHR_OW[,5,i] <- shr_ow[,2]
    COVER_MHR_OW[,5,i] <- (shr_ow2.cox<EST_MHR_OW[,5,i]+qnorm(0.975)*SE_MHR_OW[,5,i])&(shr_ow2.cox>EST_MHR_OW[,5,i]-qnorm(0.975)*SE_MHR_OW[,5,i])
    
    
    formula =Surv(Y, Delta ) ~ Z
  
    # 1- estimate subgroup RMST -main
    rmst_ow.main[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=main.eh, weight="overlap", subgroup=G, tau=365, alpha=.05, 
                  xaxismin=0, xaxismax=365)
    rmst_ipw.main[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=main.eh, weight="IPW", subgroup=G, tau=365, alpha=.05, 
                  xaxismin=0, xaxismax=365)
    
    # 2- estimate subgroup RMST -RFs
    rmst_ow.rf[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=rf.eh, weight="overlap", subgroup=G, tau=365, alpha=.05, 
                                xaxismin=0, xaxismax=365)
    rmst_ipw.rf[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=rf.eh, weight="IPW", subgroup=G, tau=365, alpha=.05, 
                                 xaxismin=0, xaxismax=365)
    
    # 3- estimate subgroup RMST -GBM
    rmst_ow.gbm[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=gbm.eh, weight="overlap", subgroup=G, tau=365, alpha=.05, 
                              xaxismin=0, xaxismax=365)
    rmst_ipw.gbm[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=gbm.eh, weight="IPW", subgroup=G, tau=365, alpha=.05, 
                               xaxismin=0, xaxismax=365)
    
    # 4- estimate subgroup RMST -pLASSO
    rmst_ow.lasso[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=lasso.eh, weight="overlap", subgroup=G, tau=365, alpha=.05, 
                              xaxismin=0, xaxismax=365)
    rmst_ipw.lasso[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=lasso.eh, weight="IPW", subgroup=G, tau=365, alpha=.05, 
                               xaxismin=0, xaxismax=365)
    
    # 5- estimate subgroup RMST - true
    rmst_ow.true[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=true.eh, weight="overlap", subgroup=G, tau=365, alpha=.05, 
                                      xaxismin=0, xaxismax=365)
    rmst_ipw.true[[i]] = akm_rmst_sub(time=sample_data$Y, status=sample_data$Delta, group=as.factor(sample_data$Z), ps=true.eh, weight="IPW", subgroup=G, tau=365, alpha=.05, 
                                      xaxismin=0, xaxismax=365)
    }
 
  
  GASD_OW <- apply(simplify2array(GASD_OW), 1:3, mean)
  GASD_IPW <- apply(simplify2array(GASD_IPW), 1:3, mean)
  
  MAIN_OW_RMST = apply(simplify2array(rmst_ow.main), 1:2, mean)
  MAIN_IPW_RMST = apply(simplify2array(rmst_ipw.main), 1:2, mean)

  RF_OW_RMST = apply(simplify2array(rmst_ow.rf), 1:2, mean)
  RF_IPW_RMST = apply(simplify2array(rmst_ipw.rf), 1:2, mean)
  
  GBM_OW_RMST = apply(simplify2array(rmst_ow.gbm), 1:2, mean)
  GBM_IPW_RMST = apply(simplify2array(rmst_ipw.gbm), 1:2, mean)
  
  LASSO_OW_RMST = apply(simplify2array(rmst_ow.lasso), 1:2, mean)
  LASSO_IPW_RMST = apply(simplify2array(rmst_ipw.lasso), 1:2, mean)
  
  TRUE_OW_RMST = apply(simplify2array(rmst_ow.true), 1:2, mean)
  TRUE_IPW_RMST = apply(simplify2array(rmst_ipw.true), 1:2, mean)
  
  MAIN_OW_RMST_SD = apply(simplify2array(rmst_ow.main), 1:2, sd)
  MAIN_IPW_RMST_SD = apply(simplify2array(rmst_ipw.main), 1:2, sd)
  
  RF_OW_RMST_SD = apply(simplify2array(rmst_ow.rf), 1:2, sd)
  RF_IPW_RMST_SD = apply(simplify2array(rmst_ipw.rf), 1:2, sd)
  
  GBM_OW_RMST_SD = apply(simplify2array(rmst_ow.gbm), 1:2, sd)
  GBM_IPW_RMST_SD = apply(simplify2array(rmst_ipw.gbm), 1:2, sd)
  
  LASSO_OW_RMST_SD = apply(simplify2array(rmst_ow.lasso), 1:2, sd)
  LASSO_IPW_RMST_SD = apply(simplify2array(rmst_ipw.lasso), 1:2, sd)
  
  TRUE_OW_RMST_SD = apply(simplify2array(rmst_ow.true), 1:2, sd)
  TRUE_IPW_RMST_SD = apply(simplify2array(rmst_ipw.true), 1:2, sd)
  
  
  
  RMST_OW_BIAS =sweep(cbind( MAIN_OW_RMST[1,],RF_OW_RMST[1,],GBM_OW_RMST[1,],LASSO_OW_RMST[1,],TRUE_OW_RMST[1,]),1,p_RATO2)
  RMST_IPW_BIAS =sweep(cbind( MAIN_IPW_RMST[1,],RF_IPW_RMST[1,],GBM_IPW_RMST[1,],LASSO_IPW_RMST[1,],TRUE_IPW_RMST[1,]),1,p_RATE2)

  RMST_OW_MCSD =cbind( MAIN_OW_RMST_SD[1,],RF_OW_RMST_SD[1,],GBM_OW_RMST_SD[1,],LASSO_OW_RMST_SD[1,],TRUE_OW_RMST_SD[1,])
  RMST_IPW_MCSD =cbind( MAIN_IPW_RMST_SD[1,],RF_IPW_RMST_SD[1,],GBM_IPW_RMST_SD[1,],LASSO_IPW_RMST_SD[1,],TRUE_IPW_RMST_SD[1,])
  
  
  tmp.list= lapply(1:nsim,function(x){(cbind(rmst_ipw.main[[x]][1,], rmst_ipw.rf[[x]][1,],rmst_ipw.gbm[[x]][1,],rmst_ipw.lasso[[x]][1,],rmst_ipw.true[[x]][1,])-as.matrix(p_RATE2)[,rep(1,5)])^2})
  RMST_IPW_RMSE= sqrt(apply(simplify2array(tmp.list), 1:2, mean))
  tmp.list2= lapply(1:nsim,function(x){(cbind(rmst_ow.main[[x]][1,], rmst_ow.rf[[x]][1,],rmst_ow.gbm[[x]][1,],rmst_ow.lasso[[x]][1,],rmst_ow.true[[x]][1,])-as.matrix(p_RATO2)[,rep(1,5)])^2})
  RMST_OW_RMSE= sqrt(apply(simplify2array(tmp.list2), 1:2, mean))
  
  cov.list= lapply(1:nsim,function(x){(cbind(rmst_ipw.main[[x]][3,], rmst_ipw.rf[[x]][3,],rmst_ipw.gbm[[x]][3,],rmst_ipw.lasso[[x]][3,],rmst_ipw.true[[x]][3,])<as.matrix(p_RATE2)[,rep(1,5)]) &  (cbind(rmst_ipw.main[[x]][4,], rmst_ipw.rf[[x]][4,],rmst_ipw.gbm[[x]][4,],rmst_ipw.lasso[[x]][4,],rmst_ipw.true[[x]][4,])>as.matrix(p_RATE2)[,rep(1,5)]) })
  RMST_IPW_CRATE = apply(simplify2array(cov.list), 1:2, mean)
  cov.list2= lapply(1:nsim,function(x){(cbind(rmst_ow.main[[x]][3,], rmst_ow.rf[[x]][3,],rmst_ow.gbm[[x]][3,],rmst_ow.lasso[[x]][3,],rmst_ow.true[[x]][3,])<as.matrix(p_RATO2)[,rep(1,5)]) &  (cbind(rmst_ow.main[[x]][4,], rmst_ow.rf[[x]][4,],rmst_ow.gbm[[x]][4,],rmst_ow.lasso[[x]][4,],rmst_ow.true[[x]][4,])>as.matrix(p_RATO2)[,rep(1,5)]) })
  RMST_OW_CRATE = apply(simplify2array(cov.list2), 1:2, mean)
  

  BIAS_IPW = abs(sweep(apply(EST_MHR_IPW, 1:2, mean), 1, shr_ipw2.cox))
  RMSE_IPW = sqrt(apply(sweep(EST_MHR_IPW, 1:2, matrix(rep(shr_ipw2.cox, 5), ncol=5))^2, 1:2, mean))
  MC_SD_IPW = apply(EST_MHR_IPW, 1:2, sd)
  EST_SD_IPW = apply(SE_MHR_IPW, 1:2, mean)
  CRATE_IPW = apply(COVER_MHR_IPW, 1:2, mean)
  
  BIAS_OW = abs(sweep(apply(EST_MHR_OW, 1:2, mean), 1, shr_ow2.cox))
  RMSE_OW = sqrt(apply(sweep(EST_MHR_OW, 1:2, matrix(rep(shr_ow2.cox, 5), ncol=5))^2, 1:2, mean))
  MC_SD_OW = apply(EST_MHR_OW, 1:2, sd)
  EST_SD_OW = apply(SE_MHR_OW, 1:2, mean)
  CRATE_OW = apply(COVER_MHR_OW, 1:2, mean)
  
  colnames(BIAS_IPW) <- colnames(BIAS_OW) <-colnames(RMSE_IPW) <-colnames(RMSE_OW) <-
    colnames(MC_SD_IPW) <-colnames(MC_SD_OW) <-colnames(EST_SD_IPW) <-colnames(EST_SD_OW) <-
    colnames(CRATE_IPW) <-colnames(CRATE_OW) <- colnames(RMST_IPW_RMSE) <- colnames(RMST_OW_RMSE) <-
    colnames(RMST_IPW_CRATE) <-colnames(RMST_OW_CRATE)<- colnames(RMST_IPW_BIAS)<-colnames(RMST_OW_BIAS)<-
    colnames(GASD_OW) <- colnames(GASD_IPW) <-
    c("Main","RFs","GBM", "pLASSO","True")
  
  
  rownames(BIAS_IPW) <- rownames(BIAS_OW) <-rownames(RMSE_IPW) <-rownames(RMSE_OW) <-
    rownames(MC_SD_IPW) <-rownames(MC_SD_OW) <-rownames(EST_SD_IPW) <-rownames(EST_SD_OW) <-
    rownames(CRATE_IPW) <-rownames(CRATE_OW) <- c("Overall","G1=0", "G1=1", "G2=0", "G2=1")
  

  all_results[[j]] <-list(BIAS_IPW=BIAS_IPW,BIAS_OW=BIAS_OW,RMSE_IPW=RMSE_IPW,RMSE_OW=RMSE_OW,
                          MC_SD_IPW=MC_SD_IPW, MC_SD_OW=MC_SD_OW, EST_SD_IPW=EST_SD_IPW,EST_SD_OW=EST_SD_OW,
                          CRATE_IPW=CRATE_IPW,CRATE_OW=CRATE_OW, RMST_OW_BIAS=RMST_OW_BIAS,RMST_IPW_BIAS=RMST_IPW_BIAS,
                          RMST_IPW_RMSE=RMST_IPW_RMSE,RMST_OW_RMSE=RMST_OW_RMSE,
                          RMST_IPW_CRATE=RMST_IPW_CRATE, RMST_OW_CRATE=RMST_OW_CRATE,
                          RMST_OW_MCSD=RMST_OW_MCSD,RMST_IPW_MCSD=RMST_IPW_MCSD, GASD_OW=GASD_OW, GASD_IPW=GASD_IPW
                          )
  
  
}

save(all_results, file=paste("all_simu_results",Sys.Date(),".Rdata",sep=""))
