#### Caution: This code takes a long time to run. 
### The true causal effects are saved in the "_true_val_censor_100_2021-09-13_.Rdata" file if you want to save time.
#### library

load.lib<-c("survival", "survminer","survRM2","dplyr", "tidyr","mvtnorm",  "parallel","gridExtra")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
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




  ## Arguments:
  # n - number of individuals
  # kappa = controls strength of interaction terms in the PS model
  # taus = conditional treatment effect
  #####################################################################


  
  all_params <- expand.grid(ns = c(1000000),
                            gammas= c(0.5,1),
                            kappas = c(-1/2, -1/4),
                            taus = c(0, -log(2)),
                            PHs=c(TRUE, FALSE)
  )
  nsim <- 100; 
  TRUEALL <- NULL
  for(j in 1:nrow(all_params)){
   n =  all_params[j, "ns"]; kappa=all_params[j, "kappas"];
   gamma =  all_params[j, "gammas"];  tau =  all_params[j, "taus"];
   PH =all_params[j, "PHs"]
   print(c(n,kappa,gamma, tau, PH))
   
   if (j<=4){p_RATE2 <-p_RATO2 <- shr_ipw2.cox <-shr_ow2.cox <- rep(0, 5)}else{
      shr_ipw.cox <- shr_ipw2.cox <- shr_ow.cox<- shr_ow2.cox <- p_RATE <-p_RATE2<- p_RATO<- p_RATO2 <- matrix(NA,5, nsim) 
      ## 1) Generate data
      for (i in 1:nsim){
        set.seed(i*10+ 123)
        
        data1 = sim_data(n, tau,gamma, kappa, PH=PH, censoring=TRUE)
        data1 = as.data.frame(data1)
        data1[,"ID"] =seq(1,nrow(data1))
        
        data1 <- data1[ , -which(names(data1) %in% c("Y","Delta"))]
        
        # subgroup variable
        subgroup= paste0("G",1:2)
        G = data1[,subgroup]
        
        Z = as.matrix(data1$Z) #treatment
        #set ordered group
        facz<-as.factor(Z)
        ncate<-nlevels(facz)
        
        #set the treatment label
        dic<-levels(facz)
        znum<-as.numeric(facz) #numeric z
        # ture ps
        e.h <- data1[,"pZ"]
        e.h <- cbind(1-e.h, e.h)
        
        main.form <- as.formula(paste("Z~",paste(paste("X",1:10,sep="", collapse="+"),paste("G",1:2,sep="", collapse="+"),sep="+")))
        
        
        # truncated time
        truncate <- 365
        
        # superpopulation with two potential outcomes
        #wide to long
        data_long <- gather(data1, status, time, Survival_time0:Survival_time1, factor_key=TRUE)
        
        data_long$status <- ifelse(data_long$status=="Survival_time0", "Control", "Treated") #treatment
        data_long$Y =pmin(data_long$time, data_long$Censor_time)
        data_long$Delta =ifelse( data_long$time < data_long$Censor_time, 1, 0)
   
        
       
        e.h_long <- rbind(e.h, e.h)
        Z_long = as.matrix(data_long$status) #treatment
        #set ordered group
        facz<-as.factor(Z_long)
        #set the treatment label
        znum_long <- as.numeric(facz) #numeric z
        z_long<-znum_long -1 #z to fit the model
        
        data_long_sorted <- data_long[ order(data_long$ID),]
   
        
        # true overall MHR - approach 1: using true propensity score weights
        data1$time <- data1$Z *data1$Survival_time1+(1-data1$Z) *data1$Survival_time0
        
       
        
        # true overall MHR - approach 2: using two potential outcomes, weighted by tilting function h
        # true S-HR - IPW
        formula =Surv(Y, Delta) ~ status
     
        # for IPW, everyone has the same weight, so I manually set the ps=0.5
        shr_ipw2.cox[,i] <- smhr(formula, ps=matrix(0.5, nrow = nrow(data_long), ncol = 2), mweight="IPW",data= data_long, subgroup)[,1]

        # true S-HR - OW
        # for OW, the weight should be h=e(1-e), corresponding to the overlap population
        # so I manually set the ps=1-e(1-e), resulting in w=1-ps= e(1-e).
        shr_ow2.cox[,i] <- smhr(formula,ps=1-e.h_long*(1-e.h_long),mweight="overlap", data=data_long, subgroup)[,1]

        
       
        # approach2:
        colnames(e.h_long) <- c("Control","Treated")
         
        p_R <- akm_rmst_sub(time=data_long$Y, status=data_long$Delta, group=as.factor(data_long$status), ps=matrix(0.5,nrow(data_long),2), weight="IPW", subgroup=data_long[,subgroup], tau=365, alpha=.05,
                     xaxismin=0, xaxismax=365)
        p_RATE2[,i] <- p_R[1,]
        
        
        p_RA <-akm_rmst_sub(time=data_long$Y, status=data_long$Delta, group=as.factor(data_long$status), ps=1-e.h_long*(1-e.h_long), weight="overlap", subgroup=data_long[,subgroup], tau=365, alpha=.05,
                     xaxismin=0, xaxismax=365)
        
        p_RATO2[,i] <- p_RA[1,]

      }
   
   shr_ipw2.cox=round(apply(shr_ipw2.cox, 1, mean),2)
   shr_ow2.cox=round(apply(shr_ow2.cox, 1, mean),2) 
   p_RATE2=round(apply(p_RATE2, 1, mean),2) 
   p_RATO2=round(apply(p_RATO2, 1, mean) ,2) 
   }
   TRUEALL[[j]] <- list( shr_ipw2.cox =shr_ipw2.cox,  shr_ow2.cox=shr_ow2.cox , p_RATE2=p_RATE2, p_RATO2=p_RATO2)

  }
  
  
  save(TRUEALL, file=paste("_true_val_censor",nsim,Sys.Date(),".Rdata",sep="_"))
  
  
 