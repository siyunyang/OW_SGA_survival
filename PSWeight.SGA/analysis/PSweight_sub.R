# source("tilt.R")
# source('pt_est_sub.R')
# source('PSmethod_sub.R')

script.dir <- dirname(sys.frame(1)$ofile)
source(paste0(script.dir,"/tilt.R"))
source(paste0(script.dir,"/pt_est_sub.R"))
source(paste0(script.dir,"/PSmethod_sub.R"))

PSweight_sub <- function(ps.formula=NULL,ps.estimate=NULL,subgroup=NULL,trtgrp=NULL,zname=NULL,yname,data,R=50,weight="overlap",method="glm",opt=FALSE, ps.control_sub=list()){
   
  ## Arguments:
  # ps.formula(formula for non-subgroup effect only) an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the propensity score model to be fitted. Additional details of model specification are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
  # ps.estimate an optional matrix or data frame containing estimated (generalized) propensity scores for
  #' each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the
  # subgroup name of subgroup variables of data
  # trtgrp an optional character defining the "treated" population for estimating the average treatment effect among the treated (ATT). Only necessary if \code{weight = "treated"}. This option can also be used to specify the treatment (in a two-treatment setting) when a vector argument is supplied for \code{ps.estimate}. Default value is the last group in the alphebatic order.
  # zname an optional character specifying the name of the treatment variable in \code{data}.
  # yname an optional character specifying the name of the outcome variable in \code{data}.
  # data an optional data frame containing the variables in the propensity score model. If not found in data, the variables are taken from \code{environment(formula)}.
  # R number of bootstrap samples, default is 50
  # weight a character or vector of characters including the types of weights to be used.
  # method 'glm' 'gbm' ...
  # opt a logical of whether to use optimal stopping for gbm
  
  ## Returns: 
  
  #  propensity: estimated propensity
  #  muhat: muhat: point estimate
  #  muboot: the bootstrap result
  #  group: the original group level
  #  trtgrp: treatment group if using att
  
  #####################################################################
  
  #####NO subgroup specified, abort the program########################################################################
    if(is.null(subgroup)){
      stop("No subgroup specified")
    }
  
  
    #####Extract zname, treatment label and relabel the treatment to 0/1##################################################
    if(is.null(ps.estimate)){
      ps.formula<-as.formula(ps.formula)
      zname<-all.vars(ps.formula)[1]
      
      #option to include all covariates
      if(all.vars(ps.formula)[2]=='.'){
        #exclude some variables
        exclude_var<-c(subgroup,zname,yname)
        all_var<-colnames(data)
        include_var<-setdiff(all_var,exclude_var)
        
        formulatmp<- paste(zname,'~')
        formulatmp<- paste(formulatmp,paste(include_var,collapse = "+"),'+')
        formulatmp<- paste(formulatmp,paste(subgroup,collapse = "+"))
        
        for(i in subgroup){
          #interaction terms
          inter_tmp<-paste0(include_var, ':', i)
          inter_tmp<-paste("+",paste(inter_tmp,collapse = " + "))
          formulatmp<-paste(formulatmp,inter_tmp)
        }
        ps.formula<-as.formula(formulatmp)
      }
    }
    
    #set ordered group to match group later
    facz<-as.factor(unlist(data[zname]))
    
    
    #set the treatment label
    oldlevel<-levels(facz)
    newlevel<-oldlevel
    
    z<-as.numeric(facz)-1
    n<-length(z) #total obs
    
    #set group for ATT
    if(weight=="treated"){
      if(is.null(trtgrp)){
        trtgrp<-oldlevel[2]
      }else{
        newlevel<-rev(unique(c(trtgrp,oldlevel)))
        facz<-factor(unlist(data[zname]),levels = newlevel)
        z<-as.numeric(facz)-1
      }
    }
    
    matchlevel<-match(oldlevel,newlevel) #match level for ATT
    
    #reassign value as numeric
    y<-unlist(data[yname])
    data[zname]<-z
    data_p<-data[,colnames(data)!=yname] #data without outcome

  

    ######## Prepare subgroup matrix for estimation#############################################################################
    submatrix<-as.matrix(data[,subgroup,drop=FALSE])
    submatrix<-cbind(1,submatrix) #include the overall effect
    colnames(submatrix)[1]<-"Overall"
    
    submatrix_level<-NULL
    level_name<-c()
    sub_n<-c()
    
    # matrix for subgroup
    for(r in 1:ncol(submatrix)){ 
      leveltmp <- sort(unique(submatrix[,r]))
      sub_n<-c(sub_n,length(leveltmp))
      for(g in leveltmp){
        submatrix_level<-cbind(submatrix_level,c(submatrix[,r]==g)*1)
        level_name<-c(level_name,paste0(colnames(submatrix)[r],"=",g))
      }
    }
    colnames(submatrix_level)<-level_name
    colnames(submatrix_level)[1]<-"Overall"
    
  
    ####### Estimate ps ########################################################################################################

    if(is.null(ps.estimate)){

      e.h <- do.call(PSmethod_sub,c(list(ps.formula = ps.formula, method=method, data=data_p, weight=weight, opt=opt),ps.control_sub))
      e.h<-c(e.h[,2])
    }else{
      #the name for the propensity score
      if(length(ps.estimate)==n){
        e.h<-c(ps.estimate)
      }else if(!setequal(colnames(ps.estimate),newlevel)){
        ps.estimate<-as.matrix(ps.estimate)
        e.h<-c(ps.estimate[,2])
        warning("wrong column name set for ps.estimate, treatment set as: ",newlevel[1], " , ", newlevel[2])
      }else{
        ps.estimate<-ps.estimate[,match(newlevel,colnames(ps.estimate))]
        e.h<-c(ps.estimate[,2])
      }
    }
    
    ###### Genearate point estimate and evaluate uncertainty through bootstrap ##############################################
    if(weight=="entropy"){
      e.h<-pmax(e.h,1e-6)
      e.h<-pmin(e.h,1-1e-6)
    }
    
    #tilting function
    ftilt<-tiltbin(weight = weight)
    
    #calculate point
    muhat<-pt_est_sub(psest=e.h,y=y,z=z,submatrix_level=submatrix_level,ftilt=ftilt,newlevel=newlevel)  
    

    #bootstrap
    mu1.boot<-NULL
    mu0.boot<-NULL
    for(i in 1:R){
      if(i %% 50==0){
        message("bootstrap ", i, " samples")
      }
      # estimate ps
      samp.b<-sample(n,n,replace = TRUE)
      data.b<-data_p[samp.b,]
      y.b<-y[samp.b]
      z.b<-z[samp.b]
      submatrix_level.b<-submatrix_level[samp.b,]      
      if(is.null(ps.estimate)&method=='glm'){
        fit.b <- glm(formula = ps.formula, data=data.b, family = binomial(link = "logit"))
        e.b <- as.numeric(fit.b$fitted.values)
      }else{
        e.b<-e.h[samp.b]
      }
     
      
      #calculate point
      muhat.b<-pt_est_sub(psest=e.b,y=y.b,z=z.b,submatrix_level=submatrix_level.b,ftilt=ftilt,newlevel=newlevel)  
      mu1.boot<-rbind(mu1.boot,c(muhat.b[,2]))
      mu0.boot<-rbind(mu0.boot,c(muhat.b[,1]))
    }
    rownames(mu1.boot)<-NULL
    rownames(mu0.boot)<-NULL
    
    muboot<-list(mu0.boot,mu1.boot)
    names(muboot)<-newlevel
    e.h<-cbind(1-e.h,e.h)
    colnames(e.h)<-newlevel
    names(sub_n)<-c("Overall",subgroup)
    
    #match back to original levels
    e.h<-e.h[,matchlevel]
    muhat<-muhat[,matchlevel]
    muboot<-muboot[matchlevel]
    
    out<-list(propensity=e.h, muhat=muhat, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp,sub_n=sub_n)
    class(out)<-'PSweight_s'
    out
  }
  
    
    
    


