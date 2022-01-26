# The tilting function for binary treatment
tiltbin<-function(weight="overlap"){
  if(weight=="overlap"){
    return(function(x){x*(1-x)})
  }else if(weight=="IPW"){
    return(function(x){rep(1,length(x))})
  }else if(weight=="matching"){
    return(function(x){pmin(x,1-x)})
  }else if(weight=="entropy"){
    return(function(x){-(x*log(x)+(1-x)*log(1-x))})
  }else if(weight=="treated"){
    return(function(x){x})
  }
}


# S-MHR
smhr <- function(formula, ps, mweight, data, subgroup, method=c("efron","breslow","exact","discrete"), bootstrap=FALSE, R=50)
{ formula <- as.formula(formula)
  zname <- all.vars(formula)[3]
  Z = as.matrix(data[,zname]) #treatment
  #set ordered group
  facz<-as.factor(Z)
  ncate<-nlevels(facz)

  #set the treatment label
  dic<-levels(facz)
  znum<-as.numeric(facz) #numeric z
  z<-znum -1 #z to fit the model
  
  # tilting function
  trt <- 2
  ftilt<-tiltbin(weight = mweight)
  
  #evaluated tilting function
  tilt.h <- ftilt(c(ps[,trt]))
  
  allwt<-(1/ps)*tilt.h
  n <- dim(data)[1]
  wt<-rep(0,n)
  for(i in 1:trt){
    wt[znum==i]<-allwt[znum==i,i]
  }
  ######## Prepare subgroup matrix for estimation #############################################################################
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
  
  mu.est<-c()
  for(g in 1:ncol(submatrix_level)){
    print(paste0("estimate subgroup:",level_name[g]))
    subtmp<-c(submatrix_level[,g])
    subdata <- data[subtmp==1, ]
    
    #normalize weights within subgroup
    wt_n <- wt*z*subtmp/sum(wt*z*subtmp)+ wt*(1-z)*subtmp/sum(wt*(1-z)*subtmp)
    subwt <- wt_n[subtmp==1]
    subdata$subwt <- subwt
    if(method %in% c("efron","breslow","exact")){
    # point estimate
    res.cox <- coxph(formula=formula,data = subdata, weights=subwt, robust=TRUE, method = method)
    coef <- round(res.cox$coefficients , 2)
    }else{
      res.cox <- clogit(formula=formula, data = subdata, weights= subwt, method="efron")
      coef <- round(res.cox$coefficients , 2)
      se <- round(summary(res.cox)$coefficients[4] ,2)
    }
    if(bootstrap){
      #bootstrap
      mu.boot<-NULL
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }
        #print(i)
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data[samp.b,]
         
        subtmp.b <- submatrix_level[samp.b, g]  
        subdata.b <- as.data.frame( data.b[subtmp.b == 1, ])
        #print(dim(subdata.b))
        
        #e.b<-ps[samp.b]
        wt.b <- wt[samp.b]
        #print(length(wt.b))
        #wt_n.b <- wt.b*z*subtmp/sum(wt*z*subtmp)+ wt*(1-z)*subtmp/sum(wt*(1-z)*subtmp)
        subwt.b <- c(wt.b[subtmp.b==1])
        #print(length(subwt.b))
        subdata.b$subwt.b <- subwt.b
        #calculate point
        if(method %in% c("efron","breslow","exact")){
        res.cox.b <- coxph(formula=formula, data = subdata.b, weights=subwt.b, method = method )
        muhat.b<-round(res.cox.b$coefficients ,2)}else{
          res.cox.b <-  clogit(formula=formula, data = subdata.b, weights=subwt.b, method="efron")  
        }
        mu.boot<-c(mu.boot,muhat.b)
      }
      se <- sd(mu.boot)
      
    }else{se <- round(summary(res.cox)$coefficients[4] ,2)}
    
    mu.est<-rbind(mu.est,c(coef, se))
    
    
  }
  
  
  rownames(mu.est)<-colnames(submatrix_level)
  colnames(mu.est)<-c("coef", "se")
  if (bootstrap){return(  list(mu.boot=mu.boot, mu.est=mu.est))}else{return(mu.est)}
  
  summary(warnings())
}

# S-KM plot
skmplot <- function(formula, ps, mweight, data, subgroup, xupper, conf.int)
{  
    zname <- all.vars(formula)[3]
    Z = as.matrix(data[,zname]) #treatment
    #set ordered group
    facz<-as.factor(Z)
    ncate<-nlevels(facz)
    
    #set the treatment label
    dic<-levels(facz)
    znum<-as.numeric(facz) #numeric z
    #z<-znum_long -1 #z to fit the model
    z<- znum -1 #z to fit the model
    
    # tilting function
    trt <- 2
    ftilt<-tiltbin(weight = mweight)
    
    #evaluated tilting function
    tilt.h<-ftilt(c(ps[,trt]))
    
    allwt<-(1/ps)*tilt.h
    n <- dim(data)[1]
    wt<-rep(0,n)
    for(i in 1:trt){
      wt[znum==i]<-allwt[znum==i,i]
    } 
  ######## Prepare subgroup matrix for estimation#############################################################################
  library(gridExtra)
  submatrix<-as.matrix(data[,subgroup,drop=FALSE])
  
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
  plot = list()
  PROB = NULL
  for(g in 1:ncol(submatrix_level)){
    #pname <- paste0("p",g)
    subtmp<-c(submatrix_level[,g])
    #normalize weights within subgroup
    wt_n <- wt*z*subtmp/sum(wt*z*subtmp)+ wt*(1-z)*subtmp/sum(wt*(1-z)*subtmp)
    
    subdata <- data[subtmp==1, ]
    subwt <- wt_n[subtmp==1]
    f1 <- surv_fit(formula, data = subdata, weights=subwt)
    fit.prob <- summary(f1, times = seq(0, xupper), extend=T)$surv
    
    #RMST
    #rmst2(subdata$Y, subdata$Delta, subdata$Z, tau= 365)
    
    PROB <- cbind(PROB, fit.prob)
    ggsurv <- ggsurvplot(f1, data = subdata, xlim= c(0,xupper), title=paste0(level_name[g]), legend.title = "",legend = c(0.8, 0.98),surv.median.line = "hv",
                legend.labs = c("Control", "Treated"),conf.int = conf.int, palette = "grey", font.title=c(12, "plain","black"),
                font.x=c(10, "plain","black"), font.y=c(10, "plain","black"),font.legend=c(10, "plain","black"),font.tickslab=c(10, "plain","black"))
    tmedian <-surv_median(f1)[,'median'] 
    names(tmedian) <- surv_median(f1)[,'strata']
    
    plot[[g]] <- ggsurv$plot+ ggplot2::annotate("text", x =tmedian, y = 0, label = paste0(tmedian), size = 3.5)

     #print(plot[[g]])
    #assign(pname, g)
    # p<-ggsurvplot(f1, data = subdata, xlim= c(0,xupper), title=paste0(level_name[g]), legend.title = "",
    # legend.labs = c("Control", "Treated"),conf.int = conf.int, palette = "grey")
    # print(p)          
  }
  colnames(PROB) <- level_name
  p <- do.call(grid.arrange, c(plot, nrow=4))
  return (list(p=p, PROB=PROB))
}

