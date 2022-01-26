summary.PSweight_s<-function(object,contrast=NULL,type='DIF',het=FALSE,...){
  
  muhat<-object$muhat
  muboot<-object$muboot
  muboot0<-muboot[[1]]
  muboot1<-muboot[[2]]
  trtgrp<-object$trtgrp
  group<-object$group
  sub_n<-object$sub_n
  
  sub_group<-rownames(muhat)
  
  
  
  #groups
  ngrp<-length(group)
  
  
  #transform contrast into matrix
  if(is.vector(contrast)){
    contrast<-t(contrast)
  }
  
  #error message for wrong contrast
  if(!is.null(contrast)){
    if(dim(contrast)[2]!=ngrp){
      cat('Contract length not equal to treatment groups, please check','\n')
      cat('\n')
      contrast<-NULL
    }
    cat('\n')
  }
  
  
  
  
  #if all contrast
  if(is.null(contrast)){
    ncst<-ngrp*(ngrp-1)/2
    contrasttmp<-matrix(0,ngrp-1,ngrp)
    contrasttmp[,1]<--1
    contrasttmp[1:(ngrp-1),2:ngrp]<-diag(1,ngrp-1,ngrp-1)
    contrast<-contrasttmp
    if(ngrp>2){
      for(i in 1:(ngrp-2)){
        ntmp<-nrow(contrasttmp)-1
        contrasttmp<-cbind(0,contrasttmp)[1:ntmp,1:ngrp]
        contrast<-rbind(contrast,contrasttmp)
      }
    }
  }
  
  rownames(contrast)<-paste('Contrast',1:nrow(contrast))
  colnames(contrast)<-group
  
  #define the bootstrap p value function
  pvalboot <- function(x,est) {
    esttmp<-abs(est)
    x<-abs(x-mean(x))
    return(mean(x>esttmp))
  }
  
 
  #calculate interval and test statistics for subgroup
  int_pval<-function(muboottmp,muhattmp,type=type,contrast=contrast){
    if(type=='DIF'){
      samp<-muboottmp%*%t(contrast)

      est.h<-c(contrast%*%muhattmp)
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025,na.rm=TRUE))
      ucl<-apply(samp,2,function(x) quantile(x,0.975,na.rm=TRUE))
      p.value<-c()
      dimcon<-dim(samp)[2]
      for(j in 1:dimcon){
        esttmp<-c(est.h)[j]
        p.value<-c(p.value,pvalboot(samp[,j],esttmp))
      }
      
    }else if(type=='RR'){
      samp<-log(muboottmp)%*%t(contrast)
      tranmuhat<-log(muhattmp)
      tranest<-c(contrast%*%tranmuhat)
      est.h<-tranest
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025))
      ucl<-apply(samp,2,function(x) quantile(x,0.975))
      p.value<-c()
      dimcon<-dim(samp)[2]
      for(j in 1:dimcon){
        esttmp<-c(est.h)[j]
        p.value<-c(p.value,pvalboot(samp[,j],esttmp))
      }
      
    }else if(type=='OR'){
      samp<-log(muboottmp/(1-muboottmp))%*%t(contrast)
      tranmuhat<-log(muhattmp/(1-muhattmp))
      tranest<-c(contrast%*%tranmuhat)
      est.h<-tranest
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025))
      ucl<-apply(samp,2,function(x) quantile(x,0.975))
      p.value<-c()
      dimcon<-dim(samp)[2]
      for(j in 1:dimcon){
        esttmp<-c(est.h)[j]
        p.value<-c(p.value,pvalboot(samp[,j],esttmp))
      }
      
    }else{
      cat('type not found')
      cat('\n')
    }
    estimatestmp<-cbind(est.h,se.h,lcl,ucl,p.value)
    colnames(estimatestmp)<-c("Estimate","Std.Error","Lower.CL","Upper.CL","p.value")
    rownames(estimatestmp)<-rownames(contrast)
    return(estimatestmp)
    
  }
  
  estimates<-list()
  for(g in sub_group){
    muhattmp<-muhat[g,]
    muboottmp<-cbind(muboot0[,g],muboot1[,g])
    estimates[[g]]<-int_pval(muboottmp,muhattmp,type = type,contrast = contrast)
  }
  
  
  
  
  #test for hetrogenous effect
  het_eval<-NULL
  het_grp<-function(muboottmp,muhattmp){
    samp<-muboottmp%*%c(1,-1)
    est.h<-c(c(1,-1)%*%muhattmp)
    return(list(samp=samp,est.h=est.h))
  }
  
  if(het){
    het_eval<-list()

    #produce the point estimate and bootstrap sample for heterogeneous test
    het_est<-c()
    het_samp<-NULL
    for(g in sub_group[-1]){
      muhattmp<-muhat[g,]
      muboottmp<-cbind(muboot0[,g],muboot1[,g])
      het_tmp<-het_grp(muboottmp,muhattmp)
      het_est<-c(het_est,c(het_tmp$est.h))
      het_samp<-cbind(het_samp,het_tmp$samp)
    }
    
    
    count=1
    for(g1 in 1:(length(sub_n)-1)){ #remove the overall
      ngrp<-sub_n[g1+1]
      nametmp<-names(sub_n)[g1+1]
      idxtmp<-count:(count+ngrp-1)  
      hetmutmp<-het_est[idxtmp]
      hetsamptmp<-het_samp[,idxtmp]
      subnametmp<-sub_group[1+idxtmp]
      
      count<-count+ngrp
      
      #specify the contrast for hetrogenous test
      ncst<-ngrp*(ngrp-1)/2
      contrasttmp<-matrix(0,ngrp-1,ngrp)
      contrasttmp[,1]<--1
      contrasttmp[1:(ngrp-1),2:ngrp]<-diag(1,ngrp-1,ngrp-1)
      contrast_het<-contrasttmp
      if(ngrp>2){
        for(i in 1:(ngrp-2)){
          ntmp<-nrow(contrasttmp)-1
          contrasttmp<-cbind(0,contrasttmp)[1:ntmp,1:ngrp]
          contrast_het<-rbind(contrast_het,contrasttmp)
        }
      }
    
      het_eval[[g1]]<-int_pval(hetsamptmp,hetmutmp,type = 'DIF',contrast = contrast_het)
      
      ones<-apply(contrast_het,1,function(x) which(x==1))
      nones<-apply(contrast_het,1,function(x) which(x==-1))
      rownames(het_eval[[g1]])<- paste(subnametmp[ones],'-',subnametmp[nones],": ")
    }
    names(het_eval)<-names(sub_n)[-1]
  }
  
  
  out<-list(estimates=estimates,contrast=contrast,group=group,trtgrp=trtgrp,het_eval=het_eval)
  class(out)<-'PSweightsum_s'
  out
}