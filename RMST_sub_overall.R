# --- RMST Using Adjusted KM ---
# Time is the time to event
# Status is 0 if censored, 1 if event
# Group should be a factor treatment variable
# ps is the estimated propensity score

# Weight is a character specifying the type of weight to be used, "IPW" specifies the inverse probability of treatment weights for estimating the
# average treatment effect among the combined population. "treated" specifies the weights for estimating the average treatment effect among the treated.
#"overlap" specifies the (generalized) overlap weights for estimating the average treatment effect among the overlap population, or population at clinical
# equipoise. "matching" specifies the matching weights for estimating the average treatment effect among the matched population (ATM). "entropy" specifies
# the entropy weights for the average treatment effect of entropy weighted population (ATEN). Default is "overlap".

# subgroup is a matrix containing subgroup indicators or a vector of string indicating subgroup names
# Tau is a user-specified truncation point. 
# If not specified, the default will be the minimum of the each groups' last event time 

akm_rmst_sub <- function(time, status, group, ps, weight="overlap", subgroup, data=NULL, tau=NULL, alpha=.05, 
                     xaxismin=0, xaxismax=max(time)){
  
  if(sum(time<0)>0){print("Error: times must be positive.")
  }else{
    if(sum(ps<=0)>0){print("Error: ps must be greater than 0.")
    }else{
      if(class(subgroup)!="character" & sum(status!=0 & status!=1)>0){print("Error: status must be a vector of 0s and/or 1s.")
      }else{
        if (class(subgroup)!="character"){
        # create subgroup
        zname <- "group"
        #set ordered group
        facz <- group}else{
          zname <- group
          facz<-as.factor(unlist(data[zname]))
          ncate<-nlevels(facz) #number of categories
        }
        ncate<-nlevels(facz)
        
        #set the treatment label
        dic<-levels(facz)
        znum<-as.numeric(facz) #numeric z
        #z<-znum_long -1 #z to fit the model
        z<- znum -1 #z to fit the model
        
        # tilting function
        trt <- 2
        ftilt<-tiltbin(weight = weight)
        
        #evaluated tilting function
        tilt.h<-ftilt(c(ps[,trt]))
        
        allwt<-(1/ps)*tilt.h
        n <- length(z)
        wt<-rep(0,n)
        for(i in 1:trt){
          wt[znum==i]<-allwt[znum==i,i]
        } 
        if (class(subgroup)!="character"){
            data <- data.frame(time, status, group, wt, subgroup)
            data <- data[!is.na(data$group) & !is.na(data$time),]
        }else{
              
            data <-cbind( select(data, all_of(time), all_of(status), all_of(group), all_of(subgroup) ), wt)
            data <- data[!is.na(data[,group]) & !is.na(data[,time]),]
            data[,group] <- as.factor(data[,group])
            }
        
        #data <- data[order(group),] 
        ######## Prepare subgroup matrix for estimation#############################################################################
        library(gridExtra)
        if( class(subgroup)=="character"){
          submatrix<-as.matrix(data[,subgroup,drop=FALSE])
          
        }else{
        #submatrix<-as.matrix(subgroup)
        submatrix<-as.matrix(subgroup)}       
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
        
        plot = list()
        RESULTS = NULL
        for(g in 1:ncol(submatrix_level)){
          #pname <- paste0("p",g)
          subtmp<-c(submatrix_level[,g])
          #normalize weights within subgroup
          wt_n <- wt*z*subtmp/sum(wt*z*subtmp)+ wt*(1-z)*subtmp/sum(wt*(1-z)*subtmp)
          
          subdata <- data[subtmp==1, ]
          subwt <- wt[subtmp==1]
          subdata$wt <- subwt
          #--- If tau not specified, use minimum tau from all groups ---
          j=length(unique(subdata[,group]))
          
          if(is.null(tau)){
            taui = rep(999, j)
            for (i in (1:j)){
              groupval <- (levels(subdata[,group])[i])
              dat_group <- subdata[which(subdata[,group]==(groupval)),]
              taui[i] <- max(dat_group[,time][dat_group$status==1])
            }
            tau <- min(taui)
          }
          
          #--- Calculate AKM RMST in each group ---
          rmst <- rep(999, length(1:j))
          groupval <- rep(999, length(1:j))
          rmst_var <- rep(999, length(1:j))
          rmst_se <- rep(999, length(1:j))
          # p <- plot(NULL, xlim=c(xaxismin, xaxismax), ylim=c(0,1), xlab='Time',ylab='Adjusted Survival Probability')
          # title(main=paste(weight,'Adjusted Kaplan-Meier'))
          # 
          for (i in 1:j){
            groupval[i] <- levels(subdata[,group])[i]
            dat_group <- subdata[which(subdata[,group]==(groupval[i])),]
            #--- AKM ---
            # Based on 'adjusted.KM' function from {IPWsurvival} package
            # Author: F. Le Borgne and Y. Foucher
            tj <- c(0,sort(unique(dat_group[,time][ dat_group[,status] ==1])))
            dj <- sapply(tj, function(x){sum(dat_group$wt[dat_group[,time]==x & dat_group[,status]==1])})
            yj <- sapply(tj, function(x){sum(dat_group$wt[dat_group[,time]>=x])})
            st <- cumprod(1-(dj/yj))
            m <- sapply(tj, function(x){sum((dat_group$wt[dat_group[,time]>=x])^2)})
            mj <- ((yj^2)/m)
            #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
            ft <- data.frame(tj, yj, dj, st, i, mj)
            
            #--- RMST ---
            # Based on 'rmst1 function' from {survRM2} package
            # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
            rtime <- ft$tj<=tau
            tj_r <- sort(c(ft$tj[rtime],tau))
            st_r <- ft$st[rtime]
            yj_r <- ft$yj[rtime]
            dj_r <- ft$dj[rtime]
            time_diff <- diff(c(0, tj_r))
            areas <- time_diff * c(1, st_r)
            rmst[i] <- sum(areas)
            
            #--- Variance ---
            mj_r <- ft$mj[rtime]
            var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
            #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
            var_r <- c(var_r,0)
            rmst_var[i] <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
            rmst_se[i] <- sqrt(rmst_var[i])
            
            #--- Plot AKM ---
            # lines(ft$tj, ft$st,type="s", col=(i+2), lwd=2)
          }
  
  
      #--- Add legend and tau to plot ---
      # abline(v=tau, col=1, lty=3, lwd=2)
      # legend('bottomleft', paste("Group", groupval), lty=rep(1, j), lwd=rep(2, j), col=3:(j+2), 
      #        cex=.75, bty ="n", inset = c(0, 0))
      # mtext(paste0(level_name[g]), side = 3)
      # 
      #--- Compare RMST between groups and compile output---
      #results <- data.frame(groupval,rmst,rmst_var,rmst_se,tau)
      
      # Based on 'rmst2 function' from {survRM2} package
      # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
      
      #--- RMST Difference ---
      rmst_diff <- rmst[2] - rmst[1]
      rmst_diff_se <- sqrt(sum(rmst_var))
      rmst_diff_low <- rmst_diff - qnorm(1-alpha/2)*rmst_diff_se
      rmst_diff_upp <- rmst_diff + qnorm(1-alpha/2)*rmst_diff_se
      results <-round(c(rmst_diff,rmst_diff_se,rmst_diff_low, rmst_diff_upp),3)
      
      # cat("\n\n\n")
      # cat(paste('RMST calculated up to tau =',round(results$tau[1],3)))
      # cat("\n\n\n")
      # 
      # cat ("Restricted Mean Survival Time (RMST) per Group \n\n")
      RESULTS <- cbind(RESULTS, results)
       # print(round(results[c(2,4)],3))
      # cat("\n\n")
      # 
      # cat ("Restricted Mean Survival Time (RMST) Differences \n\n")
      # colnames(output_diff) <- c("Groups", "Est.", "SE", "CIL", "CIU", "p")
      # rownames(output_diff) <- c(output_diff$Groups)
      # print(round(output_diff[c(2,3,4,5,6)],3))
      # cat("\n\n")
      # 
        }
  colnames(RESULTS) <-  colnames(submatrix_level)
  rownames(RESULTS) <- c("RMSTDIFF", "SE", "LL", "UL")
  
      }
    }
  }
  return(RESULTS)}