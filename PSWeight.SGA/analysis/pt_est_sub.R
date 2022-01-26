# The point estimating function for binary case in subgroups
## a function to calculate the subgroup
pt_est_sub<-function(psest,y,z,submatrix_level,ftilt,newlevel){
  mu.est<-c()
  for(g in 1:ncol(submatrix_level)){
    subtmp<-c(submatrix_level[,g])
    # point estimate
    
    mu1est <- sum(z*y*ftilt(psest)/psest*subtmp) / sum(z*ftilt(psest)/psest*subtmp)
    mu0est <- sum((1-z)*y*ftilt(psest)/(1-psest)*subtmp) / sum((1-z)*ftilt(psest)/(1-psest)*subtmp)
    mu.est<-rbind(mu.est,c(mu0est,mu1est))
  }
  
  rownames(mu.est)<-colnames(submatrix_level)
  colnames(mu.est)<-newlevel
  
  return(mu.est)
}

