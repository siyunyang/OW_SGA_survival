print.PSweightsum_s <- function(x,...) {
  group<-x$group
  trtgrp<-x$trtgrp
  contrast<-x$contrast
  estimates<-x$estimates
  het_eval<-x$het_eval
  
  prtgroup<-paste(x$group,collapse = ', ')
  cat('Original group value: ',prtgroup, "\n")
  if(!is.null(x$trtgrp)){
    cat(paste('Treatment group value: ',trtgrp),'\n')
  }
  cat('\n')
  cat('contrast','\n')
  print(contrast)
  cat('\n')
  cat('Use bootstrap sample for inference','\n')
  cat('\n')
  cat('Result by subgroups: ','\n')
  print(estimates)
  if(!is.null(het_eval)){
    cat('Test for hetrogenous effect: ','\n')
    print(het_eval)
  }
}

