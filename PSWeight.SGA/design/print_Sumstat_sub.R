print.SumStat_sub<-function(x,...) {
  #weights<-paste(x$group,collapse = ', ')
  L<-length(x)
  if (names(x)[L]=="trim"){
    cat(as.numeric(rowSums(x$trim))[1],"cases trimmed, ",as.numeric(rowSums(x$trim))[2],"cases remained","\n")
    cat("\n")
    cat("trimmed result by trt group:","\n")
    cat(capture.output(x$trim),sep="\n")
    cat("\n")
  }
  weights<-names(x$ps.weights)[-(1:2)]
  if ("treated" %in% weights){
    cat(paste('trt group for PS model is: ',x$trtgrp), "\n")
  }
  cat(paste('weights estimated for: ',paste(x$ps.weights,collapse=" ")), "\n")
  cat(paste('subgroup and size: '), "\n")
  print(x$nsubg)
}

