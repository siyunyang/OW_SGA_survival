print.PSweight_s <- function(x,...) {
  group<-paste(x$group,collapse = ', ')
  cat('Original group value: ',group, "\n")
  if(!is.null(x$trtgrp)){
    cat(paste('Treatment group value: ',x$trtgrp),'\n')
  }
  cat('\n')
  cat('Point estimate by subgroup:','\n')
  ptest<-round(x$muhat,4)
  print(ptest)
}
