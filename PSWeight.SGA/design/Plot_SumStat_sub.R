library(ggplot2)
library(plyr)

#print function

plot.SumStat_sub <- function(x,varlist=NULL,base=FALSE,plotsub=FALSE,...){

  vif<-x$vif
  ASD<-x$ASD
  subgroup<-x$subgroup
  
  if(base){
    vif<-x$vif_bs
    ASD<-x$ASD_bs
  }
  nsubg<-x$nsubg
  VIFW <- round(vif,2)
  
  
  if(!is.null(varlist)){
    ASD<-ASD[varlist,]
  }else{
    ASD<-ASD[!grepl('\\:',rownames(ASD)),]
    #print not subgroup as confounder
    if(!plotsub){
      namelist<-c()
      for (i in subgroup){
        namelist<-c(namelist,which(i==rownames(ASD)))
      }
      namelist<-unique(namelist)
      ASD<-ASD[-namelist,]
    }
  }
  

  # connect-S plot
  nv<-length(colnames(ASD))
  ncov<-length(rownames(ASD))
  mydata <- data.frame(ASMD = c(ASD),
                       Subgroups = rep(colnames(ASD), each=nrow(ASD)),
                       VarName = rep(rownames(ASD), nv)) 
  
  mydata$VarName <- factor(mydata$VarName, levels=unique(mydata$VarName))
  mydata$Subgroups <- factor(mydata$Subgroups, levels=unique(mydata$Subgroups))
  
  toplot<- mutate(mydata, lcolor = ifelse(ASMD <=0.1 |is.na(ASMD), "0",
                                          ifelse((ASMD <=0.15), "1",
                                                 ifelse((ASMD <=0.20), "2","3"))))
  toplot$lcolor <- factor(toplot$lcolor, levels=c("0","1","2","3"))
  
  group.colors <- c("0" = "White", "1" = "grey85", "2" ="grey50","3"="black")
  
  g1 <- ggplot(toplot, aes(x =VarName , y =Subgroups , fill=lcolor)) +
    geom_dotplot(binaxis = "y",binwidth = 0.75, stackdir = "center")+labs(x="Covariate", y="Subgroup"  )+
    scale_fill_manual(values=group.colors,name="ASMD",  labels=c("<0.1", "0.1-0.15", "0.15-0.20", ">0.20"), drop=FALSE)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text( angle=-45,vjust = 0, hjust=0.3),axis.text = element_text( size = 14 ),
          axis.title = element_text( size = 20, face = "bold" ),legend.position="bottom",legend.title = element_text(size=18),
          legend.text=element_text(size=18), legend.key.size = unit(2,"line"))+
    annotate("text", x = rep(ncov+1,nv), y = seq(1, nv,1), label = format(nsubg,zero.print = T),size = 5)+
    annotate("text", rep(ncov+2,nv), y = seq(1, nv,1), label = format(VIFW,zero.print = T),size = 5)+
    annotate("text", x = c(ncov+1,ncov+2), y = rep(nv+1,2), label = c("Size", "VIF"),size = 5)+
    expand_limits(x= c(0, length(levels(toplot$VarName)) + 3),y= c(0, length(levels(toplot$Subgroups)) + 1.5) )

  splot <- g1+ theme(panel.border = element_rect(linetype = "solid", fill = NA))
  
  splot
}


