#########################################################
#### #1) Data genearating function ######################
#########################################################
library(mvtnorm)
sim_data <- function(n=1000000, tau,gamma, kappa, PH=TRUE, censoring=TRUE){

X <-round( rmvnorm(n, mean = rep(0, 10), diag(rep(1, 10))),2)
colnames(X) <- paste("X",1:10,sep="")
G <- t(apply( X[,7:8], 1, function(x) x>qnorm(c(0.8,0.6))))
colnames(G) <- paste("G",1:2,sep="")
#ps model for interactions between X and R
ps.form <- as.formula(paste("~X", paste(paste(paste("G[,", 1:2, sep=""), "]*X", sep=""), collapse=" + "), sep="+"))
prop_X <- model.matrix(ps.form) [,c(-12,-13)]
#Fix the colnames
colnames(prop_X) <- c("Intercept", paste("X", 1:10, sep=""),
                      c(sapply(paste("G", 1:2, sep=""), function(x) paste(x,paste("X", 1:10 , sep="_"), sep="*"))))
alpha <- c(log(1.25), log(1.25), log(1.5), log(1.5), log(1.75), log(1.75), log(2.5), log(2.5),0,0)*gamma
alpha1 <-alpha2 <- kappa*alpha
trt_prop <- 0.5
#itcpt <- log(trt_prop/(1-trt_prop))-c(alpha1, alpha2)*c(0.2,0.4,0.2,0.4)
itcpt <- 0
theta <- c(itcpt, alpha, alpha1, alpha2)
out.form <- as.formula("~X")
out_X <- model.matrix(out.form)
beta <- c(0,0,0,log(1.5), log(1.5), log(1.75), log(1.75), log(2), log(2), log(1.25),log(1.25))
LP0 <- out_X %*% matrix(beta)
LP1 <- out_X %*% matrix(beta)+ tau
lin_termZ <- prop_X %*% matrix(theta)
pZ <- exp(lin_termZ)/(1+exp(lin_termZ))
Z <- rbinom(length(pZ), 1, pZ)


#proportional hazard
if(PH){
random_simu <- runif(n)
v =2
lambda = 0.00002
Survival_time0 = round((-log(random_simu) / (lambda * exp(LP0))) ^ (1 / v))
Survival_time1 = round((-log(random_simu) / (lambda * exp(LP1))) ^ (1 / v))
Survival_time = Survival_time1*Z + Survival_time0*(1-Z)
}else{
# # AFT but not PH; log-normal
# random_simu = rnorm(n, 0,sd = 0.81)
# Survival_time0 =round( exp(LP0 + random_simu))
# Survival_time1 = round( exp(LP1 + random_simu))
# Survival_time = Survival_time1*Z + Survival_time0*(1-Z)
  random_simu <- runif(n)
  #lambda = 0.00002
  Survival_time0 = round((-log(random_simu) / (0.0000003  * exp(LP0))) ^ (1 / 3))
  Survival_time1 = round((-log(random_simu) / (0.00002 * exp(LP1))) ^ (1 / 2))
  Survival_time = Survival_time1*Z + Survival_time0*(1-Z)
  
}
# Censor time
if (censoring){
Censor_time =round( runif(n, min = 0, max = 445))
# mean(Survival_time > Censor_time)
}else{Censor_time = max(Survival_time)+100}

# Observed survival
Censored = Survival_time > Censor_time
Y = pmin(Survival_time, Censor_time)
Delta = 1 - Censored
tab <- cbind(pZ, Survival_time0, Survival_time1, Censor_time, Y, Delta, Z, X, G)
colnames(tab)[c(1:3,5:6)] <- c("pZ", "Survival_time0", "Survival_time1","Y","Delta")
return(tab)

}