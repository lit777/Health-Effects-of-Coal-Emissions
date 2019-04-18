#------ Load the main dataset
load("Dataset_2005_2012_noMedicare.RData")
load("adj.mat.RData")
load("matrix_median.RData")

library(rootSolve)
library(matrixcalc)
library(data.table)
library(mvtnorm)
library(Rcpp)
library(mnormt)
library(gtools)
library(numDeriv)
library(corpcor)

sourceCpp("CARBayes.cpp")
source("commonfunctions.R")

# ------ Define 7 regions by state names
Northeast = c("ME", "NH", "VT", "NY", "PA", "DE", "NJ", "MD", "DC", "VA", "MA", "CT", "RI")
IndustrialMidwest = c("WV", "OH", "KY", "IN", "IL", "WI", "MI")
Southeast = c("FL", "GA", "SC", "NC", "TN", "AL", "MS", "AR","LA")
UpperMidwest = c("MN", "IA", "MO", "KS", "NE", "SD", "ND")


#------ Narrowing down to 3 EPA regions
Master <- subset(Master, State.zip %in% c(IndustrialMidwest, Northeast, Southeast))
Master <- subset(Master, State.zip %in% c(IndustrialMidwest, Northeast, Southeast))
Master <- subset(Master, !(State.zip %in% c("ME", "FL", "AR", "LA", "MS","VT","NH","WI")))
Master <- subset(Master, !is.na(as.numeric(Master$Ter)))
Master <- subset(Master, !is.na(PM2005))
Master <- subset(Master, !is.na(PctHighSchool))
Master <- subset(Master, TotPop > 100)
Master <- subset(Master, Tot_den_for_death_MA_FFS2005 > 100)
Master <- subset(Master, Person_year_FFS2005 > 100)

Master$X <- ifelse(Master$exp2005 < summary(Master$exp2005)[4], 1,0 )
Master$TTEMP5 <- Master$TTEMP5 - 273.15


#------ Preprocessing using Propensity Scores
ps <- glm(X~TTEMP5+
            smokerate2005+
            mean_age2005+
            log(MedianHHInc+1)+
            Female_rate2005+
            Black_rate2005+
            PctHighSchool+
            PctUrban+
            PctHisp+
            log(TotPop+1), family=binomial, data=Master)
pair <- cbind(ps$fitted.values, Master)
setnames(pair, "V1", "ps")
range1 <- summary(subset(pair, X==1)$ps)[c(1,6)]
range0 <- summary(subset(pair, X==0)$ps)[c(1,6)]

range <- c(max(range1[1],range0[1]),min(range1[2],range0[2]))

pair <- subset(pair, ps < range[2] & ps > range[1])



Data <- pair
m0 <- log(as.numeric(subset(Data, X==0)$PM2005))
m1 <- log(as.numeric(subset(Data, X==1)$PM2005))
y00 <- subset(Data, X==0)$AMI2005   #<- AMI rates from the Medicare data
y11 <- subset(Data, X==1)$AMI2005   #<- AMI rates from the Medicare data

offset0 <- subset(Data, X==0)$Person_year_FFS2005
offset1 <- subset(Data, X==1)$Person_year_FFS2005

y0 <- (y00)/offset0 # <- Rate outcomes
y1 <- (y11)/offset1 # <- Rate outcomes
y0 <- log(y0+1)
y1 <- log(y1+1)

zip0 <- subset(Data, X==0)$zip
zip1 <- subset(Data, X==1)$zip

x0 <- as.matrix(subset(Data, X==0)[,list(1,
                                          TTEMP5,
                                          smokerate2005,
                                          mean_age2005,
                                          MedianHHInc/10000,
                                          Female_rate2005,
                                          Black_rate2005,
                                          PctHighSchool,
                                          PctUrban,
                                          PctHisp,
                                          TotPop/10000),ps])
x1 <- as.matrix(subset(Data, X==1)[,list(1,
                                          TTEMP5,
                                          smokerate2005,
                                          mean_age2005,
                                          MedianHHInc/10000,
                                          Female_rate2005,
                                          Black_rate2005,
                                          PctHighSchool,
                                          PctUrban,
                                          PctHisp,
                                          TotPop/10000),ps])



#------ Prepare the zipcode adjecent matrix
#adj.mat1 <- ifelse(adj.mat=="TRUE", 1L, 0L) 
#adj.mat1 <- adj.mat1[c(which(Data2$X==1),which(Data2$X==0)), c(which(Data2$X==1),which(Data2$X==0))]
#islands <- which(apply(adj.mat1, 1, sum) == 1)
#islands.zip <- c(subset(Data2, X==1)$zip,subset(Data2, X==0)$zip)[islands]
#W <- common.Wcheckformat(adj.mat1)

n0 <- dim(x0)[1]
n1 <- dim(x1)[1]
dim.cov<- dimx <-dim(x0)[2]
P <- 4 # number of marginal distributions


#------- Metropolis-hastings for the marginals
metropolis <- function(h, X, Y,R,invR, W,  beta0_pre, beta_pre, sigma_pre, alpha_pre, mu_beta0_pre,
                       sigma_beta0_pre, alpha_beta0_pre, alpha_sigma_pre, K, eps,del2, del3,
                       index, cov1, cov2, ran_pre,zz, tau2,proposal.sd.phi){
  YY <- Y; XX <- X;
  
  # Estimate of the intercept
  mu_pop <- 0
  
  # Define AR variables
  acc1 <- acc2 <- 1
  acc3 <- rep(1,dim(X)[2])
  
  acc6 <- rep(1,dim(X)[2])
  
  ###### Proposal distributions for hyper-parameters
  prop6 <- max(rgamma(1,del2*alpha_beta0_pre, del2),0.1^10) # alpha_beta0
  prop7 <- max(rgamma(1,del3*alpha_sigma_pre, del3),0.1^10) # alpha_sigma
  prop8 <- rnorm(1,mean=mu_beta0_pre, 0.1) # mu_beta0
  prop9 <- runif(1, 0, var(YY)) # sigma_beta0
  
  rat <- logprior(YY, XX, beta0_pre, beta_pre, sigma_pre,  mu_beta0=prop8, sigma_beta0=prop9, alpha_beta0=prop6, alpha_sigma=prop7, mu_pop=mu_pop)+
    dgamma(alpha_sigma_pre, del3*prop7, del3, log=TRUE)-
    logprior(YY, XX, beta0_pre, beta_pre, sigma_pre,  mu_beta0_pre, sigma_beta0_pre, alpha_beta0_pre, alpha_sigma_pre, mu_pop=mu_pop)-
    dgamma(prop7,del3*alpha_sigma_pre, del3, log=TRUE)
  
  if(log(runif(1))>rat){
    prop6 <- alpha_beta0_pre
    prop7 <- alpha_sigma_pre
    prop8 <- mu_beta0_pre
    prop9 <- sigma_beta0_pre
    acc1=0
  }
  

  prop3=beta_pre
  pprop3=beta_pre
  Cov2 = 2.38^2/((dim.cov)*2)*(cov2/2+diag(0.1^10,dim.cov))
  pprop3=(rmnorm(1,mean=beta_pre, Cov2))
  
  h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, ran_pre+XX%*%(prop3), sqrt(sigma_pre)))))
  h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, ran_pre+XX%*%(pprop3), sqrt(sigma_pre)))))
  
  rat <- logCopula(h,index,h.cur2,R,zz) + loglik(YY,XX,prop2,pprop3,sigma_pre, ran_pre) +
    logprior(YY,XX,prop2,pprop3,sigma_pre,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop=mu_pop) +
    dmnorm(prop3, (pprop3), Cov2,log=TRUE) -
    logCopula(h,index,h.cur1,R,zz) - loglik(YY,XX,prop2,prop3,sigma_pre, ran_pre) -
    logprior(YY,XX,prop2,prop3,sigma_pre,mu_beta0=prop8,sigma_beta0=prop9,alpha_beta0=prop6,alpha_sigma=prop7,mu_pop=mu_pop) -
    dmnorm((pprop3), prop3, Cov2,log=TRUE)
  if(is.na(rat)){
    pprop3=beta_pre
    acc6[1]=0
  }else{
    if (log(runif(1))>rat) {
      pprop3=beta_pre
      acc6[1]=0
    }else{prop3=pprop3}}
  
  
  prop4 <- pprop4 <- sigma_pre
  log_prop4 <- log_pprop4 <- log(sigma_pre)
  log_pprop4 <- rnorm(1, log(sigma_pre), 0.1)
  
  h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, ran_pre+XX%*%(prop3), sqrt(exp(log_prop4))))))
  h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, ran_pre+XX%*%(prop3), sqrt(exp(log_pprop4))))))
  
  rat <- logCopula(h,index,h.cur2,R,zz) + loglik(YY,XX,prop2,prop3,exp(log_pprop4), ran_pre) +
    dgamma(exp(log_pprop4), 5, scale=0.2, log=TRUE)+log(abs(exp(log_pprop4)))+
    dnorm(log_prop4, log_pprop4,0.1,log=TRUE)-
    logCopula(h,index,h.cur1,R,zz) - loglik(YY,XX,prop2,prop3,exp(log_prop4), ran_pre) -
    dgamma(exp(log_prop4), 5, scale=0.2, log=TRUE)-log(abs(exp(log_prop4)))-dnorm(log_pprop4, log_prop4,0.1,log=TRUE)
  
  if(log(runif(1))>rat){
    pprop4 <- sigma_pre
  }else{
    prop4 <- exp(log_pprop4)
  }
  
  # Spatial random effects
  beta.offset <- XX%*%(prop3)
  temp1 <- gaussiancarupdate(Wtriplet=W$W.triplet, Wbegfin=W$W.begfin, W$W.triplet.sum, 
                             nsites=length(YY), phi=ran_pre, tau2=tau2, y=YY, phi_tune=proposal.sd.phi, 
                             rho=rho, offset=beta.offset, h=h, R=R, invR = invR, index=index-1, sd=sqrt(prop4))
  phi <- temp1[[1]]
  phi[islands] <- phi[islands] - mean(phi[islands])
  
  temp2 <- quadform(W$W.triplet, W$W.triplet.sum, W$n.triplet, length(YY), phi, phi, rho)
  tau2.posterior.scale <- temp2 + prior.tau2[2] 
  tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
  
  # return acceptance indicators and new values
  return(c(prop2,prop3,prop4,prop7, phi, tau2,temp1[[2]]))
}



#---------- MCMC Execution
fit.y0 <- lm(y0~x0-1)
fit.y1 <- lm(y1~x1-1)
fit.m0 <- lm(m0~x0-1)
fit.m1 <- lm(m1~x1-1)

MCMC=30000 # <- the number of MCMC iterations
#---------- Setting Initials
para.y1=matrix(nrow=2,ncol=1+dim.cov+2+length(y0)+length(y1)+1+1)
para.y1[1,] <- para.y1[2,] <- c(1,coef(fit.y1),1,1, rep(0, length(y0)), rep(0, length(y1)), 0.01,1)
para.y0=matrix(nrow=2,ncol=1+dim.cov+2+length(y0)+length(y1)+1+1)
para.y0[1,] <- para.y0[2,] <- c(1,coef(fit.y0),1,1, rep(0, length(y0)), rep(0, length(y1)), 0.01,1)
para.m1 <- matrix(nrow = 2, ncol = 1+dim.cov+2+length(y0)+length(y1)+1+1) # Y(0)
para.m1[1,] <- para.m1[2,] <- c(1,coef(fit.m1),1,1, rep(0, length(y0)), rep(0, length(y1)), 0.01,1)
para.m0 <- matrix(nrow = 2, ncol = 1+dim.cov+2+length(y0)+length(y1)+1+1) # m1(0)
para.m0[1,] <- para.m0[2,] <- c(1,coef(fit.m0),1,1, rep(0, length(y0)), rep(0, length(y1)), 0.01,1)
para.C=list(prop1=c(0.1,0,0,0,0,0.1), rrho=rep(0,4))
para.C=list(para.C, para.C)
R=matrix(c(1,cor(y1,m1),0,0,cor(y1,m1),1,0,0,0,0,1,cor(y0,m0),0,0,cor(y0,m0),1),byrow=TRUE, nrow=4)
h=rmnorm(n0+n1, mean=rep(0,4),R)

cov2.m1 <- vcov(lm(m1~x1-1))
cov2.m0 <- vcov(lm(m0~x0-1))
cov2.y1 <- vcov(lm(y1~x1-1))
cov2.y0 <- vcov(lm(y0~x0-1))

y1.all <- c(y1,rep(round(mean(y1)),n0))
y0.all <- c(rep(round(mean(y0)),n1), y0)
m1.all <- c(m1,rep(mean(m1),n0))
m0.all <- c(rep(mean(m0),n1), m0)

# Indices of parameters
ind1 <- (K+1):(2*K) # intercept
ind2 <- (2*K+1):(2*K+dim.cov-1) # regression coefficient except the intercept
ind3 <- (2*K+dim.cov+1-1):(3*K+dim.cov-1) # variance
ind4 <- 3*K+dim.cov+1-1 # alpha
ind5 <- ind4+1 # alpha_beta0
ind6 <- ind5+1 # alpha_sigma
ind7 <- ind6+1 # mu_beta0
ind8 <- ind7+1 # sigma_beta0

# MCMC Quantities
prior.tau2 <- c((n1+n0)*5, 0.01) #The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape,scale) prior for tau2.
tau2.posterior.shape <- prior.tau2[1] + 0.5 * (n1+n0-length(islands)) 
rho <- 1
proposal.sd.phi.y1 <- 0.01
proposal.sd.phi.y0 <- 0.01
proposal.sd.phi.m1 <- 0.01
proposal.sd.phi.m0 <- 0.01

tt <- 0
ttt <- 2
# Placeholders for the quantities of interest
TE <- NIE <- NDE <- matrix(nrow=2000, ncol=n0+n1)
PS <- matrix(nrow=2000, ncol=4)
para.sub.y0 <- matrix(nrow=MCMC,ncol=dim.cov)
para.sub.y1 <- matrix(nrow=MCMC,ncol=dim.cov)
para.sub.m0 <- matrix(nrow=MCMC,ncol=dim.cov)
para.sub.m1 <- matrix(nrow=MCMC,ncol=dim.cov)

pb <- txtProgressBar(min = 0, max = MCMC, style = 3)
for (ttt in 3:MCMC){
  
  t <- 2
  invR <- solve(R)
  cov.y0 <- 1
  if(t > 200){
    cov2.y0 <- cov(para.sub.y0[3:(ttt-1),(1):(dim.cov)])
  }
  para.y0[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=y0.all, R=R, invR=invR, W=W, 
                            beta0_pre=para.y0[t-1,1], beta_pre=para.y0[t-1,(1+1):(1+dim.cov)], sigma_pre=para.y0[t-1,(1+dim.cov+1)],
                            alpha_beta0_pre=0, alpha_sigma_pre=para.y0[t-1,(1+dim.cov+2)],
                            mu_beta0_pre=0, sigma_beta0_pre=0, K=K, eps=0.1, del2=10, del3=30, index=3, cov1=cov.y0, cov2=cov2.y0,ran_pre=para.y0[t-1, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)], zz=0, tau2=para.y0[t-1, (1+dim.cov+2+n0+n1+1)],proposal.sd.phi=proposal.sd.phi.y0)
  
  y0.all[1:(n1)] <- rnorm(n1,para.y0[t, (1+dim.cov+2+1):(1+dim.cov+2+n1)]+x1%*%para.y0[t,(1+1):(1+dim.cov)], sqrt(para.y0[t,(1+dim.cov+1)]))
  h[,3] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(y0.all, para.y0[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+rbind(x1,x0)%*%para.y0[t,(1+1):(1+dim.cov)], sqrt(para.y0[t,(1+dim.cov+1)])))))
  
  
  cov.y1 <- 1
  if(ttt > 200){
    cov2.y1 <- cov(para.sub.y1[3:(ttt-1),(1):(dim.cov)])
  }
  para.y1[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=y1.all, R=R,invR=invR, W=W,
                            beta0_pre=para.y1[t-1,1], beta_pre=para.y1[t-1,(1+1):(1+dim.cov)], sigma_pre=para.y1[t-1,(1+dim.cov+1)],
                            alpha_beta0_pre=0, alpha_sigma_pre=para.y1[t-1,(1+dim.cov+2)],
                            mu_beta0_pre=0, sigma_beta0_pre=0, K=K, eps=0.1, 
                            del2=10, del3=30, index=1, cov1=cov.y1, cov2=cov2.y1,  ran_pre=para.y1[t-1, (1+dim.cov+2+1):(1+dim.cov+2+length(y0)+length(y1))], zz=1, tau2=para.y1[t-1, (1+dim.cov+2+length(y0)+length(y1)+1)],proposal.sd.phi=proposal.sd.phi.y1)
  
  y1.all[(n1+1):(n0+n1)] <- rnorm(n0,para.y1[t, (1+dim.cov+2+1+length(y1)):(1+dim.cov+2+n0+n1)]+x0%*%para.y1[t,(1+1):(1+dim.cov)], sqrt(para.y1[t,(1+dim.cov+1)]))
  h[,1] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(y1.all, para.y1[t, (1+dim.cov+2+1):(1+dim.cov+2+n1+n0)]+rbind(x1,x0)%*%para.y1[t,(1+1):(1+dim.cov)], sqrt(para.y1[t,(1+dim.cov+1)])))))
  
  
  cov.m0 <- 1
  if(ttt > 200){
    cov2.m0 <- cov(para.sub.m0[3:(ttt-1),(1):(dim.cov)])
  }
  para.m0[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0.all, R=R,invR=invR, W=W,
                            beta0_pre=para.m0[t-1,1], beta_pre=para.m0[t-1,(1+1):(1+dim.cov)], sigma_pre=para.m0[t-1,(1+dim.cov+1)],
                            alpha_beta0_pre=0, alpha_sigma_pre=para.m0[t-1,(1+dim.cov+2)],
                            mu_beta0_pre=0, sigma_beta0_pre=0, K=K, eps=0.1, del2=10, del3=30, index=4, cov1=cov.m0, cov2=cov2.m0,  ran_pre=para.m0[t-1, (1+dim.cov+2+1):(1+dim.cov+2+length(y0)+length(y1))], zz=0, tau2=para.m0[t-1, (1+dim.cov+2+n0+n1+1)],proposal.sd.phi=proposal.sd.phi.m0)
  
  m0.all[1:(n1)] <- rnorm(n1,para.m0[t, (1+dim.cov+2+1):(1+dim.cov+2+n1)]+x1%*%para.m0[t,(1+1):(1+dim.cov)], sqrt(para.m0[t,(1+dim.cov+1)]))
  h[,4] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(m0.all, para.m0[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+rbind(x1,x0)%*%para.m0[t,(1+1):(1+dim.cov)], sqrt(para.m0[t,(1+dim.cov+1)])))))
  
  
  cov.m1 <- 1
  if(ttt > 200){
    cov2.m1 <- cov(para.sub.m1[3:(ttt-1),(1):(dim.cov)])
  }
  
  para.m1[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1.all, R=R, invR= invR, W=W,
                            beta0_pre=para.m1[t-1,1], beta_pre=para.m1[t-1,(1+1):(1+dim.cov)], sigma_pre=para.m1[t-1,(1+dim.cov+1)],
                            alpha_beta0_pre=0, alpha_sigma_pre=para.m1[t-1,(1+dim.cov+2)],
                            mu_beta0_pre=0, sigma_beta0_pre=0, K=K, eps=0.1,
                            del2=10, del3=30, index=2, cov1=cov.m1, cov2=cov2.m1, ran_pre=para.m1[t-1, (1+dim.cov+2+1):(1+dim.cov+2+length(y0)+length(y1))], zz=1, tau2=para.m1[t-1, (1+dim.cov+2+n0+n1+1)],proposal.sd.phi=proposal.sd.phi.m1)
  
  m1.all[(n1+1):(n0+n1)] <- rnorm(n0,para.m1[t, (1+dim.cov+2+1+n1):(1+dim.cov+2+n0+n1)]+x0%*%para.m1[t,(1+1):(1+dim.cov)], sqrt(para.m1[t,(1+dim.cov+1)]))
  h[,2] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(m1.all, para.m1[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+rbind(x1,x0)%*%para.m1[t,(1+1):(1+dim.cov)], sqrt(para.m1[t,(1+dim.cov+1)])))))
  
  accept.y1 <- c(para.y1[t, (1+dim.cov+2+n0+n1+2)],n1+n0)
  accept.y0 <- c(para.y0[t, (1+dim.cov+2+n0+n1+2)],n1+n0)
  accept.m1 <- c(para.m1[t, (1+dim.cov+2+n0+n1+2)],n1+n0)
  accept.m0 <- c(para.y0[t, (1+dim.cov+2+n0+n1+2)],n1+n0)
  
  
  para.C[[t]]=metropolisC(h=h,rho=para.C[[t-1]]$prop1)
  Prop1 <- para.C[[t]]$prop1
  R = matrix(c(1,Prop1[1:3],Prop1[1],1,Prop1[4:5],Prop1[2],Prop1[4],1,Prop1[6],Prop1[3],Prop1[c(5,6)],1),4,4,byrow=TRUE)
  

  para.sub.y0[(ttt),(1):(dim.cov)] <- para.y0[t,(1+1):(1+dim.cov)]
  para.sub.y1[(ttt),(1):(dim.cov)] <- para.y1[t,(1+1):(1+dim.cov)]
  para.sub.m0[(ttt),(1):(dim.cov)] <- para.m0[t,(1+1):(1+dim.cov)]
  para.sub.m1[(ttt),(1):(dim.cov)] <- para.m1[t,(1+1):(1+dim.cov)]
  
  COR <- para.C[[t-1]] <- para.C[[t]]
  para.y0[t-1,] <- para.y0[t, ]
  para.y1[t-1,] <- para.y1[t, ]
  para.m1[t-1,] <- para.m1[t, ]
  para.m0[t-1,] <- para.m0[t, ]
  
  
  
  #  Sys.sleep(0.001)
  #  setTxtProgressBar(pb, t)
  #  print(ttt)
  
  
  #----- Conditional distribution of Y{1;M(1)}
  F.y1m1 <- function(o1){
    mm1 <- C.sample[o1, ]
    M1.dist <- qnorm(pmin(1-0.1^15,pmax(pnorm(mm1, mean = para.m1[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+s.covariate%*%para.m1[t,(1+1):(1+dim.cov)], sd = sqrt(para.m1[t-1,(1+dim.cov+1)])), 0.1^15)))
    con.joint <- pnorm(rnorm((n1+n0)*500, Cor01[1,2] %*% solve(Cor01[2,2]) %*% t(M1.dist), sd=sqrt(Cor01[1,1]-Cor01[1,2]%*%solve(Cor01[2,2])%*%t(Cor01[1,2]))))
    result <- qnorm(pmin(1-0.1^15, pmax(0.1^15, con.joint )),mean = para.y1[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+rbind(x1,x0)%*%para.y1[t,(1+1):(1+dim.cov)], sd = sqrt(para.y1[t-1,(1+dim.cov+1)]) )
    result <- matrix(exp(result)-1, nrow=500, byrow=T)
    return(result)
  }
  
  #----- Conditional distribution of Y_{0;M(0)}
  F.y0m0 <- function(o1){
    mm1 <- C.sample[o1, ]
    M1.dist <- qnorm(pmin(1-0.1^15,pmax(pnorm(mm1, mean = para.m0[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+s.covariate%*%para.m0[t,(1+1):(1+dim.cov)], sd = sqrt(para.m0[t-1,(1+dim.cov+1)])), 0.1^15)))
    con.joint <- pnorm(rnorm((n1+n0)*500,Cor01[3,4] %*% solve(Cor01[4,4]) %*% t(M1.dist), sd=sqrt(Cor01[3,3]-Cor01[3,4]%*%solve(Cor01[4,4])%*%t(Cor01[3,4]))))
    result <- qnorm(pmin(1-0.1^15, pmax(0.1^15, con.joint )),mean = para.y0[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+rbind(x1,x0)%*%para.y0[t,(1+1):(1+dim.cov)], sd = sqrt(para.y0[t-1,(1+dim.cov+1)]) )
    result <- matrix(exp(result)-1, nrow=500, byrow=T)
    return(result)
  }
  
  index <- seq(10001, MCMC, by =10) # <- burn-in: 10000; thining 10
  if(ttt %in% index){
    tt <- tt + 1
    s.covariate <- rbind(x1,x0)
    
    Cor01 <- matrix(c(1,COR$prop1[1:3],
                      COR$prop1[1],1,COR$prop1[4:5],
                      COR$prop1[2],COR$prop1[4],1,
                      COR$prop1[6],COR$prop1[3],
                      COR$prop1[c(5,6)],1),4,4,byrow=TRUE)
        
    s.m11 <- m1.all
    s.m10 <- m0.all
    C.sample <- t(cbind(s.m10,s.m11))
    
    #----- Values of potential outcomes
    y10 <- colMeans(F.y1m1(1))
    y11 <- colMeans(F.y1m1(2))
    y00 <- colMeans(F.y0m0(1))
    
    #----- Causal Mediation Effects
    TE[tt,] <- y11 - y00
    NIE[tt, ] <- y11 - y10
    NDE[tt, ] <- y10 - y00
    
    
    M1.dist <- qnorm(pmin(1-0.1^15,pmax(pnorm(s.m11, mean = para.m1[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+s.covariate%*%para.m1[t,(1+1):(1+dim.cov)], sd = sqrt(para.m1[t-1,(1+dim.cov+1)])), 0.1^15)))
    M0.dist <- qnorm(pmin(1-0.1^15,pmax(pnorm(s.m10, mean = para.m0[t, (1+dim.cov+2+1):(1+dim.cov+2+n0+n1)]+s.covariate%*%para.m0[t,(1+1):(1+dim.cov)], sd = sqrt(para.m0[t-1,(1+dim.cov+1)])), 0.1^15)))

    #----- Conditional distribution of Y_{1;M(1),M(0)}
    con.joint <- pnorm(rnorm(n1+n0,t(Cor01[1,c(2,4)] %*% solve(Cor01[c(2,4),c(2,4)]) %*% t(cbind(M1.dist,M0.dist))), sd=sqrt(Cor01[1,1]-Cor01[1,c(2,4)]%*%solve(Cor01[c(2,4),c(2,4)])%*%t(t(Cor01[1,c(2,4)])))))
    
    #----- Conditional distribution of Y_{0;M(1),M(0)}
    con.joint <- pnorm(rnorm(n1+n0,t(Cor01[3,c(2,4)] %*% solve(Cor01[c(2,4),c(2,4)]) %*% t(cbind(M1.dist,M0.dist))), sd=sqrt(Cor01[3,3]-Cor01[3,c(2,4)]%*%solve(Cor01[c(2,4),c(2,4)])%*%t(t(Cor01[3,c(2,4)])))))
    
    y.sample <- cbind(y11,y00,exp(s.m10), exp(s.m11))
    
    c <- 0.5 #<- threshold
    
    #----- Define Strata
    abs <- which(abs(y.sample[,3]-y.sample[,4]) <= c)
    neg <- which((y.sample[,3]-y.sample[,4]) > c)
    
    #----- Estimating Principal Causal Effects
    Pr1 <- length(abs)
    Pr2 <- length(neg)
    mean.y1 <- mean(y.sample[abs,1]); mean.y0 <- mean(y.sample[abs,2]);
    EDE <- mean.y1-mean.y0
    
    mean.y1 <- mean(y.sample[neg,1]); mean.y0 <- mean(y.sample[neg,2]);
    EAE <- mean.y1-mean.y0
    
    #----- Principal Causal Effects
    PS[tt,] <- c(EDE,EAE,Pr1,Pr2)
  }
}
save.image("MCMC_data.RData")
