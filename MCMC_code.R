#------ Load the main dataset
load("Dataset_2005_2012_noMedicare.RData") # Note that the dataset doesn't contain the Medicare data (i.e. the outcome)
library(rootSolve)
library(matrixcalc)
library(data.table)
library(mvtnorm)
library(gtools)
library(numDeriv)
library(corpcor)

#------ Define 3 EPA regions by state names
Northeast = c("ME", "NH", "VT", "NY", "PA", "DE", "NJ", "MD", "DC", "VA", "MA", "CT", "RI")
IndustrialMidwest = c("WV", "OH", "KY", "IN", "IL", "WI", "MI")
Southeast = c("FL", "GA", "SC", "NC", "TN", "AL", "MS", "AR","LA")

#------ Restrict the dataset
Master <- subset(Master, State.zip %in% c(IndustrialMidwest, Northeast, Southeast))
Master <- subset(Master, !(State.zip %in% c("ME", "FL", "AR", "LA", "MS","VT","NH","WI")))
Master <- subset(Master, !is.na(as.numeric(Master$TerrainElevation)))
Master <- subset(Master, !is.na(PM2005))
Master <- subset(Master, !is.na(PctHighSchool))
Master <- subset(Master, TotPop > 100)
Master <- subset(Master, Tot_den_for_death_MA_FFS2005 > 100)
Master <- subset(Master, Person_year_FFS2005 > 100)

Master$TEMP2005 <- Master$TEMP2005 - 273.15 # kelvin to celsius

#------ Dichotomize the exposure level
Master$X <- ifelse(Master$Exposure2005 < summary(Master$Exposure2005)[4], 1,0 )  # 


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

#------ New dataset for the treated (1) vs. control (0)
Data<-pair
m0 <- as.numeric(subset(Data, X==0)$PM2005)
m1 <- as.numeric(subset(Data, X==1)$PM2005)
y0 <- (1000*subset(Data, X==0)$AMI2005/subset(Data, X==0)$Person_year_FFS2005) # outcome from the Medicare data
y1 <- (1000*subset(Data, X==1)$AMI2005/subset(Data, X==1)$Person_year_FFS2005) # outcome from the Medicare data


x0 <- as.matrix(subset(Data2, X==0)[,list(1,TEMP2005,
                                          smokerate2005,
                                          mean_age2005,
                                          log(MedianHHInc+1),
                                          Female_rate2005,
                                          Black_rate2005,
                                          PctHighSchool,
                                          PctUrban,
                                          PctHisp,
                                          log(TotPop+1)),ps])
x1 <- as.matrix(subset(Data2, X==1)[,list(1,TEMP2005,
                                          smokerate2005,
                                          mean_age2005,
                                          log(MedianHHInc+1),
                                          Female_rate2005,
                                          Black_rate2005,
                                          PctHighSchool,
                                          PctUrban,
                                          PctHisp,
                                          log(TotPop+1)),ps])


n0 <- dim(x0)[1] # Num. observations in the control
n1 <- dim(x1)[1] # Num. observations in the treated

dim.cov<-dim(x0)[2] # Num. of Covariates


#------- define parameters
P <- 4 # number of marginal distributions



#------- Metropolis-hastings for the correlation matrix R
metropolisC=function(h,rho){
  
  prop1<-rho
  interval1 <- NULL
  interval2 <- NULL

  seq <- 1:6
  
  for(q in seq){
    prop1[q] <- 1
    RHO1 <- det(matrix(c(1,prop1[1:3],prop1[1],1,prop1[4:5],prop1[2],prop1[4],1,prop1[6],prop1[3],prop1[c(5,6)],1),4,4,byrow=TRUE))
    prop1[q] <- 0
    RHO0 <- det(matrix(c(1,prop1[1:3],prop1[1],1,prop1[4:5],prop1[2],prop1[4],1,prop1[6],prop1[3],prop1[c(5,6)],1),4,4,byrow=TRUE))
    prop1[q] <- -1
    RHO_1 <- det(matrix(c(1,prop1[1:3],prop1[1],1,prop1[4:5],prop1[2],prop1[4],1,prop1[6],prop1[3],prop1[c(5,6)],1),4,4,byrow=TRUE))
    Fun <- function(r){(RHO1+RHO_1-2*RHO0)/2*r^2 + (RHO1-RHO_1)/2*r+RHO0}
    
    interval1 <- min(multiroot(Fun, start=c(-1,1))$root)
    interval2 <- max(multiroot(Fun, start=c(-1,1))$root)
    prop1[q] <- runif(1, rho[q]-0.1,rho[q]+0.1)
    
    COR=matrix(c(1,prop1[1:3],prop1[1],1,prop1[4:5],prop1[2],prop1[4],1,prop1[6],prop1[3],prop1[c(5,6)],1),4,4,byrow=TRUE)
    RHO=matrix(c(1,rho[1:3],rho[1],1,rho[4:5],rho[2],rho[4],1,rho[6],rho[3],rho[c(5,6)],1),4,4,byrow=TRUE)
    
    if(is.positive.definite(COR)){
      rat= sum(dmnorm(h, mean=rep(0,4), varcov=COR, log=TRUE))+
        sum(dunif(prop1[q], interval1, interval2,log=TRUE))+
        sum(dunif(rho[q], prop1[q]-0.1,prop1[q]+0.1,log=TRUE))-
        sum(dmnorm(h, mean=rep(0,4), varcov=RHO, log=TRUE))-
        sum(dunif(rho[q], interval1, interval2,log=TRUE))-
        sum(dunif(prop1[q], rho[q]-0.1,rho[q]+0.1,log=TRUE))
      if(is.na(rat)){
        prop1[q]=rho[q]
      }else{
        if (log(runif(1))>rat) {
          prop1[q]=rho[q]
        }else{rho[q]=prop1[q]}}}else{
          prop1[q]=rho[q]
        }
  }
  return(list(prop1=prop1, acc1=acc1))
}


#------- log likelihood
loglik <- function(Y, X, beta, sigma){
  lik<-sum(dnorm(Y, X%*%beta, sqrt(sigma), log=TRUE) )
  return(lik)
}

#------- log priors of the parameters
logprior <- function(Y, X,  beta){
  beta_part <- log(dmnorm((beta), rep(0,dim(X)[2]), diag(100,dim(X)[2]))) 
  return(beta_part)
}

#------ log Copula evaluation
logCopula=function(h,index,h.cur,R,zz=1){
  if(index %in% c(1,2)){
    h[1:(n1),index] <- h.cur; 
    hR <- h[1:n1,]%*%chol(solve(R))
    result <- 0.5*sum(rowSums(h[1:n1,]*h[1:n1,])-rowSums(hR*hR))                 
  }else{
    h[(1+n1):(n0+n1),index] <- h.cur;
    hR <- h[(1+n1):(n0+n1),]%*%chol(solve(R))
    result <- 0.5*sum(rowSums(h[(1+n1):(n0+n1),]*h[(1+n1):(n0+n1),])-rowSums(hR*hR))                 
  }
  return(result)
}


#------- Metropolis-hastings for the marginals
metropolis <- function(h, X, Y,R, beta_pre, cov2, zz, index){
  # if zz=1, this is for the data under Z=1; if zz=0, this is for the data under Z=0
  if(zz==1){
    YY <- Y[(1:n1)]; XX <- X[(1:n1),];
  }else{
    YY <- Y[(n1+1):(n1+n0)]; XX <- X[(n1+1):(n1+n0),];
  }
  
  
  prop3=beta_pre
  pprop3=beta_pre
  Cov2 = 2.38^2/((dim.cov)*2)*(cov2/2+diag(0.1^10,dim.cov))

  pprop3=(rmnorm(1,mean=beta_pre, Cov2)) # Proposal distribution for the coefficients
  
  h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, XX%*%(prop3), sqrt(sigma_pre)))))
  h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, XX%*%(pprop3), sqrt(sigma_pre)))))
  
  rat <- logCopula(h,index,h.cur2,R,zz) + loglik(YY,XX,pprop3,sigma_pre) +
    logprior(YY,XX,pprop3) +
    dmnorm(prop3, (pprop3), Cov2,log=TRUE) -
    logCopula(h,index,h.cur1,R,zz) - loglik(YY,XX,prop3,sigma_pre) -
    logprior(YY,XX,prop3) -
    dmnorm((pprop3), prop3, Cov2,log=TRUE)
  if(is.na(rat)){
    pprop3=beta_pre
  }else{
    if (log(runif(1))>rat) {
      pprop3=beta_pre
    }else{prop3=pprop3}
  }

  
  prop4 <- pprop4 <- sigma_pre
  log_prop4 <- log_pprop4 <- log(sigma_pre)
  
  log_pprop4[W] <- rnorm(1, log(sigma_pre), 0.1) # Proposal distribution for the variance
    
  h.cur1 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, XX%*%(prop3), sqrt(exp(log_prop4))))))
  h.cur2 <- qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(YY, XX%*%(prop3), sqrt(exp(log_pprop4))))))

  rat <- logCopula(h,index,h.cur2,R,zz) + loglik(YY,XX,prop3,exp(log_pprop4)) +
      dgamma(exp(log_pprop4), 5, scale=0.2, log=TRUE)+log(abs(exp(log_pprop4)))+
      dnorm(log_prop4, log_pprop4,0.1,log=TRUE)-
      logCopula(h,index,h.cur1,R,zz) - loglik(YY,XX,prop3,exp(log_prop4)) -
      dgamma(exp(log_prop4), 5, scale=0.2, log=TRUE)-log(abs(exp(log_prop4)))-dnorm(log_pprop4, log_prop4,0.1,log=TRUE)
    
  if(log(runif(1))>rat){
    pprop4[W] <- sigma_pre[W]
  }else{
    prop4[W] <- exp(log_pprop4[W])
  }
  
  return(c(prop3,prop4))
}

#-----------------------
#------ Run MCMC -------
#-----------------------
MCMC=30000 # Num. of Iterations

fit.y0 <- lm(y0~as.matrix(x0)-1)
fit.y1 <- lm(y1~as.matrix(x1)-1)
fit.m0 <- lm(m0~as.matrix(x0)-1)
fit.m1 <- lm(m1~as.matrix(x1)-1)

#------- Initial Values
para.y1=matrix(nrow=MCMC,ncol=dim.cov+1)
para.y1[1,] <- para.y1[2,] <- c(coef(fit.y1),1)

para.y0=matrix(nrow=MCMC,ncol=dim.cov+1)
para.y0[1,] <- para.y0[2,] <- c(coef(fit.y0),1)

para.m1 <- matrix(nrow = MCMC, ncol = dim.cov+1) 
para.m1[1,] <- para.m1[2,] <- c(coef(fit.m1),1)

para.m0 <- matrix(nrow = MCMC, ncol = dim.cov+1) 
para.m0[1,] <- para.m0[2,] <- c(coef(fit.m0),1)

para.C=list(prop1=c(cor(y1,m1),0,0,0,0,cor(y0,m0)), rrho=rep(0,4))
para.C=list(para.C, para.C)
R=matrix(c(1,cor(y1,m1),0,0,cor(y1,m1),1,0,0,0,0,1,cor(y0,m0),0,0,cor(y0,m0),1),byrow=TRUE, nrow=4)
h=rmnorm(n0+n1, mean=rep(0,4),R)

cov2.m1 <- vcov(lm(m1~x1-1))
cov2.m0 <- vcov(lm(m0~x0-1))

cov2.y1 <- vcov(lm(y1~x1-1))
cov2.y0 <- vcov(lm(y0~x0-1))

y1.all <- c(y1,rep(mean(y1),n0))
y0.all <- c(rep(mean(y0),n1), y0)

m1.all <- c(m1,rep(mean(m1),n0))
m0.all <- c(rep(mean(m0),n1), m0)


pb <- txtProgressBar(min = 0, max = MCMC, style = 3)
for (t in 3:MCMC){
  SEQ <- seq(54,MCMC, by=25)

  # sampling Y(0)
  if(t > 200){
    cov2.y0 <- cov(para.y0[3:(t-1),(1):(dim.cov)])
  }

  para.y0[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=y0.all, R=R,
                            beta_pre=para.y0[t-1,(1):(dim.cov)], sigma_pre=para.y0[t-1,(dim.cov+1)],index=3, cov2=cov2.y0, zz=0)
  
  y0.all[1:(n1)] <- rnorm(n1,x1%*%para.y0[t,(1):(dim.cov)], sqrt(para.y0[t-1,(dim.cov+1)]))
  
  h[,3] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(y0.all, rbind(x1,x0)%*%para.y0[t,(1):(dim.cov)], sqrt(para.y0[t-1,(dim.cov+1)])))))
  
  
  # sampling Y(1)
  if(t > 200){
    cov2.y1 <- cov(para.y1[3:(t-1),(1):(dim.cov)])
  }
  
  para.y1[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=y1.all, R=R,
                            beta_pre=para.y1[t-1,(1):(dim.cov)], sigma_pre=para.y1[t-1,(dim.cov+1)],index=1, cov2=cov2.y1, zz=1)
  
  y1.all[(n1+1):(n0+n1)] <- rnorm(n0,x0%*%para.y1[t,(1):(dim.cov)], sqrt(para.y1[t-1,(dim.cov+1)]))
  
  h[,1] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(y1.all, rbind(x1,x0)%*%para.y1[t,(1):(1+dim.cov)], sqrt(para.y1[t-1,(1+dim.cov+1)])))))
  
  # sampling M(0)
  if(t > 200){
    cov2.m0 <- cov(para.m0[3:(t-1),(1):(dim.cov)])
  }
  
  para.m0[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m0.all, R=R,
                            beta_pre=para.m0[t-1,(1):(dim.cov)], sigma_pre=para.m0[t-1,(dim.cov+1)], index=4,cov2=cov2.m0, zz=0)
  
  m0.all[1:(n1)] <- rnorm(n1,x1%*%para.m0[t,(1):(dim.cov)], sqrt(para.m0[t-1,(dim.cov+1)]))
  
  h[,4] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(m0.all, rbind(x1,x0)%*%para.m0[t,(1):(1+dim.cov)], sqrt(para.m0[t-1,(dim.cov+1)])))))
  

  # sampling M(1)
  if(t > 200){
    cov2.m1 <- cov(para.m1[3:(t-1),(1):(dim.cov)])
  }
  
  para.m1[t,] <- metropolis(h=h, X=rbind(x1,x0), Y=m1.all, R=R, 
                            beta_pre=para.m1[t-1,(1):(dim.cov)], sigma_pre=para.m1[t-1,(dim.cov+1)], index=2, cov2=cov2.m1, zz=1)
  
  m1.all[(n1+1):(n0+n1)] <- rnorm(n0,x0%*%para.m1[t,(1):(dim.cov)], sqrt(para.m1[t-1,(dim.cov+1)]))
  
  h[,2] <-  qnorm(pmin(1-0.1^15,pmax(0.1^15,pnorm(m1.all, rbind(x1,x0)%*%para.m1[t,(1):(dim.cov)], sqrt(para.m1[t-1,(dim.cov+1)])))))

  
  # sampling R
  para.C[[t]]=metropolisC(h=h,rho=para.C[[t-1]]$prop1)
  Prop1 <- para.C[[t]]$prop1
  R = matrix(c(1,Prop1[1:3],Prop1[1],1,Prop1[4:5],Prop1[2],Prop1[4],1,Prop1[6],Prop1[3],Prop1[c(5,6)],1),4,4,byrow=TRUE)
  
  Sys.sleep(0.001)
  setTxtProgressBar(pb, t)
}

save.image("MCMC_dataset.RData")
