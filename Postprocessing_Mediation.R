# Run this code on clusters

library(data.table)
library(mnormt)
library(gtools)
library(numDeriv)
library(matrixcalc)
library(corpcor)
library(rootSolve)

#------ Load MCMC samples and Data
load("MCMC_dataset.RData")
load("Dataset_2005_2012_noMedicare.RData")


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


#----- Extract MCMC samples from each marginal distribution and R
data.y0 <- para.y0[,1:(dim.cov+1)]
data.y1 <- para.y1[,1:(dim.cov+1)]
data.m0 <- para.m0[,1:(dim.cov+1)]
data.m1 <- para.m1[,1:(dim.cov+1)]

COR <- para.C

#----- Parameters of the models under Z=1
gamma11 <- data.y1[,(1):(dim.cov)]
beta1_11 <- data.m1[,(1):(dim.cov)]
s1_1 <- data.m1[,dim.cov+1]
w1_1 <- data.y1[,dim.cov+1]

#----- Parameters of the models under Z=0
gamma01 <- data.y0[,(1):(dim.cov)]
beta1_01 <- data.m0[,(1):(dim.cov)]
s1_0 <- data.m0[,dim.cov+1]
w1_0 <- data.y0[,dim.cov+1]


Thin <- 20  # Thining
Burn0 <- 10000  # Extra burn-in periods for Models under Z=0
Burn1 <- 10000   # Extra burn-in periods for Models under Z=1

#----- Set the number of post-processing steps, N
n.iter <- 1000

##---------------------------------------------------------------
## 'Main' function on each cluster (parallel)
##---------------------------------------------------------------

#----- Index of iterations (on parallel)
j <- process # process number from your parallel code (should be from 1:n.iter)

#----- The number of covariates samples (with replacement) for the empirical distribution
size <- n1+n0
s.covariate <- rbind(x1,x0)

#----- Index of posterior samples
index0 <- Thin * j + Burn0; index1 <- Thin * j + Burn1

#----- Construct the correlation matrix
Cor01 <- matrix(c(1,COR[[index0]]$prop1[1:3],
                  COR[[index0]]$prop1[1],1,COR[[index0]]$prop1[4:5],
                  COR[[index0]]$prop1[2],COR[[index0]]$prop1[4],1,
                  COR[[index0]]$prop1[6],COR[[index0]]$prop1[3],
                  COR[[index0]]$prop1[c(5,6)],1),4,4,byrow=TRUE)

#---- Sampling mediators conditional on covariates
s.m10 <- rnorm(n0+n1, mean = s.covariate%*%beta1_01[index0,], sd = sqrt(s1_0[index0]))
s.m10[(1+n1):(n0+n1)] <- m0
s.m11 <- rnorm(n0+n1, mean = s.covariate%*%beta1_11[index1,], sd = sqrt(s1_1[index1]))
s.m11[(1):(n1)] <- m1

C.sample <- t(cbind(s.m10,s.m11))

#----- Conditional distribution of Y{1;M(1)}
F.y1m1 <- function(o1){
  mm1 <- C.sample[o1, ]
  M1.dist <- qnorm(pmin(1-0.1^15,pmax(pnorm(mm1, mean = s.covariate%*%beta1_11[index1,], sd = sqrt(s1_1[index1])), 0.1^15)))
  con.joint <- pnorm(rnorm(n1+n0, Cor01[1,2] %*% solve(Cor01[2,2]) %*% t(M1.dist), sd=sqrt(Cor01[1,1]-Cor01[1,2]%*%solve(Cor01[2,2])%*%t(Cor01[1,2]))))
  result <- qnorm(pmin(1-0.1^15, pmax(0.1^15, con.joint )),mean = s.covariate%*%gamma11[index1,], sd = sqrt(w1_1[index1]) )
  return(result)
}

#----- Conditional distribution of Y_{0;M(0)}
F.y0m0 <- function(o1){
    mm1 <- C.sample[o1, ]
    M1.dist <- qnorm(pmin(1-0.1^15,pmax(pnorm(mm1, mean = s.covariate%*%beta1_01[index1,], sd = sqrt(s1_0[index1])), 0.1^15)))
    con.joint <- pnorm(rnorm(n1+n0,Cor01[3,4] %*% solve(Cor01[4,4]) %*% t(M1.dist), sd=sqrt(Cor01[3,3]-Cor01[3,4]%*%solve(Cor01[4,4])%*%t(Cor01[3,4]))))
    result <- qnorm(pmin(1-0.1^15, pmax(0.1^15, con.joint )),mean = s.covariate%*%gamma01[index1,], sd = sqrt(w1_0[index1]) )
    return(result)
}


#----- Values of potential outcomes
y10 <- F.y1m1(1)
y11 <- F.y1m1(2)
y00 <- F.y0m0(1)

#----- Causal Effects
TE <- y11 - y00
NIE <- y11 - y10
NDE <- y10 - y00

result<-cbind(TE, NIE, NDE)



