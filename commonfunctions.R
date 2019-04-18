#### expit() function
expit <- function(x){exp(x)/(1+exp(x))}


#### log-likelihood (Normal)
loglik <- function(Y, X, beta0, beta, sigma, ran){
  lik<-sum(dnorm(Y, ran+X%*%beta, sqrt(sigma), log=TRUE) )
  return(lik)
}

#### log-likelihood (Poisson)
loglik_log <- function(Y, X, beta0, beta, off, ran){
  lik<-sum(dpois(Y, exp(ran+X%*%beta)*off,log=TRUE) )
  return(lik)
}


#### log priors of hyper-parameters
logprior <- function(Y, X, beta0, beta, sigma, 
                     mu_beta0, sigma_beta0, alpha_beta0, alpha_sigma, mu_pop){
  
  beta_part <- log(dmnorm((beta), rep(0,dim(X)[2]), diag(100,dim(X)[2]))) # Coefficients (other than the intercept) parameters

  return(beta_part)
}

#### log-Copula
logCopula=function(h,index,h.cur,R,zz=1){
  
  h[,index] <- h.cur;
  hR <- h%*%chol(solve(R))
  result <- 0.5*sum(rowSums(h*h)-rowSums(hR*hR))
  return(result)
}



#### Check the W matrix
common.Wcheckformat <- function(W)
{
  #### Check W is a matrix of the correct dimension
  if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
  n <- nrow(W)
  if(ncol(W)!= n) stop("W is not a square matrix.", call.=FALSE)    
  
  #### Check validity of inputed W matrix
  if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
  if(min(W)<0) stop("W has negative elements.", call.=FALSE)
  if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
  if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
  
#### Create the triplet form
W.triplet <- c(NA, NA, NA)
for(i in 1:n)
  {
    for(j in 1:n)
    {
      if(W[i,j]>0)
      {
        W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
      }else{}
    }
  }
  W.triplet <- W.triplet[-1, ]     
  n.triplet <- nrow(W.triplet) 
  W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
  n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
  
  #### Create the start and finish points for W updating
  W.begfin <- array(NA, c(n, 2))     
  temp <- 1
  for(i in 1:n)
  {
    W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
    temp <- temp + n.neighbours[i]
  }
  
  #### Return the critical quantities
  results <- list(W=W, W.triplet=W.triplet, n.triplet=n.triplet, W.triplet.sum=W.triplet.sum, n.neighbours=n.neighbours, W.begfin=W.begfin, n=n)
  return(results)   
}


##### metropolis-hastings for the correlation matrix
metropolisC=function(h,rho){
  
  library(mnormt)
  acc1<-rep(7)
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

# Adjust acceptance rates
common.accceptrates1 <- function(accept, sd, min, max)
{
  #### Update the proposal standard deviations
  rate <- 100 * accept[1] / accept[2]
  
  if(rate > max)
  {
    sd <- sd + 0.1 * sd
  }else if(rate < min)              
  {
    sd <- sd - 0.1 * sd
  }else
  {
  }
  return(sd)
}


