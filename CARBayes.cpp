#include <RcppArmadillo.h>
using namespace Rcpp;


// This is a modified version of one of Rcpp functions in CARBayes R package (Duncan Lee, 2013):

// [[Rcpp::depends(RcppArmadillo)]]
  
  // [[Rcpp::export]]
  List poissoncarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                          NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                          double tau2, const NumericVector y, const double phi_tune, 
                          double rho, NumericVector offset, NumericMatrix h, NumericMatrix R,  NumericMatrix invR,  const int index)
  {
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi, proposal_var;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    double hcur1, hcur2, mu=0, sigma=1;
    NumericVector hi(4);
    NumericVector invRj(4);
    
      
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
      // Calculate prior variance
      priorvardenom = rho * (Wtripletsum[j] - 1) + 1 - rho;
      priorvar = tau2 / priorvardenom;
        
      // Calculate the prior mean
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      sumphi = 0;
        if(rowstart == rowend) {
            sumphi = 0;
        } else {
      for(int l = (rowstart+1); l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        }
      priormean = rho * sumphi / priorvardenom; 
      
      // propose a value
      proposal_var = priorvar * phi_tune;
      propphi = rnorm(1, phinew[j], pow(proposal_var, 0.5))[0];
      
      // Accept or reject it
      // Full conditional ratio
      hcur1 = R::qnorm(R::ppois(y[j], exp(phinew[j]+offset[j]), TRUE, FALSE), mu, sigma, TRUE, FALSE);
      hcur2 = R::qnorm(R::ppois(y[j], exp(propphi+offset[j]), TRUE, FALSE), mu, sigma, TRUE, FALSE);
        hi = h(j,_); //check this out hi=h(j,index)
        hi[index] = 0;
        invRj = invR(index,_);
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2);
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      lpold = offset[j] + phinew[j];
      lpnew = offset[j] + propphi;
        oldlikebit = y[j] * lpold - exp(lpold) + 0.5 * (1 - 1/R(index,index)) * pow(hcur1,2) - hcur1 * sum(invR(index,_) * hi);
      newlikebit = y[j] * lpnew - exp(lpnew) + 0.5 * (1 - 1/R(index,index)) * pow(hcur2,2) - hcur2 * sum(invR(index,_) * hi);
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      
      // Acceptace or reject the proposal
      if(runif(1)[0] <= acceptance) 
      {
        phinew[j] = propphi;
        accept = accept + 1;
      }
      else
      { 
      }    
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
  }
  

  // [[Rcpp::export]]
  double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                  NumericVector phi, NumericVector theta, double rho)
  {
    // Compute a quadratic form for the random effects
    // Create new objects 
    double tau2_posteriorscale;
    double tau2_quadform = 0, tau2_phisq = 0;
    int row, col;
    
    
    // Compute the off diagonal elements of the quadratic form
    for(int l = 0; l < n_triplet; l++)
    {
      row = Wtriplet(l,0) - 1;
      col = Wtriplet(l,1) - 1;
      if(row == col) {
        tau2_quadform = tau2_quadform + 0;
      } else {
        tau2_quadform = tau2_quadform + phi[(Wtriplet(l,0) - 1)] * theta[(Wtriplet(l,1) - 1)] * Wtriplet(l,2); 
      }
    }
    
    
    // Compute the diagonal elements of the quadratic form          
    for(int l = 0; l < nsites; l++)
    {
      tau2_phisq = tau2_phisq + phi[l] * theta[l] * (rho * (Wtripletsum[l] -1) + 1 - rho);    
    }
    
    
    // Compute the quadratic form
    tau2_posteriorscale = 0.5 * (tau2_phisq - rho * tau2_quadform);
    
    
    // Return the simulated value
    return tau2_posteriorscale;
  }
  

 
  
  // [[Rcpp::export]]
  List gaussiancarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                  NumericVector Wtripletsum, const int nsites, NumericVector phi, double tau2, const NumericVector y, const double phi_tune, 
                                  double rho, NumericVector offset, NumericMatrix h, NumericMatrix R,  NumericMatrix invR,  const int index,
                                  double sd)
  {
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0, rowstart=0, rowend=0;
    double acceptance, sumphi, proposal_var;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    double hcur1, hcur2, mu=0, sigma=1;
    NumericVector hi(4);
    NumericVector invRj(4);
    

    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
      // Calculate prior variance
      priorvardenom = rho * (Wtripletsum[j] - 1) + 1 - rho;
      priorvar = tau2 / priorvardenom;
      
      // Calculate the prior mean
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      sumphi = 0;
      if(rowstart == rowend) {
        sumphi = 0;
      } else {
        for(int l = (rowstart+1); l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
      }
      priormean = rho * sumphi / priorvardenom; 
      

      // propose a value
      proposal_var = priorvar * phi_tune;
      propphi = rnorm(1, phinew[j], pow(proposal_var, 0.5))[0];
      
      // Accept or reject it
      // Full conditional ratio
      hcur1 = R::qnorm(R::pnorm(y[j], (phinew[j]+offset[j]), sd, TRUE, FALSE), mu, sigma, TRUE, FALSE);
      hcur2 = R::qnorm(R::pnorm(y[j], (propphi+offset[j]), sd, TRUE, FALSE), mu, sigma, TRUE, FALSE);
      hi = h(j,_); //check this out hi=h(j,index)
      hi[index] = 0;
      invRj = invR(index,_);
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2);
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      lpold = offset[j] + phinew[j];
      lpnew = offset[j] + propphi;
      oldlikebit = - 0.5 * pow((y[j] - lpold), 2) / pow(sd, 2) + 0.5 * (1 - 1/R(index,index)) * pow(hcur1,2) - hcur1 * sum(invR(index,_) * hi);
      newlikebit = - 0.5 * pow((y[j] - lpnew), 2) / pow(sd, 2) + 0.5 * (1 - 1/R(index,index)) * pow(hcur2,2) - hcur2 * sum(invR(index,_) * hi);
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      
      
      // Acceptace or reject the proposal
      if(runif(1)[0] <= acceptance) 
      {
        phinew[j] = propphi;
        accept = accept + 1;
      }
      else
      { 
      }    
    }
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
  }
  
  
  
