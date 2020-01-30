library(Rcpp); library(tidyverse); library(MCMCpack); 
library(RcppArmadillo);library(mgcv);library(MCMCpack)

BNPtesting <- function(y, 
                       prior_par,
                       NSIM = 20000, 
                       burn_in = 100, 
                       thinning = 1,
                       verbose=T, 
                       verbose_step=10, 
                       SM=.5, 
                       optThresh=.44, 
                       batch=50,
                       sed=12345){

  #####################################
  # Extraction of the hyperparam. from the input list
  ###################################### 
  # hyperparameters for P^*_1 = NIGmix(m1,V1,a1,b1)
  # m1 is estimated by data, the value here works as a starting value
  m1 <- prior_par$m1
  V1 <- prior_par$V1
  a1 <- prior_par$a1
  b1 <- prior_par$b1
  ###################################### 
  # hyperparameters for P^*_0 = N(m0,s0) * IG(a0,b0)
  m0 <- prior_par$m0
  s0 <- prior_par$s0
  a_0 <- prior_par$a_0
  b_0 <- prior_par$b_0
  
  ###################################### 
  # hyperparameters for \rho \sim Beta(a,b)
  a_rho  <- prior_par$a_rho
  b_rho  <- prior_par$b_rho
  
  ###################################### 
  # hyperparameters for two PY processes
  sigma0   <- prior_par$sigma0
  theta0   <- prior_par$theta0
  sigma1   <- prior_par$sigma1
  theta1   <- prior_par$theta1
  
  ###################################### 
  # hyperparameters Abs-NLP
  kappaNLP <- prior_par$kappaNLP
  s1       <- prior_par$s1
  
  ###########################################################
  # Containers
  ###########################################################
  
  n  <- length(y)
  zi <- numeric(n)
  ThetaCon <- array(0, dim=c(n,2,NSIM))
  LabCon   <- array(0, dim=c(n,2,NSIM))
  ZiCon    <- matrix(0, NSIM, n)
  RhoCon   <- TauCon <- m1Con <- numeric(NSIM)
  
  ###########################################################
  # Initialization of the algorithm
  ###########################################################
  # 
  # Obtain a starting configuration for the clustering and the parameters
  sy                               <- sort(y,ind=T)
  km                               <- kmeans(abs(sy$x),2)
  kmc                              <- km$cluster
  kmc[(km$cluster==km$cluster[1])] <- 0
  kmc[(km$cluster!=km$cluster[1])] <- 1
  zi[sy$ix]                        <- kmc 
  
  mean_k  <- variance_k    <- sample( x = quantile(y)[2:4], size = n, replace = T )
  mean_k[zi==1]  <- mean_k[zi==1]+runif(1)
  umk     <- unique(mean_k)
  
  for(q in 1:length(umk)){
    variance_k[which(mean_k==umk[q])] <- rgamma(1,a_0,b_0)
  }
  rho     <- rbeta(1,a_rho,b_rho)
  Data    <- cbind(yy = y, mu = mean_k, sig = variance_k, zi = zi, lab0=0, lab1=0)
  
  # Vector for Adaptive MH
  acc_m1 <- numeric(NSIM*thinning + burn_in)
  
  
  ###########################################################
  # GIBBS loop
  ###########################################################
  for(sim in 2:(NSIM*thinning + burn_in)){
    
    
    # Find unique combinations of mu,sigma and zi in the dataset and code them:
    Unique_dataset <- uniquecombs(Data[,2:4])
    
    
    ########### UPDATE Zi and Rho ##########################
    
    tempLIST <-   pk_ratio(yyy = Data[,1], 
                           TH  = Data[,2], 
                           zi  = Data[,4],
                           
                           a_rho = a_rho, 
                           b_rho = b_rho,
                           
                           sigma0 = sigma0, sigma1 = sigma1, 
                           theta0 = theta0, theta1 = theta1,
                           
                           a0 = a_0, b0 = b_0,
                           a1 = a1, b1 = b1, 
                           m1 = m1, V1 = V1, 
                           m0 = m0, s0 = s0,
                           
                           n   = n,
                           K   = nrow(Unique_dataset),
                           uTH = Unique_dataset[,1],
                           zz  = Unique_dataset[,3])
    
    
    rho       <- tempLIST$rho
    Data[,4]  <- tempLIST$zi
    
    ########### Update Theta  ###############################
    
    Data[,2:3] <- Get_THETA_updated(data = Data, 
                                    sig0 = sigma0, 
                                    theta0 = theta0, 
                                    sig1 = sigma1, 
                                    theta1 = theta1,
                                    a0=a_0, b0 = b_0,
                                    m0 = m0,s0 = s0,
                                    m1 = m1,V1 = V1,
                                    a1 = a1,b1 = b1)
    
    ##################################
    # Once theta is updated, we update the correspondig label.
    # z partitions the observations into Null (0) and Alternative Process (1)
    # lab0 > 0 gives the cluster for the observations | z=0
    # lab1 > 0 gives the cluster for the observations | z=1
    # Needed for the acceleration step
    
    if(sum(Data[,4]==0)>1){
      UC0 <- uniquecombs(Data[Data[,4]==0,2:3])
      lab0 <- numeric(n); lab0[Data[,4]==0] <- attributes(UC0)$index
      Data[,5] <- lab0
    }
    if(sum(Data[,4]==1)>1){
      UC1 <- uniquecombs(Data[Data[,4]==1,2:3])
      lab1 <- numeric(n); lab1[Data[,4]==1] <- attributes(UC1)$index
      Data[,6] <- lab1
    }
    
    
    ########### Acceleration Step  ###########################
    
    Data[,2:3] <- Accel_step(Data = Data, 
                             m0 = m0,
                             s0 = s0, 
                             m1 = m1, V1 = V1,
                             a1 = a1, b1 = b1,
                             a0 = a_0, b0 =b_0)
    
    ########## Update m1 #####################################
    if(sum(Data[,4]==1)>0){
      UNI_THETA_1 <- uniquecombs(Data[Data[,4]==1,2:3])
      mUNI_THETA_1 <- as.matrix(UNI_THETA_1)
      if(ncol(mUNI_THETA_1)==1){
        mUNI_THETA_1 = t(mUNI_THETA_1)
      }
      emme1 <- m_MH_mix(previous_m = m1,
                                     V =  V1, 
                                     s1 = s1, 
                                     unique_mu1  = mUNI_THETA_1[,1], 
                                     unique_sig1 = mUNI_THETA_1[,2], 
                                     kappa = kappaNLP,
                                     sigma_m = SM)
    }else{
      emme1 <- m_MH_mix_prior(previous_m = m1,
                                           s1 = s1, 
                                           kappa = kappaNLP,
                                           sigma_m = SM)
    }
    
    m1          <- emme1[1]
    acc_m1[sim] <- emme1[2]

    if( ((sim-1) %% batch) == 0){
      SM = exp( ifelse(mean(acc_m1[(sim-batch):sim]) < optThresh, 
                       (log(SM) - min(0.01,1/sqrt(sim) )) , 
                       (log(SM) + min(0.01,1/sqrt(sim)) ) ))
    }
    
    ############################################################################
    
    if (sim > burn_in & ((sim - burn_in) %% thinning == 0)) {
      rr                      <- floor((sim - burn_in)/thinning)
      ThetaCon[,,rr]          <- Data[,2:3]
      ZiCon[rr,]              <- Data[,4]
      RhoCon[rr]              <- rho
      LabCon[,,rr]            <- Data[,5:6]
      m1Con[rr]               <- m1
    }
    # Printing the results
    if (verbose) {
      if (sim%%(verbose_step*thinning) == 0) 
        cat(paste("Sampling iteration: ", sim, " out of ",NSIM*thinning + burn_in, "\n",
                  sep = ""))}
  }
  
  out <- list(ThetaCon  = ThetaCon,
              ZiCon     = ZiCon ,
              RhoCon    = RhoCon,
              LabCon    = LabCon,
              m1Con     = m1Con,
              Y         = y, 
              prior_par = prior_par
              )
  class(out) <- "BNPtesting" #Code-name for our model
  return(out)
}

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
