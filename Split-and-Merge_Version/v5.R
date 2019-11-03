library(Rcpp); library(tidyverse); library(MCMCpack); 
library(RcppArmadillo);library(mgcv);
#Rcpp::sourceCpp(here::here("code/ULTIMAVERSIONE/v3.cpp"))


GIBBS_SAMPLER_BNPtesting_v5<- function(NSIM = 20000, 
                                        burn_in = 100, 
                                        thinning = 1,
                                        y, 
                                        prior_par, 
                                        verbose=T, 
                                        verbose_step=10, 
                                        SM=.5, 
                                        optThresh=.44, 
                                        batch=50,
                                        random=T,
                                        sed=12345,
                                        StartingSWEEP=50,
                                        SWEEPevery=50,warmstart=F,
                                        keepTheta=F,densities=F){
  
  ###########################################################
  # Useful quantities
  ###########################################################
  time1 <- Sys.time()
  n     <- length(y)
  zi    <- numeric(n)
  ###################################### hyperparameters for P^*_1 = NIGmix(m1,V1,a1,b1)
  # m1 is estimated by data
  m1 <- prior_par$m1
  V1 <- prior_par$V1
  a1 <- prior_par$a1
  b1 <- prior_par$b1
  
  
  ###################################### hyperparameters for P^*_0 = N(m0,s0) * IG(a0,b0)
  m0  <- prior_par$m0
  s0  <- prior_par$s0
  a_0 <- prior_par$a_0
  b_0 <- prior_par$b_0
  
  ###################################### hyperparameters for \rho \sim Beta(a,b)
  a_rho  <- prior_par$a_rho
  b_rho  <- prior_par$b_rho
  
  ###################################### hyperparameters for two PY processes
  sigma0   <- prior_par$sigma0
  theta0   <- prior_par$theta0
  sigma1   <- prior_par$sigma1
  theta1   <- prior_par$theta1
  
  ###################################### hyperparameters Abs-NLP
  kappaNLP <- prior_par$kappaNLP
  s1       <- prior_par$s1
  
  ###########################################################
  # Containers
  ###########################################################
  
  ThetaCon <- array(0, dim=c(n,2,NSIM))
  LabCon   <- array(0, dim=c(n,2,NSIM))
  ZiCon    <- matrix(0, NSIM, n)
  RhoCon   <- m1Con <- numeric(NSIM)
  
  ##############################################################
  ygrid = seq(min(y)-2,max(y)+2,length.out=1000)
  d0 = d1 = d01 = numeric(1000)
  postpred <- cbind(d0,d1,d01)
  ##############################################################
  
  
  ###########################################################
  # Initialization of the algorithm
  ###########################################################
  if(!random) {set.seed(sed)}
  
  
  ######################################3
  
  
  #####################################
  if(warmstart){
  sy  <- sort(y,ind=T)
  km  <- kmeans(abs(sy$x),2)
  kmc <- km$cluster
  kmc[(km$cluster==which.min(km$centers))] <- 0
  kmc[(km$cluster!=which.min(km$centers))] <- 1
  zi[sy$ix] <- kmc
  }else{
   zi <- sample(0:1,n,replace = T) 
  }
  mean_k  <- variance_k    <- sample( x = quantile(y)[2:4], size = n, replace = T )
  mean_k[zi==1]  <- mean_k[zi==1]+runif(1)
  umk     <- unique(mean_k)
  
  for(q in 1:length(umk)){
    variance_k[which(mean_k==umk[q])] <- rgamma(1,a_0,b_0)
  }
  
  rho         <- rbeta(1,a_rho,b_rho)
  Data        <- cbind(yy = y, mu = mean_k, sig = variance_k, zi = zi, lab0=0, lab1=0)
  
  acc_m1 <- numeric(NSIM*thinning + burn_in)
   
  #cat(paste("Starting the Gibbs loop"))
  ###########################################################
  # GIBBS loop
  ###########################################################
  for(sim in 2:(NSIM*thinning + burn_in)){
    
    
    # Find unique combinations of mu,sigma and zi in the dataset:
    Unique_dataset <- uniquecombs(Data[,2:4])
    
    
    ########### UPDATE Zi and Rho ##########################
    
    tempLIST <-   pk_ratio(yyy = Data[,1],                          # indipendente dalle etichette di gruppo, dati i theta precedenti decide se assegnare obs a processo 0 o 1
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
    Data[,4]  <- tempLIST$zi                                              #### Attento, ora certe etichette di gruppo possono essere meaningless se una obs viene riassegnata a processo alternativo
    
    ########### Update Theta  ###############################
    if(sim < StartingSWEEP || sim %% SWEEPevery == 0){
     
     Data[,2:3] <- Get_THETA_updated(data = Data,                          ##### totalmente indipendente da etichette di gruppo, Data[,5:6]. Rialloca le osservazioni DATO il processo e valori unici di theta.
                                      sig0 = sigma0,                       # |
                                      theta0 = theta0,                     # |
                                      sig1 = sigma1,                       # |
                                      theta1 = theta1,                     # |
                                      a0 = a_0, b0 = b_0,                  # |
                                      m0 = m0,  s0 = s0,                   # |
                                      m1 = m1,  V1 = V1,                   # |
                                      a1 = a1,  b1 = b1)                   # |
     ##################################                                    # v
     Data <- RESET_LAB0LAB1_GivenTHETA(Data)                               # Qui riassegno le labels dati i valori unici di theta
     ########################################################## 
     # Data[,2:3] <- Accel_step(Data = Data,
     #                          m0 = m0,  s0 = s0,
     #                          a0 = a_0, b0 =b_0,
     #                          m1 = m1, V1 = V1,
     #                          a1 = a1, b1 = b1)
     
     }else{
     #########################################################
     # lo split and merge invece dipende essenzialmente dalle etichette 
     # quindi sistemiamo prima dati i nuovi valori del processo
     Data <- RESET_LAB0LAB1_GivenTHETA(Data)  
     Data[,2:3] <- Accel_step(Data = Data,
                              m0 = m0,  s0 = s0,
                              a0 = a_0, b0 =b_0,
                              m1 = m1, V1 = V1,
                              a1 = a1, b1 = b1)
     #Data <- RESET_LAB0LAB1_GivenTHETA(Data)

     Data <- SPLIT_MERGE_PROCESS0(Data = Data,
                                  sigma0 = sigma0,theta0 = theta0,
                                  a_0 = a_0,b_0 = b_0)
     Data <- SPLIT_MERGE_PROCESS1(Data = Data,
                                  sigma1 = sigma1,theta1 = theta1,
                                  a1 = a1,b1 = b1,m1 = m1,V1 = V1)
    
     }
    
    
    
    ############ Acceleration Step  ###########################")
    
    Data[,2:3] <- Accel_step(Data = Data,
                               m0 = m0,  s0 = s0,
                               a0 = a_0, b0 =b_0,
                               m1 = m1, V1 = V1,
                               a1 = a1, b1 = b1)
      
      ##############################################################
      if(densities & sim > burn_in){
      postpred <- posterior_predictive2(Data = Data,postpred = postpred,NSIM = NSIM,ygrid = ygrid,rho = rho,
                                        t0 = theta0,t1 = theta1,s0 = sigma0,s1 = sigma1,a0 = a_0,b0 = b_0,
                                        a1 =  a1,b1 = b1, m1 = m1,V1 = V1 )
      }
      
      
      
      
      
      
  # cat("########## Update m1 #####################################")
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
    
    
    m1     <- emme1[1]
    acc_m1[sim] <- emme1[2]
    
    
    if( ((sim-1) %% batch) == 0){
    #  cat("previous SM:", SM,"...",min(0.01,1/sqrt(sim) ),"\n")
      SM = exp( ifelse( mean(acc_m1[(sim-batch):sim]) < optThresh,
                       (log(SM) - min(0.01,1/sqrt(sim) )) ,
                       (log(SM) + min(0.01,1/sqrt(sim)) ) ))
     # cat("ACCRATE:",mean(acc_m1[(sim-batch):sim]),"NEW SM:",SM,"\n")
      }

    
#    if((sim%%500==0) & sim > burn_in) plot(apply(ZiCon[1:sim,],2,mean)~y) 
    
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
  time2 <- Sys.time()
  dx <- diff(ygrid)[1]
  postpred <- apply(postpred,2,function(x) x/(sum(x)*dx)  )
  
  if(keepTheta){
  out   <- list(ThetaCon = ThetaCon,
                mppi    = apply(ZiCon,2,mean)   ,
                RhoCon   = RhoCon,
                m1Con = m1Con,
                m1_accrate = mean(acc_m1),
                #LabCon   = LabCon,
                Y=y, prior_par=prior_par,
                Elapsed_time = time2-time1,
                PostPred=cbind(ygrid,postpred))
    }else{
  out   <- list(mppi    = apply(ZiCon,2,mean)   ,
                RhoCon   = RhoCon,
                m1Con = m1Con,
                m1_accrate = mean(acc_m1),
                Y=y, prior_par=prior_par,
                Elapsed_time = time2-time1,
                PostPred=cbind(ygrid,postpred))
    
  }
  
  class(out) <- "BNPtesting" #Code-name for our model
  return(out)
}

###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
############## BRUTTA SPLIT MERGE






# # coviene distinguere fra indice di posizione (pos) e indice di gruppo (lab)
# all_lab0old <- all_lab0 <- Data[,5]
# all_lab1old <- all_lab1 <- Data[,6]
# all_pos0old <- all_pos0 <- which(Data[,5]>0)
# all_pos1old <- all_pos1 <- which(Data[,6]>0)
# 
# yy0      <- Data[all_pos0,1];   yy1       <- Data[all_pos1,1]
# n0       <- length(all_pos0);   n1        <- length(all_pos1)
# ij0      <- sample(all_pos0,2); ij1       <- sample(all_pos1,2)
# lab_i0   <- all_lab0[ij0[1]];   lab_i1    <- all_lab1[ij1[1]]
# lab_j0   <- all_lab0[ij0[2]];   lab_j1    <- all_lab1[ij1[2]]
# 
# 
# ############################################################################
# # PROCESS ZERO
# ############################################################################
# 
# 
# #######################################################
# ############################# SPLIT MOVE
# if(lab_i0==lab_j0){
#   current_lab               <- lab_i0
#   pos_with_current_lab      <- which(all_lab0 == current_lab)
#   fake_alllab0              <- all_lab0; fake_alllab0[ij0] <- -999
#   pos_with_current_lab_noij <- which(fake_alllab0==current_lab)
#   
#   ##############################################################
#   #if{codizione se i due indici sono soli}
#   ##############################################################
#   
#   LL0                       <- Split0(yy = Data[,1],indici_position_original_noij = pos_with_current_lab_noij,
#                                       i_ind1 = ij0[1],j_ind2 = ij0[2],sig = sigma0,a = a_0,b = b_0)
#   # Data[LL0$Si,5] <- current_lab ...should be useless
#   all_lab0[LL0$Sj] <- max(Data[,5])+1
#   n_Kappa            <- table(all_lab0[all_pos0])
#   n_Kappa_old        <- table(all_lab0old[all_pos0])
#   # qui serve MHratio per comparare all_lab0 con old vettore Data[,5]
#   acc_sp0 <- log_MH0_SPLIT(yy0 = yy0,labels0 = all_lab0old[all_pos0],
#                            labels0new = all_lab0[all_pos0],
#                            nk = n_Kappa_old,nknew = n_Kappa,
#                            a0 = a_0,b0 = b_0,t0 = theta0,s0 = sigma0,
#                            log_prob_split = LL0$cumlogprob)
#   
#   if(acc_sp0){ Data[,5] <- all_lab0 }
#   
#   #######################################################
#   ############################# MERGE MOVE
# }else{
#   pos_with_current_labs     <- which(all_lab0old == lab_i0 | all_lab0old == lab_j0)
#   fake_alllab0              <- all_lab0old;      fake_alllab0[ij0] <- -999
#   pos_with_current_labs_noij <- which(fake_alllab0 == lab_i0 | fake_alllab0 == lab_j0)
#   
#   
#   LL0 <- Merge0(yy = Data[,1],indici_position_original_noij = pos_with_current_labs_noij,
#                 indici_labels_original = all_lab0,
#                 i_ind1 = ij0[1], j_ind2 = ij0[2],
#                 i_lab1 = lab_i0, j_lab2 = lab_j0, sig = sigma0, a = a_0,b = b_0)
#   all_lab0[all_lab0==lab_j0] <- lab_i0
#   all_lab0                   <- reset(all_lab0)-1
#   n_Kappa                    <- table(all_lab0[all_pos0])
#   n_Kappa_old                <- table(all_lab0old[all_pos0])
#   # qui serve MHratio per comparare all_lab0 con old vettore Data[,5]
#   acc_sp0 <- log_MH0_MERGE(yy0        = yy0,
#                            labels0    = all_lab0old[all_pos0],
#                            labels0new = all_lab0[all_pos0],
#                            nk         = n_Kappa_old,
#                            nknew      = n_Kappa,
#                            a0 = a_0, b0 = b_0, t0 = theta0, s0 = sigma0,
#                            log_prob_split = LL0$cumlogprob)
#   
#   if(acc_sp0){
#     Data[,5] <- all_lab0
#   }
# }
# 
# ############################################################################
# # PROCESS ONE
# ############################################################################
# 
# #######################################################
# ############################# SPLIT MOVE
# if(lab_i1==lab_j1){
#   current_lab               <- lab_i1
#   pos_with_current_lab      <- which(all_lab1==current_lab)
#   fake_alllab1              <- all_lab1; fake_alllab1[ij1] <- -999
#   pos_with_current_lab_noij <- which(fake_alllab1==current_lab)
#   
#   ##############################################################
#   #        if{codizione se i due indici sono soli}
#   ##############################################################
#   
#   LL1                       <- Split1(yy = Data[,1],indici_position_original_noij = pos_with_current_lab_noij,
#                                       i_ind1 = ij1[1],j_ind2 = ij1[2],sig = sigma1,m1 = m1,V1 = V1,a1 =  a1,b1 = b1)
#   
#   all_lab1[LL1$Sj]   <- max(Data[,6])+1
#   
#   n_Kappa            <- table(all_lab1[all_pos1])
#   n_Kappa_old        <- table(all_lab1old[all_pos1])
#   acc_sp1 <- log_MH1_SPLIT(yy1 = yy1,labels1 = all_lab1old[all_pos1],
#                            labels1new = all_lab1[all_pos1],
#                            nk = n_Kappa_old,nknew = n_Kappa,
#                            m1 = m1,V1 = V1,
#                            a1 = a1,b1 = b1,t1 = theta1,s1 = sigma1,
#                            log_prob_split = LL1$cumlogprob)
#   
#   if(acc_sp1){
#     Data[,6] <- all_lab1
#   }
#   
#   #######################################################
#   ############################# MERGE MOVE
# }else{
#   pos_with_current_labs     <- which(all_lab1old == lab_i1 | all_lab1old == lab_j1)
#   fake_alllab1              <- all_lab1old;      fake_alllab1[ij1] <- -999
#   pos_with_current_labs_noij <- which(fake_alllab1 == lab_i1 | fake_alllab1 == lab_j1)
#   
#   
#   LL1 <- Merge1(yy = Data[,1],indici_position_original_noij = pos_with_current_labs_noij,
#                 indici_labels_original = all_lab1,
#                 i_ind1 = ij1[1], j_ind2 = ij1[2],
#                 i_lab1 = lab_i1, j_lab2 = lab_j1, sig = sigma1, m1 = m1, V1 = V1, a1 = a1, b1 = b1)
#   all_lab1[all_lab1==lab_j1] <- lab_i1
#   all_lab1                   <- reset(all_lab1)-1
#   n_Kappa                    <- table(all_lab1[all_pos1])
#   n_Kappa_old                <- table(all_lab1old[all_pos1])
#   # qui serve MHratio per comparare all_lab0 con old vettore Data[,5]
#   acc_sp1 <- log_MH1_MERGE(yy1        = yy1,
#                            labels1    = all_lab1old[all_pos1],
#                            labels1new = all_lab1[all_pos1],
#                            nk         = n_Kappa_old,
#                            nknew      = n_Kappa,
#                            m1=m1,V1=V1,a1 = a1, b1 = b1, t1 = theta1, s1 = sigma1,
#                            log_prob_split = LL1$cumlogprob)
#   if(acc_sp1){ Data[,6] <- all_lab1}
# }