########################### R HELPERS

reset <- function(x){
  z <- y <- x
  sy <- unique(sort(y))
  for(i in 1:length(unique(x))){
    z[which(x==sy[i])] <- i
  }
  return(z)
}
############################################################################################
############################################################################################
############################################################################################
RESET_LAB0LAB1_GivenTHETA <- function(Data){
  n <- nrow(Data)
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
  return(Data)
}

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

SPLIT_MERGE_PROCESS0 <- function( Data,
                                  sigma0,
                                  theta0,
                                  a_0,
                                  b_0)   
{
  
  all_lab0old <- all_lab0 <- Data[,5]
  all_pos0old <- all_pos0 <- which(Data[,5]>0)
  yy0         <- Data[all_pos0,1]; 
  n0          <- length(all_pos0);   
  ij0         <- sample(all_pos0,2);
  lab_i0      <- all_lab0[ij0[1]];
  lab_j0      <- all_lab0[ij0[2]];
  
      ############################# SPLIT MOVE
if(lab_i0==lab_j0){
  current_lab               <- lab_i0
  #pos_with_current_lab      <- which(all_lab0 == current_lab)
  fake_alllab0              <- all_lab0; fake_alllab0[ij0] <- -999
  pos_with_current_lab_noij <- which(fake_alllab0==current_lab)
  
  ##############################################################
  #if{codizione se i due indici sono soli}
  ##############################################################
  
  LL0                       <- Split0(yy = Data[,1],indici_position_original_noij = pos_with_current_lab_noij,
                                      i_ind1 = ij0[1],j_ind2 = ij0[2],sig = sigma0,a = a_0,b = b_0)
  all_lab0[LL0$Sj]   <- max(Data[,5])+1
  n_Kappa            <- table(all_lab0[all_pos0])
  n_Kappa_old        <- table(all_lab0old[all_pos0])
  # qui serve MHratio per comparare all_lab0 con old vettore Data[,5]
  acc_sp0 <- log_MH0_SPLITnew(yy0 = yy0,
                           labels0 = all_lab0old[all_pos0],
                           labels0new = all_lab0[all_pos0],
                           nk = n_Kappa_old,
                           nknew = n_Kappa,
                           
                           labsplit = c(lab_i0,max(Data[,5])+1),
                           labmerge = lab_i0,
                           
                           a0 = a_0,b0 = b_0,t0 = theta0,s0 = sigma0,
                           log_prob_split = LL0$cumlogprob)
  if(acc_sp0){ Data[,5] <- all_lab0 }
  
  #######################################################
  ############################# MERGE MOVE
}else{
  
  #pos_with_current_labs      <- which(all_lab0old == lab_i0 | all_lab0old == lab_j0)
  fake_alllab0               <- all_lab0old;      fake_alllab0[ij0] <- -999
  pos_with_current_labs_noij <- which(fake_alllab0 == lab_i0 | fake_alllab0 == lab_j0)
  
  
  LL0 <- Merge0(yy = Data[,1],indici_position_original_noij = pos_with_current_labs_noij,
                indici_labels_original = all_lab0,
                i_ind1 = ij0[1], j_ind2 = ij0[2],
                i_lab1 = lab_i0, j_lab2 = lab_j0, sig = sigma0, a = a_0,b = b_0)
  all_lab0[all_lab0==lab_j0] <- lab_i0
  all_lab0                   <- reset(all_lab0)-1
  n_Kappa                    <- table(all_lab0[all_pos0])
  n_Kappa_old                <- table(all_lab0old[all_pos0])
  # qui serve MHratio per comparare all_lab0 con old vettore Data[,5]
  acc_sp0 <- log_MH0_MERGEnew(yy0        = yy0,
                           labels0    = all_lab0old[all_pos0],
                           labels0new = all_lab0[all_pos0],
                           nk         = n_Kappa_old,
                           nknew      = n_Kappa,
                           
                           labsplit = c(lab_i0,lab_j0),
                           labmerge = ifelse(lab_i0<lab_j0,lab_i0,lab_i0-1),
                           
                           a0 = a_0, b0 = b_0, 
                           t0 = theta0, s0 = sigma0,
                           log_prob_split = LL0$cumlogprob)
  
  if(acc_sp0){
    Data[,5] <- all_lab0
  }
}
#   print(acc_sp0)
return(Data)
  }  
  
  



SPLIT_MERGE_PROCESS1 <- function( Data,
                                  sigma1,
                                  theta1,
                                  m1, V1,
                                  a1,
                                  b1){
    
    all_lab1old <- all_lab1 <- Data[,6]
    all_pos1old <- all_pos1 <- which(Data[,6]>0)
    yy1       <- Data[all_pos1,1]
    ij1       <- sample(all_pos1,2)
    lab_i1    <- all_lab1[ij1[1]]
    lab_j1    <- all_lab1[ij1[2]]
    
    
    
    if(lab_i1==lab_j1){
      current_lab               <- lab_i1
      #pos_with_current_lab      <- which(all_lab1==current_lab)
      fake_alllab1              <- all_lab1; fake_alllab1[ij1] <- -999
      pos_with_current_lab_noij <- which(fake_alllab1==current_lab)
      
      ##############################################################
      #        if{codizione se i due indici sono soli}
      ##############################################################
      
      LL1                       <- Split1(yy = Data[,1],indici_position_original_noij = pos_with_current_lab_noij,
                                          i_ind1 = ij1[1],j_ind2 = ij1[2],sig = sigma1,m1 = m1,V1 = V1,a1 =  a1,b1 = b1)
      
      all_lab1[LL1$Sj]   <- max(Data[,6])+1
      
      n_Kappa            <- table(all_lab1[all_pos1])
      n_Kappa_old        <- table(all_lab1old[all_pos1])
      acc_sp1 <- log_MH1_SPLITnew(yy1 = yy1,labels1 = all_lab1old[all_pos1],
                               labels1new = all_lab1[all_pos1],
                               nk = n_Kappa_old,
                               nknew = n_Kappa,
                               labsplit = c(lab_i1,max(Data[,6])+1),
                               labmerge = lab_i1,
                               m1 = m1,V1 = V1,
                               a1 = a1,b1 = b1,t1 = theta1,s1 = sigma1,
                               log_prob_split = LL1$cumlogprob)
      
      if(acc_sp1){
        Data[,6] <- all_lab1
      }
      
      #######################################################
      ############################# MERGE MOVE
    }else{
      #pos_with_current_labs     <- which(all_lab1old == lab_i1 | all_lab1old == lab_j1)
      fake_alllab1              <- all_lab1old;      fake_alllab1[ij1] <- -999
      pos_with_current_labs_noij <- which(fake_alllab1 == lab_i1 | fake_alllab1 == lab_j1)
      
      
      LL1 <- Merge1(yy = Data[,1],indici_position_original_noij = pos_with_current_labs_noij,
                    indici_labels_original = all_lab1,
                    i_ind1 = ij1[1], j_ind2 = ij1[2],
                    i_lab1 = lab_i1, j_lab2 = lab_j1, sig = sigma1, m1 = m1, V1 = V1, a1 = a1, b1 = b1)
      all_lab1[all_lab1==lab_j1] <- lab_i1
      all_lab1                   <- reset(all_lab1)-1
      n_Kappa                    <- table(all_lab1[all_pos1])
      n_Kappa_old                <- table(all_lab1old[all_pos1])
      # qui serve MHratio per comparare all_lab0 con old vettore Data[,5]
      acc_sp1 <- log_MH1_MERGEnew(yy1     = yy1,
                               labels1    = all_lab1old[all_pos1],
                               labels1new = all_lab1[all_pos1],
                               nk         = n_Kappa_old,
                               nknew      = n_Kappa,#####################
                               labsplit   = c(lab_i1,lab_j1),
                               labmerge   = ifelse(lab_i1<lab_j1,lab_i1,lab_i1-1),
                               
                               m1=m1,V1=V1,a1 = a1, b1 = b1, t1 = theta1, s1 = sigma1,
                               log_prob_split = LL1$cumlogprob)
      if(acc_sp1){ Data[,6] <- all_lab1}
    }
    return(Data)
  }  



###########################################################################################
posterior_predictive <- function(Data, postpred, NSIM, ygrid=ygrid, rho=rho){
tempdens0 <- tempdens1 <- numeric(1000)
ne0 <- table(Data[,5])[-1]
ne1 <- table(Data[,6])[-1]
T0 <- cbind(uniquecombs(Data[Data[,4]==0,c(2,3,5)]),ne0)
if(length(ne1)==1){
  T1 <- c(uniquecombs(Data[Data[,4]==1,c(2,3,6)]),ne1) 
  NNN1 = 1
  tempdens1 <- tempdens1 +   dnorm(ygrid, T1[1],sqrt(T1[2]))
  }else{
  T1 <- cbind(uniquecombs(Data[Data[,4]==1,c(2,3,6)]),ne1)
  NNN1 = nrow(T1)
  for(qq in 1:NNN1){
    tempdens1 <- tempdens1 + ne1[qq]/sss1 * dnorm(ygrid, T1[qq,1],sqrt(T1[qq,2]))
  }
  }
sss0 <- sum(ne0); sss1 <- sum(ne1)
for(oo in 1:nrow(T0)){
  tempdens0 <- tempdens0 + ne0[oo]/sss0 * dnorm(ygrid, T0[oo,1],sqrt(T0[oo,2]))
}
postpred[,1]  = postpred[,1]  + tempdens0/NSIM
postpred[,2]  = postpred[,2]  + tempdens1/NSIM
postpred[,3]  = postpred[,3]  + (postpred[,1]*(1-rho) + rho*postpred[,2])/NSIM

return(postpred)
}








###########################################################################################
posterior_predictive2 <- function(Data, postpred, NSIM, ygrid=ygrid, rho=rho,t0,t1,s0,s1,a0,b0,a1,b1,m1,V1){
  
  tempdens0 <- tempdens1 <- numeric(1000)
  
  ne0 <- table(Data[,5])[-1]
  dk0 <- length(ne0)
  ne1 <- table(Data[,6])[-1]
  dk1 <- length(ne1)
  enne0 <- sum(ne0); enne1 <- sum(ne1)
 
  T0 <- cbind(uniquecombs(Data[Data[,4]==0,c(2,3,5)]),ne0)
  
  if(length(ne1)==1){
    T1 <- c(uniquecombs(Data[Data[,4]==1,c(2,3,6)]),ne1) 
    NNN1 = 1
    tempdens1 <- tempdens1 +  (ne1-s1)/(enne1+t1) * dnorm(ygrid, T1[1],sqrt(T1[2]))
  }else{
    T1 <- cbind(uniquecombs(Data[Data[,4]==1,c(2,3,6)]),ne1)
    NNN1 = nrow(T1)
    for(qq in 1:NNN1){
      tempdens1 <- tempdens1 + (ne1[qq]-s1)/(enne1+t1) * dnorm(ygrid, T1[qq,1],sqrt(T1[qq,2])) 
    }
  }
  # ignoro contributo prior, che rallenta di brutto
  #tempdens1 <- tempdens1 + (t1 + dk1*s1)/(enne1+t1) * exp(log_Marginal_P1_grid(grid = ygrid,m1,V1,a1,b1))
  
  for(oo in 1:nrow(T0)){
    tempdens0 <- tempdens0 + ne0[oo]/(enne0+t0) * dnorm(ygrid, T0[oo,1],sqrt(T0[oo,2]))
  }
  #tempdens0 <- tempdens0 + (t0 + dk0*s0)/(enne0+t0) * exp(log_Marginal_P0_grid(ygrid,a0,b0))
  postpred[,1]  = postpred[,1]  + tempdens0/NSIM
  postpred[,2]  = postpred[,2]  + tempdens1/NSIM
  postpred[,3]  = postpred[,3]  + (postpred[,1]*(1-rho) + rho*postpred[,2])/NSIM
  
  return(postpred)
}



####################################################################
