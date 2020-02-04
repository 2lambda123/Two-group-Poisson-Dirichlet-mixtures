library(tidyverse)
library(mltools)
library(pROC)
library(MLmetrics)
library(modEvA)
library(latex2exp)
library(knitr)

###################################################################
  ck <- function(k,p){
    bg  <- 1-p
    ind <- (bg <= k)
    return( sum(bg*ind) )
  }
  
  ck_threshold <- function(KK,prob,alpha){
    n <- length(prob)
    ck(KK,prob)/n-alpha
  }

  bfdr <- function(k,p,alpha){
    fdr <- sum((1-p)*(p>k))/sum(p>k) 
    fdr-alpha
  }
  
  choose_threshold <- function(prob,AL){
    uniroot(function(x) bfdr(x,p = prob,alpha = AL), interval = c(min(prob)+.00001,max(prob)-.00001))$root
  }
  
  calculate_statistics <- function(tbl) { # from https://gist.github.com/kieranrcampbell/09b9d8ffdbac7582742d7a5896c3cd85
    P <- sum(tbl[,2])
    N <- sum(tbl[,1])
    TP <- tbl[2,2]
    TN <- tbl[1,1]
    FP <- tbl[2,1]
    FN  <- tbl[1,2]
    TPR <- TP / P
    FPR <- FP / N
    FDR <- FP / (TP + FP)
    ACC <- (TP + TN)/ sum(tbl)
    RECALL    <- TPR
    PRECISION <- TP/sum(tbl[2,])
    F1 <- 2/(1/RECALL+1/PRECISION)
    data.frame(P = P, N = N, TP = TP, TN = TN, FP = FP, 
               FN = FN, TPR = TPR, FPR = FPR, FDR = FDR, ACC =ACC, PREC=PRECISION, F1=F1)
  }
  
  l1 <- function(R){
    plot(R$RhoCon,type="l")
  }
  
  
  l2 <- function(R,Y,TrueZ,TH){
    par(mfrow=c(1,2))
    mppi <- colMeans(R$ZiCon)
    thr <- choose_threshold(mppi,TH)
    plot(mppi,type="l")+abline(h=thr,col=4,lwd=2,lty=2)
    poster_zi <- ifelse(mppi>thr,1,0)
    plot(Y,bg=poster_zi+2,pch=21+TrueZ)
    tab <- table(poster_zi,TrueZ)
    calculate_statistics(tab)
    par(mfrow=c(1,1))
    return(tab)
  }
  

  
  
  acc <- function(preds,actual) mean(preds==actual)
  
  compute_stat2 <- function(A,P,mppi){ #actual, #pred
    PRECISION   = Precision(y_true = A,y_pred = P,positive = 1)
    RECALL      = Recall(y_true = A,y_pred = P,positive = 1)
    SPECIFICITY = Specificity(y_true = A, y_pred = P,positive = 1)
    s=c(MCC = mcc(preds = P,actuals = A), 
        F1  = 2/(1/RECALL+1/PRECISION),
        SPECIFICITY=SPECIFICITY,
        ACC=acc(P,A),
        PRE = PRECISION,
        Auc = modEvA::AUC(obs = A, pred =mppi, simplif = T)
        )
        return(s)
  }
  
  SUMMA2 <- function(MA,TZ,TH){
    KKK <- matrix(NA,6,30)
    ROC <- array(NA,c(1001,2,30))
    for(i in 1:30){
      mppi      <- colMeans(MA[,i]$ZiCon)
      thr       <- choose_threshold(mppi,TH)
      cat(i)
      P         <- ifelse(mppi>thr,1,0)
      KKK[,i]   <- compute_stat2(A = TZ,P = P,mppi = mppi)
      ss  <- modEvA::AUC(obs = TZ, pred =mppi,interval = .001,plot = F)$th
      }
    return(list(KKK,ROC)) 
  }
  
  
