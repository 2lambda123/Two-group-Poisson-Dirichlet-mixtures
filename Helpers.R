require(latex2exp)
#######################################################################
RHS <- function(k,p,alpha){
  fdr <- sum((1-p)*(p>k))/sum(p>k) 
  fdr-alpha
}
choose_threshold <- function(prob,AL){
  uniroot(function(x) RHS(x,p = prob,alpha = AL), interval = c(min(prob)+.00001,max(prob)-.00001))$root
}
##################################################################
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
###################################################################
plot.BNPtesting <- function(x,TH=.2){
  require("tidyverse")
  mppi   <- colMeans(x$ZiCon)
  thr    <- choose_threshold(mppi,TH)
  postz  <- ifelse(mppi>thr,1,0)
  D      <- cbind(id=1:length(x$Y),y=x$Y,mppi=mppi,postz=postz)
  D      <- as.tibble(D)
  
  p1    <- ggplot()+geom_line(data=D,aes(x=y, y=mppi),col="grey")+
    geom_point(data=D,aes(x=y, y=mppi,col=as.factor(postz+1)),alpha=.5)+
    theme_bw()+
    geom_hline(yintercept=thr,alpha=1,col=4,linetype=2)+ylab(unname(TeX("$P\\left(\\gamma_i=1|data\\right)$")))+xlab("z")+
    scale_color_manual(guide=F,values=c(2,3))
  p1
}


plot.BNPtesting <- function(x,TH=.2){
  require("tidyverse")
  mppi   <- colMeans(x$ZiCon)
  thr    <- choose_threshold(mppi,TH)
  postz  <- ifelse(mppi>thr,1,0)
  D      <- cbind(id=1:length(x$Y),y=x$Y,mppi=mppi,postz=postz)
  D      <- as.tibble(D)
  
  p1    <- ggplot()+geom_line(data=D,aes(x=y, y=mppi),col="grey")+
    geom_point(data=D,aes(x=y, y=mppi,col=as.factor(postz+1)),alpha=.5)+
    theme_bw()+
    geom_hline(yintercept=thr,alpha=1,col=4,linetype=2)+ylab(unname(TeX("$P\\left( \\gamma_i=1|data\\right)$")))+xlab("z")+
    scale_color_manual(guide=F,values=c(2,3))
  p1
}


Ranking <- function(x,names,TH=.2){
  require("knitr")
  mppi   <- colMeans(x$ZiCon)
  thr    <- choose_threshold(mppi,TH)
  postz  <- ifelse(mppi>thr,1,0)
  a      <- data.frame(Observation = names[which(postz==1)], PPI = mppi[which(postz==1)])
  a      <- a[sort(a[,2],index=T,decreasing = T)$ix,]
  kable(a)
}



plot.BNPtestingSM <- function(x,TH=.2){
  require("tidyverse")
  mppi   <- x$mppi
  thr    <- choose_threshold(mppi,TH)
  postz  <- ifelse(mppi>thr,1,0)
  D      <- cbind(id=1:length(x$Y),y=x$Y,mppi=mppi,postz=postz)
  D      <- as.tibble(D)
  
  p1    <- ggplot()+geom_line(data=D,aes(x=y, y=mppi),col="grey")+
    geom_point(data=D,aes(x=y, y=mppi,col=as.factor(postz+1)),alpha=.5)+
    theme_bw()+
    geom_hline(yintercept=thr,alpha=1,col=4,linetype=2)+ylab(unname(TeX("$P\\left( \\gamma_i=1|data\\right)$")))+xlab("z")+
    scale_color_manual(guide=F,values=c(2,3))
  p1
}


Ranking.BNPtestingSM<- function(x,names,TH=.2){
  require("knitr")
  mppi   <- x$mppi
  thr    <- choose_threshold(mppi,TH)
  postz  <- ifelse(mppi>thr,1,0)
  a      <- data.frame(Observation = names[which(postz==1)], PPI = mppi[which(postz==1)])
  a      <- a[sort(a[,2],index=T,decreasing = T)$ix,]
  kable(a)
}