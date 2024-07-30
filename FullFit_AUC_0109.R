require(mvtnorm)
require(osDesign)
require(plyr)
require(MESS)
require(mvtnorm)
beta0 <- -2.02
beta_X <- c(log(0.6), log(1.6), log(0.6))
gama_Z <- log(1.5)
beta<-c(beta0,beta_X,gama_Z)
set.seed(2023)
n <- 3000
rou <- 0
ratio<-1

replicate <- 100

Result.AUC <- matrix(0, nrow = length(rou), ncol = 4)
colnames(Result.AUC) <- c("rou",  "AUC0", "Emp.sd_AUC0", "Asy.sd_AUC0")

AUC0_values <- numeric(replicate)
Asy_sd_AUC0<-numeric(replicate)
  
for (ci in 1:replicate) {
  
  #generated full data
  sigma <- matrix(c(1, rou, rou, 1), nrow=2, byrow=TRUE) 
  x1z <-rmvnorm(n, mean=c(0, 0), sigma=sigma)
  x1 <- x1z[,1]
  z <- x1z[,2]
  
  x2 <- runif(n)
  x3 <- rbinom(n, size=1, prob=0.2)
  X <- cbind(1, x1, x2, x3, z)
  y <- rbinom(n, size=1, prob=exp(X %*% beta)/(1+exp(X %*% beta)))
  dat <- data.frame(cbind(y,x1,x2,x3,z))
  dat_case <- dat[dat$y==1, ]
  dat_con <- dat[dat$y==0, ]
  data<-dat #data include the full data
  
  m <- glm(y ~ x1 + x2 + x3 + z, family=binomial, data=data)
  yhat.full <- m$fitted.values #this is the yhat （non-0-1）
  beta.full <- m$coefficients
  beta.var.full <- summary(m)$cov.unscaled
  X <- cbind(rep(1,nrow(data)), data$x1, data$x2, data$x3, data$z) #Design matrix
  data <- cbind(data, yhat.full)
  
  #caculate the Asy.SD and the AUC
  # Score for each subject: case/control
  S1_beta <- apply(cbind(X, yhat.full), 1, function(x) x[1:5]*(1-x[6]))
  S0_beta <- apply(cbind(X, yhat.full), 1, function(x) x[1:5]*(0-x[6]))
  
  # inverse of -dS/dbeta*(1/n)
  I_betabeta <- solve(t(X) %*% diag(as.vector(yhat.full*(1-yhat.full))) %*% X/n)
  
  # influence function 
  h1 <- I_betabeta %*% S1_beta
  h0 <- I_betabeta %*% S0_beta
    

  #calculate auc
  #c <- seq(0, 1, by=0.05) #possible threshold values for ROC curve
  c <- c(seq(0, 1, by=0.01)*max(yhat.full),1)
  denom_tpr.full <- sum(yhat.full)
  denom_fpr.full <- sum(1-yhat.full)
  num_tpr.full <- numeric(length(c))
  num_fpr.full <- numeric(length(c))
  for (j in 1:length(c)){
   ind <- (yhat.full>=c[j])
   num_tpr.full[j] <- sum(yhat.full[ind==1])
   num_fpr.full[j] <- sum((1-yhat.full)[ind==1])
}
  tpr.full <- num_tpr.full/denom_tpr.full
  fpr.full <- num_fpr.full/denom_fpr.full
 
  AUC0_values[ci] <- auc(fpr.full, tpr.full)

# numerical method for dFPR/dbeta and dTPR/dbeta as a function of c
  numbeta <- length(beta)
  epsilon_beta <- abs(beta)/5

  dFPR_beta <- matrix(0, nrow=numbeta, ncol=length(c))
  dTPR_beta <- matrix(0, nrow=numbeta, ncol=length(c))
  yfit1 <- numeric(nrow(data))
  yfit2 <- numeric(nrow(data))
 for (i in 1:numbeta){
   tmp1 <- beta.full
   tmp2 <- beta.full
   tmp1[i] <- beta.full[i] + epsilon_beta[i]
   tmp2[i] <- beta.full[i] - epsilon_beta[i]
   yfit1 <- exp(X %*% tmp1)/(1 + exp(X %*% tmp1))
   yfit2 <-  exp(X %*% tmp2)/(1 + exp(X %*% tmp2))
   num1_tpr <- 0
   num2_tpr <- 0
   num1_fpr <- 0
   num2_fpr <- 0
   denom1_tpr <- sum(yfit1)
   denom2_tpr <- sum(yfit2)
   denom1_fpr <- sum(1-yfit1)
   denom2_fpr <- sum(1-yfit2)
  
  for (j in 1:length(c)){
    ind1 <- (yfit1>=c[j])
    ind2 <- (yfit2>=c[j])
    num1_tpr <- sum(yfit1[ind1==1])
    num2_tpr <- sum(yfit2[ind2==1])
    num1_fpr <- sum((1-yfit1)[ind1==1])
    num2_fpr <- sum((1-yfit2)[ind2==1])
    dTPR_beta[i,j] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*epsilon_beta[i])
    dFPR_beta[i,j] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*epsilon_beta[i])
  }
}

# calculate influence function for TPR/FPR as a function of c
 mat_tpr.full <- matrix(0, nrow=nrow(data), ncol=length(c))
 mat_fpr.full <- matrix(0, nrow=nrow(data), ncol=length(c))
 for (j in 1:length(c)){
   ind <- (yhat.full>=c[j])
   mat_tpr.full[,j] <- yhat.full*(ind-tpr.full[j])/(denom_tpr.full/n)
   mat_fpr.full[,j] <- (1-yhat.full)*(ind-fpr.full[j])/(denom_fpr.full/n)
}

# calculate the influence function for tpr/fpr
 h1_tpr <- t(h1) %*% dTPR_beta + mat_tpr.full
 h0_tpr <- t(h0) %*% dTPR_beta + mat_tpr.full
 h1_fpr <- t(h1) %*% dFPR_beta + mat_fpr.full
 h0_fpr <- t(h0) %*% dFPR_beta + mat_fpr.full

#calcualte influence function for auc
 h1_auc <- 0-apply(h1_fpr, 1, function(x) auc(tpr.full, x))+apply(h1_tpr, 1, function(x) auc(fpr.full, x))
 h0_auc <- 0-apply(h0_fpr, 1, function(x) auc(tpr.full, x))+apply(h0_tpr, 1, function(x) auc(fpr.full, x))

 emp_cov_auc <- sum((h1_auc)^2*data$y, (h0_auc)^2*(1-data$y))
 var_auc <- emp_cov_auc/n^2
 Asy_sd_AUC0[ci]<-sqrt(var_auc)

}

Result.AUC[,1]<-rou
Result.AUC[,2]<-mean(AUC0_values)
Result.AUC[,3]<-sd(AUC0_values)
Result.AUC[,4]<-mean(Asy_sd_AUC0)
Result.AUC<-data.frame(Result.AUC)
print('okk')