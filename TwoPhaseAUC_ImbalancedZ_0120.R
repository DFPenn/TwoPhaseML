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

Result.AUC <- matrix(0, nrow = length(rou), ncol = 8)
colnames(Result.AUC) <- c("rou",  "AUC1", "Emp.sd_AUC1", "Asy.sd_AUC1","count_Z0","count_Z1","count_Z2","count_Z3")


# Initialize variables to store AUC values and standard deviation of AUC values for each replicate
AUC1_values <- numeric(replicate)
Asy_sd_AUC1<-numeric(replicate)
Phase2_case_counts <- numeric(replicate)
Phase2_noncase_counts <- numeric(replicate)
count_z0<-numeric(replicate)
count_z1<-numeric(replicate)
count_z2<-numeric(replicate)
count_z3<-numeric(replicate)

for (ci in 1:replicate) {
  
  #generated full data
  sigma <- matrix(c(1, rou, rou, 1), nrow=2, byrow=TRUE) 
  x1z <-rmvnorm(n, mean=c(0, 0), sigma=sigma)
  x1 <- x1z[,1]
  z1 <- x1z[,2]
  z <- cut(z1, breaks = c(-Inf, -1.2, 0, 1.2, Inf), labels = c(0, 1, 2, 3), include.lowest = TRUE)
  z<- as.numeric(as.character(z))
  
  #quantiles <- quantile(z1, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  #z<- cut(z1, breaks = quantiles, labels =c(0,1,2,3), include.lowest = TRUE)
  #z<- as.numeric(as.character(z))
  
  count_z0[ci]<-sum(z==0)
  count_z1[ci]<-sum(z==1)
  count_z2[ci]<-sum(z==2)
  count_z3[ci]<-sum(z==3)
  
  x2 <- runif(n)#用于返回0-1之间的伪随机数
  x3 <- rbinom(n, size=1, prob=0.2)
  X <- cbind(1, x1, x2, x3, z)
  y <- rbinom(n, size=1, prob=exp(X %*% beta)/(1+exp(X %*% beta)))
  dat <- data.frame(cbind(y,x1,x2,x3,z))
  dat_case <- dat[dat$y==1, ]
  dat_con <- dat[dat$y==0, ]
  data<-dat
  
  #更极端的情况1 这一组是把z=0和z=3的全部抽出来，然后z=1，2的平分
  #dat2 <- data.frame() 
  # Extract second-stage data
  #indices_z0 <- sample(which(dat$z == 0), replace = FALSE)
  #indices_z1 <- sample(which(dat$z == 1), size = (900-sum(z==0)-sum(z==3))/2, replace = FALSE)
  #indices_z2 <- sample(which(dat$z == 2), size = (900-sum(z==0)-sum(z==3))/2, replace = FALSE)
  #indices_z3 <- sample(which(dat$z == 3), replace = FALSE)
  #indices<- c(indices_z0, indices_z1, indices_z2, indices_z3)
  #dat2 <- rbind(dat2, dat[indices, ])
  
  #更极端的情况2
  #all case in Z size
  #200，20，20，200 VS 20，200，200，20
  #100，10，10，200 VS 10，100，200，10
  #10，650，190，50
  #20，200，660，200
  #5，500，390，5
  dat2 <- data.frame() 
  indices_z0 <- sample(which(dat$z == 0), size=100,replace = FALSE)
  indices_z1 <- sample(which(dat$z == 1), size=10, replace = FALSE)
  indices_z2 <- sample(which(dat$z == 2), size=20 , replace = FALSE)
  indices_z3 <- sample(which(dat$z == 3), size=100,replace = FALSE)
  indices<- c(indices_z0, indices_z1, indices_z2, indices_z3)
  dat2 <- rbind(dat2, dat[indices, ])
  
 
  #generate Phase II data SRS
  #dat2_indices <- sample(1:n, size = round(n * 0.15 * 2))#默认无放回抽样
  #dat2 <- dat[dat2_indices, ]
  
  #balanced
  #dat2 <- data.frame()
  #for (z_value in unique(z)) {
  #  z_indices <- sample(which(z == z_value), size = 225, replace = FALSE)
  #  dat2 <- rbind(dat2, dat[z_indices, ])
  #}
  
  Phase2_case_counts[ci] <- sum(dat2$y)
  Phase2_noncase_counts[ci] <- sum(1-dat2$y)
  
  #合并数据
  nonDeath <- aggregate(1-dat$y, by=list(x1=dat$x1, x2=dat$x2, x3=dat$x3, z=dat$z), FUN=sum)$x
  dat1.mle <- data.frame(aggregate(dat$y, by=list(x1=dat$x1, x2=dat$x2, x3=dat$x3, z=dat$z), FUN=sum))
  dat1.mle <- data.frame(cbind(dat1.mle, nonDeath))
  names(dat1.mle) <- c("x1","x2","x3","z", "Death","nonDeath")
  
  conts <- aggregate(1-dat2$y, by=list(x1=dat2$x1, x2=dat2$x2, x3=dat2$x3, z=dat2$z), FUN=sum)$x
  dat2.mle <- data.frame(aggregate(dat2$y, by=list(x1=dat2$x1, x2=dat2$x2, x3=dat2$x3, z=dat2$z), FUN=sum))
  dat2.mle <- data.frame(cbind(dat2.mle, conts))
  names(dat2.mle) <- c("x1","x2","x3","z","cases","conts")
  
  # Final = Phase I + Phase II
  fdat.mle <- merge(dat1.mle, dat2.mle, by=c("x1","x2","x3","z"), all=T)
  fdat.mle$cases <- ifelse(is.na(fdat.mle$cases), 0, fdat.mle$cases)
  fdat.mle$conts <- ifelse(is.na(fdat.mle$conts), 0, fdat.mle$conts)
  
  # Estimate u_i = n_i/N_i
  nn1 <- sum(dat$y)
  nn0 <- sum(1-dat$y)
  n1 <- sum(dat2$y)
  n0 <- sum(1-dat2$y)
  
  w1 <- round(n1/nn1, 3)
  w0 <- round(n0/nn0, 3)
  mu_1 <- rep(w1, nrow(fdat.mle))
  mu_0 <- rep(w0, nrow(fdat.mle))
  
  #估计系数
  #Fit Two-Phase maximum likelihood model
  case=fdat.mle$cases
  
  x=cbind(rep(1,nrow(fdat.mle)),fdat.mle$x1,fdat.mle$x2,fdat.mle$x3,fdat.mle$z)
  N=fdat.mle$cases+fdat.mle$conts
  
  mod<- tps(cbind(cases, conts) ~ x1+x2+x3+z, data=fdat.mle, nn0=nn0, nn1=nn1, group=rep(1,nrow(fdat.mle)), method="ML", cohort=T)
  beta.mle <- mod$coef
  XX2 <- cbind(rep(1,nrow(dat2.mle)), dat2.mle$x1, dat2.mle$x2, dat2.mle$x3, dat2.mle$z)
  yhat2.mle <- exp(XX2 %*% beta.mle)/(1 + exp(XX2 %*% beta.mle))
  gamma_i <- matrix(NA, 2, 1)
  
  gamma_i[1,] <- 1
  for(itr in 1:40)
  {
    mm1 <- (n1-gamma_i[1])/(nn1-gamma_i[1,])
    mm0 <- (n0+gamma_i[1,])/(nn0+gamma_i[1,])
    gamma_i[1,] <- n1-sum(yhat2.mle*mm1/(mm0*(1-yhat2.mle)+yhat2.mle*mm1))
    itr=itr+1
  }
  gamma_i[2,] <- -gamma_i[1,]
  
  Q_i <- matrix(NA, 2, 1)
  
  Q_i[1,] <- (nn1-gamma_i[1,])/(nn0+nn1)
  Q_i[2,] <- (nn0-gamma_i[2,])/(nn0+nn1)
  
  dat11.mle <- data.frame(count(dat, vars=c("y")))
  
  strata_mat <- matrix(NA, 4, 1)
  
  strata_mat[1,] <- n1
  strata_mat[2,] <- n0
  strata_mat[3,] <- nn1-n1
  strata_mat[4,] <- nn0-n0
  
  u_i <- matrix(NA, 2, 1)
  for (i in 1:2){
    u_i[i, 1] <- 1 - strata_mat[i+2, 1]/(sum(strata_mat[,1])*Q_i[i, 1])
  }
  
  
  numbeta <- length(beta)
  epsilon_beta <- abs(beta)/5
  # All things that we need to calculate based on above information obtained 
  W <- 0
  B <- matrix(NA, numbeta, 1)
  A <- matrix(NA, 2, 1)
  
  counts <- dat2.mle$cases+dat2.mle$conts
  x1 <- dat2.mle$x1
  x2 <- dat2.mle$x2
  x3 <- dat2.mle$x3
  z <- dat2.mle$z
  X <- cbind(rep(1,length(x1)),x1,x2,x3,z) 
  u <- c(u_i[1, 1], u_i[2, 1])
  ystar.mle <- u[1]*yhat2.mle/(u[1]*yhat2.mle+u[2]*(1-yhat2.mle))
  
  #W
  W <- sum(counts*(ystar.mle*(1-ystar.mle)))
  #B
  B[, 1] <- t(X) %*% (counts*ystar.mle*(1-ystar.mle))
  #A
  A[, 1] <- c(1/(strata_mat[1,1]-gamma_i[1,1])-1/(strata_mat[1,1]+strata_mat[3,1]-gamma_i[1,1]),
              1/(strata_mat[2,1]-gamma_i[2,1])-1/(strata_mat[2,1]+strata_mat[4,1]-gamma_i[2,1]))    
  
  A_0 <- apply(A, 2, sum)
  K <- 1/A_0-W
  dgamma_1 <- matrix(NA, numbeta, 1)
  dgamma_1[,1] <- -B[,1]/(1-A_0*W)
  
  
  #covariate distribution NPMLE estimates for each strata: f(x1,x2,x3,z)
  #since x1,x2 are continuous, each unique value will only be counted once
  #f(x1,x2,x3,z) is a function of beta
  f.mle <- numeric(nrow(dat2.mle))
  df <- matrix(NA, numbeta, nrow(dat2.mle)) #df(x1,x2,x3,z)/dbeta
  s_2 <- matrix(NA, numbeta, nrow(dat2.mle)) #score for R_k=1
  for (k in 1:nrow(dat2.mle)){
    y <- dat2.mle[k, "cases"]
    x1 <- dat2.mle[k, "x1"]
    x2 <- dat2.mle[k, "x2"]
    x3 <- dat2.mle[k, "x3"]
    z <- dat2.mle[k, "z"]
    counts <- dat2.mle[k, "cases"]+dat2.mle[k, "conts"]
    u <- c(u_i[1,1], u_i[2,1])
    yfit.mle <- yhat2.mle[k]
    dgamma <- dgamma_1[, 1]
    a <- A[,1]
    
    f.mle[k] <- counts/(sum(strata_mat[,1])*(u[1]*yfit.mle+u[2]*(1-yfit.mle)))
    df[,k] <- -counts/(sum(strata_mat[,1])*(u[1]*yfit.mle+u[2]*(1-yfit.mle))^2)*(c(1,x1,x2,x3,z)*yfit.mle*(1-yfit.mle)*(u[1]-u[2])+dgamma*((1-yfit.mle)*u[2]*a[2]-yfit.mle*u[1]*a[1]))
    s_2[, k] <- c(1,x1,x2,x3,z)*(y-yfit.mle) + df[, k]/f.mle[k]
  }
  
  s_1 <- matrix(NA, numbeta, nrow(dat11.mle)) #score for R_k=0
  for (k in 1:nrow(dat11.mle)){
    y <- dat1.mle[k, "Death"]
    Q <- NULL
    if (y==1){
      Q = Q_i[1, 1]
      s_1[, k] <- -dgamma_1[, 1]/sum(strata_mat[,1])/Q
    }
    else{
      Q = Q_i[2, 1]
      s_1[,k] <- dgamma_1[, 1]/sum(strata_mat[,1])/Q
    }
  }
  
  ## Influence function for odds ratios
  h_1 <- mod$covm %*% s_1*n
  h_2 <- mod$covm %*% s_2*n
  
  emp_cov <- matrix(0, numbeta, numbeta)
  counts2 <- c(strata_mat[4,], strata_mat[3,])
  for (k in 1:ncol(h_1)){
    emp_cov = emp_cov + h_1[, k] %*% t(h_1[, k])*counts2[k]
  }
  for (k in 1:ncol(h_2)){
    emp_cov = emp_cov + h_2[, k] %*% t(h_2[, k])*(dat2.mle$cases[k]+dat2.mle$conts[k])
  }
  # check emp_cov/n^2 = mod$covm
  
  
  #calculate auc
  
  beta.mle=beta.mle
  beta.var.mle=diag(mod$covm)
  h1=h_1
  h2=h_2
  yhat.mle=yhat2.mle
  fhat.mle=f.mle
  df_beta=df
  strata_mat=strata_mat
  
  # calculate auc
  denom_tpr.mle <- sum(yhat.mle*fhat.mle)
  denom_fpr.mle <- sum((1-yhat.mle)*fhat.mle)
  c <- c(seq(0, 1, by=0.01)*max(yhat.mle),1)
  num_tpr.mle <- numeric(length(c))
  num_fpr.mle <- numeric(length(c))
  for (j in 1:length(c)){
    ind <- (yhat.mle>=c[j])
    num_tpr.mle[j] <- sum((yhat.mle*fhat.mle)[ind==1])
    num_fpr.mle[j] <- sum(((1-yhat.mle)*fhat.mle)[ind==1])
  }
  tpr.mle <- num_tpr.mle/denom_tpr.mle
  fpr.mle <- num_fpr.mle/denom_fpr.mle
  AUC1_values[ci] <- auc(fpr.mle, tpr.mle)
  
  # numerical method for dFPR/dbeta and dTPR/dbeta as a function of c
  dFPR_beta <- matrix(0, nrow=numbeta, ncol=length(c))
  dTPR_beta <- matrix(0, nrow=numbeta, ncol=length(c))
  yfit1 <- numeric(nrow(dat2.mle))
  yfit2 <- numeric(nrow(dat2.mle))
  fhat1 <- numeric(nrow(dat2.mle))
  fhat2 <- numeric(nrow(dat2.mle))
  X <- cbind(rep(1,nrow(dat2.mle)), dat2.mle$x1, dat2.mle$x2, dat2.mle$x3, dat2.mle$z) #Design matrix
  for (i in 1:numbeta){
    tmp1 <- beta.mle
    tmp2 <- beta.mle
    tmp1[i] <- beta.mle[i] + epsilon_beta[i]
    tmp2[i] <- beta.mle[i] - epsilon_beta[i]
    yfit1 <- exp(X %*% tmp1)/(1 + exp(X %*% tmp1))
    yfit2 <-  exp(X %*% tmp2)/(1 + exp(X %*% tmp2))
    fhat1 <- fhat.mle + df_beta[i,]*epsilon_beta[i]
    fhat2 <- fhat.mle - df_beta[i,]*epsilon_beta[i]
    
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yfit1*fhat1)
    denom2_tpr <- sum(yfit2*fhat2)
    denom1_fpr <- sum((1-yfit1)*fhat1)
    denom2_fpr <- sum((1-yfit2)*fhat2)
    
    for (j in 1:length(c)){
      ind1 <- (yfit1>=c[j])
      ind2 <- (yfit2>=c[j])
      num1_tpr <- sum((yfit1*fhat1)[ind1==1])
      num2_tpr <- sum((yfit2*fhat2)[ind2==1])
      num1_fpr <- sum(((1-yfit1)*fhat1)[ind1==1])
      num2_fpr <- sum(((1-yfit2)*fhat2)[ind2==1])
      dTPR_beta[i,j] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*epsilon_beta[i])
      dFPR_beta[i,j] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*epsilon_beta[i])
    }
  }
  
  # calculate influence function for TPR/FPR as a function of c
  mat_tpr_h2.mle <- matrix(0, nrow=nrow(dat2.mle), ncol=length(c))
  mat_fpr_h2.mle <- matrix(0, nrow=nrow(dat2.mle), ncol=length(c))
  for (j in 1:length(c)){
    ind <- (yhat.mle>=c[j])
    mat_tpr_h2.mle[,j] <- fhat.mle*n*yhat.mle*(ind-tpr.mle[j])/denom_tpr.mle
    mat_fpr_h2.mle[,j] <- fhat.mle*n*(1-yhat.mle)*(ind-fpr.mle[j])/denom_fpr.mle
  }
  
  # calculate the influence function for tpr/fpr
  h1_tpr <- t(h1) %*% dTPR_beta
  h2_tpr <- t(h2) %*% dTPR_beta + mat_tpr_h2.mle
  h1_fpr <- t(h1) %*% dFPR_beta
  h2_fpr <- t(h2) %*% dFPR_beta + mat_fpr_h2.mle
  
  # calculate the influence function of auc
  h1_auc <- 0-apply(h1_fpr, 1, function(x) auc(tpr.mle, x))+apply(h1_tpr, 1, function(x) auc(fpr.mle, x))
  h2_auc <- 0-apply(h2_fpr, 1, function(x) auc(tpr.mle, x))+apply(h2_tpr, 1, function(x) auc(fpr.mle, x))
  
  # calcualte the variance of auc
  a <- rep(h1_auc, c(strata_mat[4, ], strata_mat[3, ]))
  b <- rep(h2_auc, (dat2.mle$cases+dat2.mle$conts))
  tmp <- c(a, b)
  var_auc <- t(tmp) %*% (tmp)/n^2
  Asy_sd_AUC1[ci]<-sqrt(var_auc)
  
}

countZ0<-sum(dat2$z==0)
countZ1<-sum(dat2$z==1)
countZ2<-sum(dat2$z==2)
countZ3<-sum(dat2$z==3)


Result.AUC[,1]<-rou
Result.AUC[,2]<-mean(AUC1_values)
Result.AUC[,3]<-sd(AUC1_values)
Result.AUC[,4]<-mean(Asy_sd_AUC1)
Result.AUC[,5] <- as.integer(mean(count_z0))
Result.AUC[,6] <- as.integer(mean(count_z1))
Result.AUC[,7] <- as.integer(mean(count_z2))
Result.AUC[,8] <- as.integer(mean(count_z3))
Result.AUC<-data.frame(Result.AUC)
print('okk')


mean(Phase2_case_counts)
mean(Phase2_noncase_counts)


library(openxlsx)
write.xlsx(Result.AUC, "E:/table/AUC_Imbalanced_Z1.xlsx", colNames = TRUE, rowNames = FALSE)
