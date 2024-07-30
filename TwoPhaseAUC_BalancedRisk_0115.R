require(mvtnorm)
require(osDesign)
require(plyr)
require(MESS)
require(mvtnorm)
beta0 <- -3.52
beta_X <- c(log(0.6), log(1.6), log(0.6))
gama_Z <- log(1.5)
beta<-c(beta0,beta_X,gama_Z)
set.seed(2023)
n <- 3000
rou <- 0
ratio<-1

replicate <- 100

Result.AUC <- matrix(0, nrow = length(rou), ncol = 5)
colnames(Result.AUC) <- c("rou", "Beta0" , "AUC3", "Emp.sd_AUC3", "Asy.sd_AUC3")

# Initialize variables to store AUC values and standard deviation of AUC values for each replicate
AUC3_values <- numeric(replicate)
Asy_sd_AUC3<-numeric(replicate)


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
  ddat <- data.frame(y,x1,x2,x3,z)
  glm.fit <- glm(y~1+x1+x2+x3,data=ddat,family = binomial(link = "logit"))
  theta <- glm.fit$coefficients
  
  prob <- exp(X[,-5] %*% theta)/(1+exp(X[,-5] %*% theta))
  Quan <- unname(quantile(prob[ddat$y==1],c(0.25,0.5,0.75)))
  #Quan
  G <- sapply(prob, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
              else x<-4)
  dat <- data.frame(cbind(y,x1,x2,x3,z,G))
  dat_case <- dat[dat$y==1, ]
  dat_con <- dat[dat$y==0, ]
  data<-dat
  
  numbeta <- length(beta)
  numG <- 4
  numbeta <- length(beta)
  epsilon_beta <- abs(beta)/5
  
  #Generate Phase II data with E-balanced Design
  #n: total number of subjects in phase I
  #beta: odds ratio
  #numG: number of stratum defined by variables x1,x2,x3
  #p0: pre-specified probability of being selected into phase II for Y=0,G=g
  #ratio: pre-spefified number of controls relative to cases
  dat2_con <- NULL
  x0 <- xtabs(~dat$y+dat$G)[2,]*ratio
  p0 <- round(x0/xtabs(~dat$y+dat$G)[1,],3)
  for (i in 1:numG){
    dat_con_i <- dat_con[dat_con$G==i, ]
    selected_con <- 0
    tot_con <- 1
    ## controls for stratum G=i
    while (tot_con <= nrow(dat_con_i) & selected_con < x0[i]){
      choice <- rbinom(1, size=1, prob=p0[i])
      if (choice==1){
        dat2_con <- rbind(dat2_con, dat_con_i[tot_con, ])
        selected_con <- selected_con + 1
      }
      tot_con <- tot_con + 1
    }
  }
  dat2 <- rbind(dat_case, dat2_con)
  
  #合并数据 formatG
  nonDeath <- aggregate(1-dat$y, by=list(G=dat$G, x1=dat$x1, x2=dat$x2, x3=dat$x3, z=dat$z), FUN=sum)$x
  dat1.mle <- data.frame(aggregate(dat$y, by=list(G=dat$G, x1=dat$x1, x2=dat$x2, x3=dat$x3, z=dat$z), FUN=sum))
  dat1.mle <- data.frame(cbind(dat1.mle, nonDeath))
  names(dat1.mle) <- c("G","x1","x2","x3","z", "Death","nonDeath")
  
  conts <- aggregate(1-dat2$y, by=list(G=dat2$G, x1=dat2$x1, x2=dat2$x2, x3=dat2$x3, z=dat2$z), FUN=sum)$x
  dat2.mle <- data.frame(aggregate(dat2$y, by=list(G=dat2$G, x1=dat2$x1, x2=dat2$x2, x3=dat2$x3, z=dat2$z), FUN=sum))
  dat2.mle <- data.frame(cbind(dat2.mle, conts))
  names(dat2.mle) <- c("G","x1","x2","x3","z","cases","conts")
  
  # Final = Phase I + Phase II
  fdat.mle <- merge(dat1.mle, dat2.mle, by=c("G","x1","x2","x3","z"), all=T)
  fdat.mle$cases <- ifelse(is.na(fdat.mle$cases), 0, fdat.mle$cases)
  fdat.mle$conts <- ifelse(is.na(fdat.mle$conts), 0, fdat.mle$conts)
  
  # Estimate u_ig = n_ig/N_ig
  nn1 <- aggregate(dat$y, by=list(G=dat$G), FUN=sum)$x
  nn0 <- aggregate(1-dat$y, by=list(G=dat$G), FUN=sum)$x
  n1 <- aggregate(dat2$y, by=list(G=dat2$G), FUN=sum)$x
  n0 <- aggregate(1-dat2$y, by=list(G=dat2$G), FUN=sum)$x
  
  w1 <- round(n1/nn1, 3)
  w0 <- round(n0/nn0, 3)
  mu_1 <- sapply(fdat.mle$G, function(x) w1[x])
  mu_0 <- sapply(fdat.mle$G, function(x) w0[x])
  
  #估计theta mleG
  case=fdat.mle$cases
  group=fdat.mle$G
  x=cbind(rep(1,nrow(fdat.mle)),fdat.mle$x1,fdat.mle$x2,fdat.mle$x3,fdat.mle$z)
  N=fdat.mle$cases+fdat.mle$conts
  
  mod <- tps(cbind(cases, conts) ~ x1+x2+x3+z, data=fdat.mle, nn0=nn0, nn1=nn1, group=fdat.mle$G, method="ML", cohort=T)
  #beta.mle <- mod$coef
  beta.mle <- mod$coef
  #XX <- cbind(rep(1,nrow(fdat.mle)), fdat.mle$x1, fdat.mle$x2, fdat.mle$x3, fdat.mle$z) # Design matrix
  #yhat.mle <- exp(XX %*% beta.mle)/(1 + exp(XX %*% beta.mle))
  XX2 <- cbind(rep(1,nrow(dat2.mle)), dat2.mle$x1, dat2.mle$x2, dat2.mle$x3, dat2.mle$z)
  yhat2.mle <- exp(XX2 %*% beta.mle)/(1 + exp(XX2 %*% beta.mle))
  gamma_ig <- matrix(NA, 2, numG)
  for(g in 1:numG)
  {
    gamma_ig[1,g] <- 1
    for(itr in 1:40)
    {
      mm1 <- (n1[g]-gamma_ig[1,g])/(nn1[g]-gamma_ig[1,g])
      mm0 <- (n0[g]+gamma_ig[1,g])/(nn0[g]+gamma_ig[1,g])
      gamma_ig[1,g] <- n1[g]-sum(yhat2.mle[which(dat2.mle$G==g)]*mm1/(mm0*(1-yhat2.mle[which(dat2.mle$G==g)])+yhat2.mle[which(dat2.mle$G==g)]*mm1))
      itr=itr+1
    }
  }
  gamma_ig[2,] <- -gamma_ig[1,]
  
  Q_ig <- matrix(NA, 2, numG)
  for(g in 1:numG)
  {
    Q_ig[1,g] <- (nn1[g]-gamma_ig[1,g])/(nn0[g]+nn1[g])
    Q_ig[2,g] <- (nn0[g]-gamma_ig[2,g])/(nn0[g]+nn1[g])
  }
  dat11.mle <- data.frame(count(dat, vars=c("y", "G")))
  
  strata_mat <- matrix(NA, 4, numG)
  for (g in 1:numG)
  {strata_mat[1,g] <- sum(1*(fdat.mle$G==g)*(fdat.mle$cases==1))
  strata_mat[2,g] <- sum(1*(fdat.mle$G==g)*(fdat.mle$conts==1))
  strata_mat[3,g] <- sum(1*(fdat.mle$G==g)*(fdat.mle$Death-fdat.mle$cases==1))
  strata_mat[4,g] <- sum(1*(fdat.mle$G==g)*(fdat.mle$nonDeath-fdat.mle$conts==1))}
  
  u_ig <- matrix(NA, 2, numG)
  for (i in 1:2){
    for (g in 1:numG){
      u_ig[i, g] <- 1 - strata_mat[i+2, g]/(sum(strata_mat[,g])*Q_ig[i, g])
      
    }
  }
  
  # All things that we need to calculate based on above information obtained 
  W_g <- numeric(numG)
  B_g <- matrix(NA, numbeta, numG)
  A_g <- matrix(NA, 2, numG)
  p_g <- numeric(numG) #p(G=g)
  for (g in 1:numG){
    #counts <- dat2.mle[dat2.mle$G==g, "counts"]
    counts <- dat2.mle[dat2.mle$G==g, "cases"] + dat2.mle[dat2.mle$G==g, "conts"]
    x1 <- dat2.mle[dat2.mle$G==g, "x1"]
    x2 <- dat2.mle[dat2.mle$G==g, "x2"]
    x3 <- dat2.mle[dat2.mle$G==g, "x3"]
    z <- dat2.mle[dat2.mle$G==g, "z"]
    X <- cbind(rep(1,length(x1)),x1,x2,x3,z) 
    yfit.mle <- yhat2.mle[dat2.mle$G==g]
    u <- c(u_ig[1, g], u_ig[2, g])
    ystar.mle <- u[1]*yfit.mle/(u[1]*yfit.mle+u[2]*(1-yfit.mle))
    
    #W_g
    W_g[g] <- sum(counts*(ystar.mle*(1-ystar.mle)))
    #B_g
    B_g[, g] <- t(X) %*% (counts*ystar.mle*(1-ystar.mle))
    #A_g
    A_g[, g] <- c(1/(strata_mat[1,g]-gamma_ig[1,g])-1/(strata_mat[1,g]+strata_mat[3,g]-gamma_ig[1,g]),
                  1/(strata_mat[2,g]-gamma_ig[2,g])-1/(strata_mat[2,g]+strata_mat[4,g]-gamma_ig[2,g]))    
    #p_g
    p_g[g] <- sum(strata_mat[, g])/sum(strata_mat)
  }
  
  A_0g <- apply(A_g, 2, sum)
  K_g <- 1/A_0g-W_g
  dgamma_1g <- matrix(NA, numbeta, numG)
  for (g in 1:numG){
    dgamma_1g[,g] <- -B_g[,g]/(1-A_0g[g]*W_g[g])
  }
  
  #covariate distribution NPMLE estimates for each strata: f(x1,x2,x3,z,g)=f(x1,x2,x3,z|g)*p(G=g)
  #since x1,x2 are continuous, each unique value will only be counted once
  #f(x1,x2,x3,z,g) is a function of beta
  f.mle <- numeric(nrow(dat2.mle))
  df_g <- matrix(NA, numbeta, nrow(dat2.mle)) #df(x1,x2,x3,z|g)/dbeta 
  df <- matrix(NA, numbeta, nrow(dat2.mle)) #df(x1,x2,x3,z)/dbeta
  s_2 <- matrix(NA, numbeta, nrow(dat2.mle)) #score for R_k=1
  for (k in 1:nrow(dat2.mle)){
    #y <- dat2.mle[k, "y"]
    y <- dat2.mle[k, "cases"]
    x1 <- dat2.mle[k, "x1"]
    x2 <- dat2.mle[k, "x2"]
    x3 <- dat2.mle[k, "x3"]
    z <- dat2.mle[k, "z"]
    #counts <- dat2.mle[k, "counts"]
    counts <- dat2.mle[k, "cases"] + dat2.mle[k, "conts"]
    g <- dat2.mle[k,"G"]
    u <- c(u_ig[1, g], u_ig[2,g])
    yfit.mle <- yhat2.mle[k]
    dgamma <- dgamma_1g[, g]
    a <- A_g[,g]
    f.mle[k] <- counts/(u[1]*yfit.mle+u[2]*(1-yfit.mle))/n
    df_g[,k] <- -counts/(sum(strata_mat[,g])*(u[1]*yfit.mle+u[2]*(1-yfit.mle))^2)*(c(1,x1,x2,x3,z)*yfit.mle*(1-yfit.mle)*(u[1]-u[2])+dgamma*((1-yfit.mle)*u[2]*a[2]-yfit.mle*u[1]*a[1]))
    df[, k] <- df_g[,k]*p_g[g]
    s_2[, k] <- c(1,x1,x2,x3,z)*(y-yfit.mle) + df_g[, k]/f.mle[k]*p_g[g]
  }
  
  s_1 <- matrix(NA, numbeta, nrow(dat11.mle)) #score for R_k=0
  for (k in 1:nrow(dat11.mle)){
    y <- dat11.mle[k, "y"]
    g <- dat11.mle[k, "G"]
    Q <- NULL
    if (y==1){
      Q = Q_ig[1, g]
      s_1[, k] <- -dgamma_1g[, g]/sum(strata_mat[,g])/Q
    }
    else{
      Q = Q_ig[2, g]
      s_1[,k] <- dgamma_1g[, g]/sum(strata_mat[,g])/Q
    }
  }
  #S_1 <- apply(rbind(s_1, c(strata_mat[4, ], strata_mat[3, ])), 2, function(x) x[1:2]*x[3])
  #S_2 <- apply(rbind(s_2, dat2.mle$cases+dat2.mle$conts), 2, function(x) x[1:2]*x[3])
  #apply(S_2, 1, sum) + apply(S_1, 1, sum)
  
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
  
  
  #calculate auc and SD
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
  AUC3_values[ci] <- auc(fpr.mle, tpr.mle)
  
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
  Asy_sd_AUC3[ci]<-sqrt(var_auc)
  
}

Result.AUC[,1]<-rou
Result.AUC[,2]<-beta0
Result.AUC[,3]<-mean(AUC3_values)
Result.AUC[,4]<-sd(AUC3_values)
Result.AUC[,5]<-mean(Asy_sd_AUC3)
Result.AUC<-data.frame(Result.AUC)
print('okk')