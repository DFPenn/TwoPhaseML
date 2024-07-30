###############################################################
###############################################################
## The file functions.R includes the main functions used in the paper:
## "Two-phase stratified sampling and analysis for predicting binary outcomes".
##
## Key Notations: 
##     Outcome: outcome information (1: case; 0: control)
##     X: standard predictors whose information is available to all subjects
##     Z: new biomarkers whose information is only available to a subset of subjects
##     Phase_ID: information of whether a subject is Phase I or Phase II (1: phase I; 2: phase II)
##     Stratum: stratification variables defined based on X
##
## Methods: 
##     Benchmark: fit logistic regression model to all subjects: logitP(Y=1|X) ~ X
#Yaqi: Actually, here benchmark means using Phase I variables as covariates.
##     MLE: fit logistic regression model to Two-phase dataset with maximum likelihood: logitP(Y=1|X,Z) ~ X+Z
##   
## Predictive accuracy measures:
##     TPR(q), FPR(p), AUC
###############################################################
###############################################################


###################################
## Main functions of estimating TPR, FPR & AUC with bootstrapping method for standard
## error calculations under two-phase study design 
##  Functions:
##      evalTwoPhase(): estimate predictive accuracy measures
##      seTwoPhase(): estimate bootstrap standard errors
##      evalTwoPhase.fun(): helper function used in seTwoPhase()
##      summaryTwoPhase(): print results
##
##  Input arguments:
##      Outcome: vector of length n 
##      X: matrix of size n x t (t standard predictors included in the model)
##      Z: matrix of size n x m (m biomarkers included in the model); some subjects have NA values
##      Stratum: vector of length n
##      Phase_ID: vector of length n
##      namesX: names of X; vector of length t
##      namesZ: names of Z; vector of length m
##      method: "Benchmark", "ML"
##      q: threshold value for TPR(q)
##      p: threshold value for FPR(p)
##      numBoot: number of bootstrap samples
##################################
evalTwoPhase <- function(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q, p, method){
  n <- length(Outcome)
  dat_list <- Dat_format(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, method)
  fit_list <- Dat_fit(method, dat_list)
  tpr <- Dat_tpr(q, method, dat_list, fit_list)
  fpr <- Dat_fpr(p, method, dat_list, fit_list)
  auc <- Dat_auc(method, dat_list, fit_list)
  return(c(tpr, fpr, auc))
}

seTwoPhase <- function(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q, p, numBoot, method){
  
  n <- length(Outcome)
  data <- data.frame(cbind(Outcome, X, Z, Stratum, Phase_ID))
  names(data) <- c("Y",namesX,namesZ,"G","Phase_ID")
  if(method=="Benchmark")
  {evalTwoPhase.boot <- boot(data, evalTwoPhase.fun, R=numBoot, namesX=namesX, namesZ=namesZ, q=q, p=p, method=method)}
  if(method=="ML")
  {evalTwoPhase.boot <- boot(data, evalTwoPhase.fun, R=numBoot, strata=Stratum, namesX=namesX, namesZ=namesZ, q=q, p=p, method=method)}
  return (apply(evalTwoPhase.boot$t, 2, sd))
}


# helper functions 
evalTwoPhase.fun <- function(data, indices, namesX, namesZ, q, p, method){
  d <- data[indices, ]
  Outcome <- d$Y
  X <- d[, namesX]
  Z <- d[, namesZ] 
  Stratum <- d$G
  Phase_ID <- d$Phase_ID
  return (evalTwoPhase(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q, p, method))
}

summaryTwoPhase <- function(eval_vec, se_vec, q, p, method){
  cat(paste("Method: ", method), "\n",
      "Results: Est(se)", "\n",
      paste("TPR(", q, "): ",round(eval_vec[1],3), " (", round(se_vec[1],3), ")", sep=""), "\n",
      paste("FPR(", p, "): ",round(eval_vec[2],3), " (", round(se_vec[2],3), ")", sep=""), "\n",
      paste("AUC: ", round(eval_vec[3],3), " (", round(se_vec[3],3), ")", sep=""))
}








##################################
## Function of formatting dataset
##################################
Dat_format <- function(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, method){
  dat <- data.frame(cbind(Outcome, X, Z, Stratum))
  names(dat) <- c("y",namesX,namesZ,"G")
  dat2 <- dat[Phase_ID==2, ]
  
  if (method=="Benchmark"){
    return (list(dat1=dat[,c("y",namesX)]))
  }
  if (method=="ML"){
    name <- names(dat)
    v1 <- length(namesX); v2 <- length(namesZ)
    for (i in 1:length(namesZ))
    {dat[,namesZ[i]] <- ifelse(is.na(dat[,namesZ[i]]), 10, dat[,namesZ[i]])}#fake values
    l <- list(G=dat$G)
    for (i in 1:(v1+v2))
    {l[[name[i+1]]] <- dat[,(1+i)]}
    nonDeath <- aggregate(1-dat[,1], by=l, FUN=sum)$x
    
    
    dat1.mle <- data.frame(aggregate(dat[,1], by=l, FUN=sum))
    dat1.mle <- data.frame(cbind(dat1.mle, nonDeath))
    names(dat1.mle) <- c("G",name[2:(v1+v2+1)], "Death","nonDeath")
    
    l2 <- list(G=dat2$G)
    for (i in 1:(v1+v2))
    {l2[[name[i+1]]] <- dat2[,(1+i)]}
    conts <- aggregate(1-dat2[,1], by=l2, FUN=sum)$x
    dat2.mle <- data.frame(aggregate(dat2[,1], by=l2, FUN=sum))
    dat2.mle <- data.frame(cbind(dat2.mle, conts))
    names(dat2.mle) <- c("G",name[2:(v1+v2+1)],"cases","conts")
    
    # Final = Phase I + Phase II
    fdat.mle <- merge(dat1.mle, dat2.mle, by=c("G",name[2:(v1+v2+1)]), all=T)
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
      return(list(dat2=dat2,fdat.mle=fdat.mle, dat2.mle=dat2.mle, dat1.mle=dat1.mle, dat=dat, w1=w1, w0=w0, mu_1=mu_1, mu_0=mu_0, nn0=nn0, nn1=nn1, n0=n0, n1=n1,namesX=namesX,namesZ=namesZ))
    }
}




##################################
## Model Fitting:
##   Dat_fit(): wrapper function
##   Dat_benchmarkFit(): function to fit the logistic regression model to full cohort with only X
##   Dat_mleFit(): function to fit the logistic regression model to Two-phase dataset with maximum likelihood with both X and Z
##################################
Dat_fit <- function(method, dat_list){
  if (method=="Benchmark"){
    return(Dat_benchmarkFit(dat_list))
  }
 if (method=="ML"){
   return(Dat_mleFit(dat_list))
 }
}

Dat_benchmarkFit <- function(dat_list){
  #data: generated by Dat_format("Benchmark")
  mod <- glm(y ~ ., family=binomial, data=dat_list$dat1)
  yhat <- mod$fitted.values
  return (list(yhat=yhat))
}

Dat_mleFit <- function(dat_list){
  # dat_list: generated by Dat_format("ML")
  fdat.mle <- dat_list$fdat.mle
  dat2.mle <- dat_list$dat2.mle
  dat1.mle <- dat_list$dat1.mle
  dat <- dat_list$dat
  dat2 <- dat_list$dat2
  w1 <- dat_list$w1
  w0 <- dat_list$w0
  mu_1 <- dat_list$mu_1
  mu_0 <- dat_list$mu_0
  nn0 <- dat_list$nn0
  nn1 <- dat_list$nn1
  n0 <- dat_list$n0
  n1 <- dat_list$n1
  case <- fdat.mle$cases
  group <- fdat.mle$G
  namesX <- dat_list$namesX; namesZ <- dat_list$namesZ
  x <- cbind(rep(1,nrow(fdat.mle)),fdat.mle[,c(namesX,namesZ)])
  N <- fdat.mle$cases+fdat.mle$conts
  fdat.mle1 <- fdat.mle[,c(namesX,namesZ,"cases","conts")]
  mod <- tps(cbind(cases, conts) ~ ., data=fdat.mle1, nn0=nn0, nn1=nn1, group=fdat.mle$G, method="ML", cohort=T)
  beta.mle <- mod$coef
  
  XX2 <- cbind(rep(1,nrow(dat2.mle)), dat2.mle[,c(namesX,namesZ)])
  yhat2.mle <- exp(as.matrix(XX2) %*% beta.mle)/(1 + exp(as.matrix(XX2) %*% beta.mle))
  X2 <- cbind(rep(1,nrow(dat2)), dat2[,c(namesX,namesZ)])
  y2.mle <- exp(as.matrix(X2) %*% beta.mle)/(1 + exp(as.matrix(X2) %*% beta.mle))
  numG <- length(unique(fdat.mle$G))
  gamma_ig <- matrix(NA, 2, numG)
  for(g in 1:numG)
  {
    gamma_ig[1,g] <- 1
    for(itr in 1:40)
    {
      mm1 <- (n1[g]-gamma_ig[1,g])/(nn1[g]-gamma_ig[1,g])
      mm0 <- (n0[g]+gamma_ig[1,g])/(nn0[g]+gamma_ig[1,g])
      datg <- dat2.mle[which(dat2.mle$G==g),]
      gamma_ig[1,g] <- n1[g]-sum((datg[,"cases"]+datg[,"conts"])*yhat2.mle[which(dat2.mle$G==g)]*mm1/(mm0*(1-yhat2.mle[which(dat2.mle$G==g)])+yhat2.mle[which(dat2.mle$G==g)]*mm1))
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
 
  strata_mat[1,] <- xtabs(~dat2$y+dat2$G)[2,]
  strata_mat[2,] <- xtabs(~dat2$y+dat2$G)[1,]
  strata_mat[3,] <- xtabs(~dat$y+dat$G)[2,]-xtabs(~dat2$y+dat2$G)[2,]
  strata_mat[4,] <- xtabs(~dat$y+dat$G)[1,]-xtabs(~dat2$y+dat2$G)[1,]
  
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
    counts <- dat2.mle[dat2.mle$G==g, "cases"] + dat2.mle[dat2.mle$G==g, "conts"]
    X <- rep(1,length(counts))
    for(ii in (2:(length(namesX)+length(namesZ)+1)))
    {assign(names(dat2.mle)[ii],dat2.mle[dat2.mle$G==g, ii])
      X<- cbind (X, dat2.mle[dat2.mle$G==g, ii])}
    
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
  df_g2 <- matrix(NA, numbeta, nrow(dat2)) #df(x1,x2,x3,z|g)/dbeta 
  df <- matrix(NA, numbeta, nrow(dat2.mle)) #df(x1,x2,x3,z)/dbeta
  f2.mle <- numeric(nrow(dat2))
  df2 <- matrix(NA, numbeta, nrow(dat2)) #df(x1,x2,x3,z)/dbeta
  s_2 <- matrix(NA, numbeta, nrow(dat2)) #score for R_k=1
  for (k in 1:nrow(dat2.mle)){
    #y <- dat2.mle[k, "y"]
    y <- dat2.mle[k, "cases"]
    xx <- 1
    for(ii in (2:(length(namesX)+length(namesZ)+1)))
    {assign(names(dat2.mle)[ii],dat2.mle[k, ii])
      xx <- c(xx,dat2.mle[k, ii])}
    
    counts <- dat2.mle[k, "cases"] + dat2.mle[k, "conts"]
    g <- dat2.mle[k,"G"]
    u <- c(u_ig[1, g], u_ig[2,g])
    yfit.mle <- yhat2.mle[k]
    dgamma <- dgamma_1g[, g]
    a <- A_g[,g]
    
    #f.mle[k] <- counts/(sum(strata_mat[,g])*(u[1]*yfit.mle+u[2]*(1-yfit.mle)))*p_g[g]
    f.mle[k] <- counts/(u[1]*yfit.mle+u[2]*(1-yfit.mle))/n
    df_g[,k] <- -counts/(sum(strata_mat[,g])*(u[1]*yfit.mle+u[2]*(1-yfit.mle))^2)*(xx*yfit.mle*(1-yfit.mle)*(u[1]-u[2])+dgamma*((1-yfit.mle)*u[2]*a[2]-yfit.mle*u[1]*a[1]))
    df[, k] <- df_g[,k]*p_g[g]
    
  }
  for (k in 1:nrow(dat2)){
    y <- dat2[k, "y"]
    xx2 <- 1
    for(ii in (2:(length(namesX)+length(namesZ)+1)))
    {assign(names(dat2)[ii],dat2[k, ii])
      xx2 <- c(xx2,dat2[k, ii])}
    
    g <- dat2[k,"G"]
    u <- c(u_ig[1, g], u_ig[2,g])
    yfit.mle <- y2.mle[k]
    dgamma <- dgamma_1g[, g]
    a <- A_g[,g]
    
    f2.mle[k] <- 1/(u[1]*yfit.mle+u[2]*(1-yfit.mle))/n
    df_g2[,k] <- -1/(sum(strata_mat[,g])*(u[1]*yfit.mle+u[2]*(1-yfit.mle))^2)*(xx2*yfit.mle*(1-yfit.mle)*(u[1]-u[2])+dgamma*((1-yfit.mle)*u[2]*a[2]-yfit.mle*u[1]*a[1]))
    df2[, k] <- df_g2[,k]*p_g[g]
    s_2[, k] <- xx2*(y-yfit.mle) + df_g2[, k]/f2.mle[k]*p_g[g]
  }
  s_1 <- matrix(NA, numbeta, nrow(dat11.mle)) #score for R_k=0
  for (k in 1:nrow(dat11.mle)){
    y <- dat11.mle[k, "y"]
    #y <- dat1.mle[k, "Death"]
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
  #h_1 <- out$cove %*% s_1*n
  #h_2 <- out$cove %*% s_2*n
  emp_cov <- matrix(0, numbeta, numbeta)
  counts2 <- c(strata_mat[4,], strata_mat[3,])
  for (k in 1:ncol(h_1)){
    emp_cov = emp_cov + h_1[, k] %*% t(h_1[, k])*counts2[k]
    #emp_cov = emp_cov + h_1[, k] %*% t(h_1[, k])*(1-fdat.mle$cases[k]-fdat.mle$conts[k])
  }
  for (k in 1:ncol(h_2)){
    #emp_cov = emp_cov + h_2[, k] %*% t(h_2[, k])*dat2.mle$counts[k]
    emp_cov = emp_cov + h_2[, k] %*% t(h_2[, k])
  }
  # check emp_cov/(n^2) = mod$cove
  
  return (list(beta.mle=beta.mle, beta.var.mle=diag(mod$cove),h1=h_1, h2=h_2, y2.mle=y2.mle, yhat.mle=yhat2.mle, f.mle=f.mle,f2.mle=f2.mle, df_beta=df,df2=df2, strata_mat=strata_mat))
}

##################################
## Calculate TPR, FPR & AUC:
##   Dat_tpr(), Dat_fpr(), Dat_auc(): wrapper functions
##   tpr_benchmark(), fpr_benchmark(), auc_benchmark(): function to calcualte pcf,pnf,auc using full cohort
##   tpr_mle(), fpr_mle(), auc_mle(): function to calculate pcf,pnf,auc uisng two-phase data with maximum likelihood
##################################
Dat_tpr <- function(q, method, dat_list, fit_list){
  if (method=="Benchmark"){
    return(tpr_benchmark(q, fit_list))
  }
  if (method=="ML"){
    return(tpr_mle(q, dat_list, fit_list))
  }
}

Dat_fpr <- function(p, method, dat_list, fit_list){
  if (method=="Benchmark"){
    return(fpr_benchmark(p, fit_list))
  }
  if (method=="ML"){
    return(fpr_mle(p, dat_list, fit_list))
  }
}

Dat_auc <- function(method, dat_list, fit_list){
  if (method=="Benchmark"){
    return(auc_benchmark(fit_list))
  }
  if (method=="ML"){
    return(auc_mle(dat_list, fit_list))
  }
}

tpr_benchmark <- function(q, fit_list){
  yhat.full <- fit_list$yhat
  c <- q
  denom_tpr.full <- sum(yhat.full)
    ind <- (yhat.full>=q)
    num_tpr.full <- sum(yhat.full[ind==1])
  tpr.full <- num_tpr.full/denom_tpr.full
  return(tpr.full)
}

fpr_benchmark <- function(p,fit_list){
  yhat.full <- fit_list$yhat
  c <- p
  denom_fpr.full <- sum(1-yhat.full)
    ind <- (yhat.full>=p)
    num_fpr.full <- sum((1-yhat.full)[ind==1])
  
  fpr.full <- num_fpr.full/denom_fpr.full
  
  return(fpr.full)
}


auc_benchmark <- function(fit_list){
  yhat.full <- fit_list$yhat
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
  auc <- auc(fpr.full, tpr.full)
  
  return(auc)
}

tpr_mle <- function(q, dat_list,fit_list){
  fdat.mle <- dat_list$fdat.mle
  dat1.mle <- dat_list$dat1.mle
  dat2.mle <- dat_list$dat2.mle
  dat2 <- dat_list$dat2
  
  #fit_list: generated by Dat_mleFit(data_list)
  fhat.mle <- fit_list$f.mle
  f2.mle <- fit_list$f2.mle
  yhat.mle <- fit_list$yhat.mle
  y2.mle <- fit_list$y2.mle
  beta.mle <- fit_list$beta.mle
  denom_tpr.mle <- sum(y2.mle*f2.mle)
  c <- q
    ind <- (y2.mle>=c)
    num_tpr.mle <- sum((y2.mle*f2.mle)[ind==1])
  tpr.mle <- num_tpr.mle/denom_tpr.mle
  return(tpr.mle)
}

fpr_mle <- function(p, dat_list,fit_list){
  fdat.mle <- dat_list$fdat.mle
  dat1.mle <- dat_list$dat1.mle
  dat2.mle <- dat_list$dat2.mle
  dat2 <- dat_list$dat2
  
  #fit_list: generated by Dat_mleFit(data_list)
  fhat.mle <- fit_list$f.mle
  f2.mle <- fit_list$f2.mle
  yhat.mle <- fit_list$yhat.mle
  y2.mle <- fit_list$y2.mle
  beta.mle <- fit_list$beta.mle
  
  denom_fpr.mle <- sum((1-y2.mle)*f2.mle)
  c <- p
    ind <- (y2.mle>=c)
    num_fpr.mle <- sum(((1-y2.mle)*f2.mle)[ind==1])
  fpr.mle <- num_fpr.mle/denom_fpr.mle
  return(fpr.mle)
}

auc_mle <- function(dat_list,fit_list){
  fdat.mle <- dat_list$fdat.mle
  dat1.mle <- dat_list$dat1.mle
  dat2.mle <- dat_list$dat2.mle
  dat2 <- dat_list$dat2
  
  #fit_list: generated by Dat_mleFit(data_list)
  fhat.mle <- fit_list$f.mle
  f2.mle <- fit_list$f2.mle
  yhat.mle <- fit_list$yhat.mle
  y2.mle <- fit_list$y2.mle
  beta.mle <- fit_list$beta.mle
 
  # calculate auc
  denom_tpr.mle <- sum(y2.mle*f2.mle)
  denom_fpr.mle <- sum((1-y2.mle)*f2.mle)
  c <- c(seq(0, 1, by=0.01)*max(yhat.mle),1)
  num_tpr.mle <- numeric(length(c))
  num_fpr.mle <- numeric(length(c))
  for (j in 1:length(c)){
    ind <- (y2.mle>=c[j])
    num_tpr.mle[j] <- sum((y2.mle*f2.mle)[ind==1])
    num_fpr.mle[j] <- sum(((1-y2.mle)*f2.mle)[ind==1])
  }
  tpr.mle <- num_tpr.mle/denom_tpr.mle
  fpr.mle <- num_fpr.mle/denom_fpr.mle
  auc <- auc(fpr.mle, tpr.mle)
  return(auc)
}





