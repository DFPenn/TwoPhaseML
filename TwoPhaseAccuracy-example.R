#############################################################
#############################################################
## The file example.R includes one example of how to estimate
## TPR, FPR & AUC as well as calculating standard errors with
## bootstrapping under Two-phase study designs 
#############################################################
#############################################################

library(osDesign)
library(mvtnorm)
library(plyr)
library(MESS)
library(boot)
library(progress)
source("/Users/yaqicao/Desktop/Postdoc in Penn/Oct,2021/10.7-RPackage/AUC-10.8/10.8-functions-Yaqi.R")



################################
### Simulate a two-phase dataset for illustration purpose
################################
set.seed(2016311033)
Prev <- "Moder" #"Moder"/Rare"
Rho <- "0.5" #"Independent"/"0.3"/"0.5"
ratio <- 2 #1/2
n <- 3000
if(Prev == "Rare")
{beta <- log(c(0.03,0.6,1.6,0.6,1.5))}
if(Prev == "Moder")
{beta <- log(c(0.13,0.6,1.6,0.6,1.5))}
numbeta <- length(beta)
if (Rho == "Independent")
{x1 <- rnorm(n)
z <- rnorm(n)}
if (Rho == "0.3")
{sigma <- matrix(c(1, 0.3, 0.3, 1), nrow=2, byrow=TRUE) #0.3*sqrt(1)*sqrt(1)
x1z <- rmvnorm(n, mean=c(0, 0), sigma=sigma)
x1 <- x1z[,1]
z <- x1z[,2]
#cor(x1,z)
}
if (Rho == "0.5")
{sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2, byrow=TRUE) 
x1z <- rmvnorm(n, mean=c(0, 0), sigma=sigma)
x1 <- x1z[,1]
z <- x1z[,2]
#cor(x1,z)
}
x2 <- runif(n)
x3 <- rbinom(n, size=1, prob=0.2)
mat <- cbind(rep(1,n), x1, x2, x3, z)
Outcome <- rbinom(n, size=1, prob=exp(mat %*% beta)/(1+exp(mat %*% beta)))
ddat <- data.frame(Outcome,x1,x2,x3,z)
glm.fit <- glm(Outcome~1+x1+x2+x3,data=ddat,family = binomial(link = "logit"))
theta <- glm.fit$coefficients

prob <- exp(mat[,-5] %*% theta)/(1+exp(mat[,-5] %*% theta))
Quan <- unname(quantile(prob[ddat$Outcome==1],c(0.25,0.5,0.75)))

Stratum <- sapply(prob, function(x) if (x<=Quan[1]) x<-1 else if (x>Quan[1] & x<=Quan[2]) x<-2 else if (x>Quan[2] & x<=Quan[3]) x<-3
            else x<-4)
dat <-data.frame(cbind(Outcome, x1, x2, x3, z, Stratum))
dat$R <- ifelse(dat$Outcome==1, 1, 0)
x1 <- xtabs(~dat$Outcome+dat$Stratum)[2,]*ratio
for (g in 1:length(unique(Stratum))){
  dat_cont_g <- dat[dat$Outcome==0 & dat$Stratum==g, ]
  id <- sample(c(1:nrow(dat_cont_g)), size=x1[g], replace=F)
  dat[dat$Outcome==0 & dat$Stratum==g, ][id, "R"] <- 1
}




################################
### Estimation of predictive accuracy measures with bootstrap standard error estimates
################################
Outcome <- dat$Outcome
x1 <- dat$x1
x2 <- dat$x2
x3 <- dat$x3
z <- dat$z
X <- cbind(x1,x2,x3,z)
#z[dat$R==0] <- NA
#Z <- z
Z <- rep(1,n)
Stratum <- dat$Stratum
Phase_ID <- dat$R + 1
namesX <- c("x1","x2","x3","x4")
namesZ <- c("z")
#Stratum <- rep(1, length(Outcome))

#Benchmark
est_benchmark <- evalTwoPhase(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q=0.2, p=0.2, method="Benchmark")
se_benchmark <- seTwoPhase(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q=0.2, p=0.2, numBoot=30, method="Benchmark")
summaryTwoPhase(est_benchmark, se_benchmark, q=0.2, p=0.2, method="Benchmark")


#####################################################

Outcome <- dat$Outcome
x1 <- dat$x1
x2 <- dat$x2
x3 <- dat$x3
z <- dat$z
X <- cbind(x1,x2,x3)
z[dat$R==0] <- NA
Z <- z
Stratum <- dat$Stratum
Phase_ID <- dat$R + 1
namesX <- c("x1","x2","x3")
namesZ <- c("z")

#ML
est_ml <- evalTwoPhase(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q=0.2, p=0.2, method="ML")

se_ml <- seTwoPhase(Outcome, X, Z, Stratum, Phase_ID, namesX, namesZ, q=0.2, p=0.2, numBoot=30, method="ML")
summaryTwoPhase(est_ml, se_ml, q=0.2, p=0.2, method="ML")


