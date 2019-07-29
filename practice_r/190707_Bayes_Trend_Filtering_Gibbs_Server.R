
library(parallel)

## read in the command line arguments
## run with: R CMD BATCH '--args p.val <- 10 q.val <- 60 off.val1 <- 0 off.val2 <- 0 Data.num <- c(1:50) NCORES=specified number of cores' RNAME.R
args <- commandArgs(TRUE)
if(length(args) > 0) 
  for(i in 1:length(args)) 
    eval(parse(text=args[[i]]))

# NCORES <- 2
cl <- makeCluster(NCORES)

parallel.Bayes.BTF <- function(Data.num, d.num=d.val, m.steps=m.val){


library(statmod)
library(invgamma)
library(mvtnorm)
library(VGAM)

load("~/src/RData/Trend3_Data.RData")

y <- matrix(data.list[[Data.num]]$y,ncol=1)
X <- matrix(data.list[[Data.num]]$X,ncol=1) 

n <- length(X)
d.order <- d.val

gibbsBTF = function(x=diag(1,n), y, D.order, max.steps = m.steps) {
  n <- nrow(x)
  p <- n
  
  XtX <- t(x) %*% x	#Time saving
  xy <- t(x) %*% y
  
  r <- 1
  delta <- 0.1 # 1.78
  
  betaSamples <- matrix(0, max.steps, p)
  sigma2Samples <- rep(0, max.steps)
  #invTau2Samples <- matrix(0, max.steps, p)
  invW2Samples <- matrix(0, max.steps, p-D.order-1)
  lambdaSamples <- rep(0, max.steps)
  #lambda2Samples <- rep(0, max.steps)
  
  beta <- drop(backsolve(XtX + diag(nrow=p), xy))
  residue <- drop(y - x %*% beta)
  sigma2 <- drop((t(residue) %*% residue) / n)
  #invTau2 <- 1 / (beta * beta)
  
  D.opt <- matrix(0,n-1,n)
  for(i in 1:(n-1)){
    
    D.opt[i,c(i,i+1)] <- c(1,-1)
    
  }
  
  if(D.order==0){
    
    D.matrix <- D.opt
    
  }else if(D.order > 0 ){
    
    D.matrix <- D.opt
    
    for(d in 1:D.order)
    
    D.matrix <- D.opt[1:c(n-d-1),1:c(n-d)] %*% D.matrix  
  }
    
    
  diff.beta <- D.matrix %*% beta
  invW2 <- as.numeric(1/(diff.beta * diff.beta))
  
  #lambda1 <- sqrt(2*p / sum(beta^2))
  lambda <- sqrt(2*(p-D.order-1)/ sum(diff.beta^2))
  
  k <- 0
  while (k < max.steps) {
    k <- k + 1
    
    if (k %% 1000 == 0) {
      cat('Iteration:', k, "\r")
    }
    
    # sample beta
    invBeta <- diag(invW2,c(p-D.order-1))
    
    invBeta <- t(D.matrix)%*%invBeta%*%D.matrix
    
    invA <- solve(XtX + invBeta)
    mean <- invA %*% xy
    varcov <- sigma2 * invA
    beta <- drop(rmvnorm(1, mean, varcov))
    betaSamples[k,] <- beta
    
    # sample sigma2
    shape <- n
    residue <- drop(y - x %*% beta)
    scale <- (t(residue) %*% residue + t(beta) %*% invBeta %*% beta)/2
    sigma2 <- 1/rgamma(1, shape, 1/scale)
    sigma2Samples[k] <- sigma2
    

    # sample w2
    diff.beta <- diff(beta)
    
    muPrime <- sqrt(lambda^2 * sigma2 / diff.beta^2)
    lambdaPrime <- lambda^2
    invW2 <- rep(0, p-D.order-1)
    for (i in seq(p-D.order-1)) {
      invW2[i] <- rinv.gaussian(1, muPrime[i], lambdaPrime)
    }
    invW2Samples[k, ] <- invW2
    
    
  
    # update lambda
    shape = r + p -D.order-1
    scale = delta + sum(1/invW2)/2
    #lambda <- sqrt(2*(p-1)/sum(1/invW2))
    lambda <- sqrt(rgamma(1, shape, 1/scale))
    lambdaSamples[k] <- lambda
    
  }
  
  result.list <- list(beta=betaSamples,sigma2=sigma2Samples,lambda=lambdaSamples)
  
  return(result.list)
  
}

est.summary <- gibbsBTF(y=y, D.order=d.order)
betas <- est.summary$beta[m.steps,]
MSE <- sum((y-betas)^2)/(n-1)

summary.result <- list(beta=est.summary$beta, sigma2= est.summary$sigma2, lambda = est.summary$lambda, MSE=MSE)

return(summary.result)
}

withGlobals <- function(FUN, ...){ 
  environment(FUN) <- list2env(list(...)) 
  FUN 
} 
#Data.num <- c(1:2)

result.summary <- parLapply(cl,Data.num,withGlobals(parallel.Bayes.BTF,d.val=d.val, m.val=m.val))

save(result.summary, file="~/src/RData/Trend3_Result.RData")

stopCluster(cl) 

