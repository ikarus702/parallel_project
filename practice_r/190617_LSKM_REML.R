library(parallel)

## read in the command line arguments
## run with: R CMD BATCH '--args p.val <- 10 q.val <- 30 off.val1 <- 0 off.val2 <- 0 Data.num <- c(1:50) NCORES=specified number of cores' RNAME.R
args <- commandArgs(TRUE)
if(length(args) > 0) 
  for(i in 1:length(args)) 
    eval(parse(text=args[[i]]))

# NCORES <- 2
cl <- makeCluster(NCORES)


n.val <- 100
xi.val1 <- 1
xi.val2 <- 1


# p=p.val
# q=q.val
# n=n.val
# xi1=xi.val1
# xi2=xi.val2
# off1=off.val1
# off2=off.val2


parallel.NGK <- function(Data.num,p=p.val,q=q.val,n=n.val,xi1=xi.val1,xi2=xi.val2,off1=off.val1,off2=off.val2){
  
  library(psych,lib.loc="~/R/lib" )
  library(mvtnorm,lib.loc="~/R/lib")
  library(matrixcalc,lib.loc="~/R/lib")  #,lib.loc="~/R/lib"
  
  ############# Setting with Xi matrix ##############
  Tot.p <- p + q 
  
  true.xi <- xi1
  true.sig.eps <- 1
  true.tau <- 10
  
  mat.I <- diag(n)
  
  
  load("NGK_PK_GP_Data.RData")
  
  
  ############### Apply NGK with Gauss ################
  
  result.summary <- list()
  est.summary <- list()
  est.prec.summary <- list()
  
  ######## Import dataset from the list
  
  y <- matrix(data.list[[Data.num]][,1],nrow=n)
  X <- matrix(0,nrow=n,ncol=1)
  z <- t(data.list[[Data.num]][,-c(1)])
  eps <- eps.list[[Data.num]]
  true.f <- matrix((unlist(y - eps)),nrow=n)
  
  scale.y <- y-mean(y)
  X <- X -mean(X)
  scale.f <- true.f - mean(true.f)
  
  #######################################
  
  diag.I <- diag(n)

  square.z <- list()
  off.square.z <- list()
  Off.index <- 0
  
  p <- Tot.p
  
  off.xi <- rep(0,c(p*(p-1)/2))
  
  
  for(j in 1:p){
    vec.z <- matrix(scale(z[j,]),ncol=1)
    square.z[[j]] <- vec.z%*%t(vec.z)
  }
  
  for(l in 2:p){
    for(m in l-c((l-1):1)){
      
      Off.index <- Off.index + 1
      off.square.z[[Off.index]] <- matrix(z[m,],ncol=1)%*%matrix(z[l,],nrow=1)
      
    }
  }
  
  
  ### step 1
  
  
  #### initial step for tau
  
  next.sign <- 0
  psi.index <- 0
  
  old.sig.eps <- 1
  opt.sig.eps <- 0
  opt.rho <- 0
  
  g.Iter <- 2
  
  
  tau.list <- exp(seq(0,3,length.out = 150))
  
  
  g.psi <- matrix(NA,2,1)
  H.psi <- diag(NA,2,2)
  
  Iter <- 100
  
  
  for(g.iter in 1:g.Iter){
    
    l <- 0 
    stop.sign <- 0
    
    save.psi <- matrix(999,length(tau.list),3)
    save.loglik <- rep(NA,length(tau.list))
    save.r2 <- rep(NA,length(tau.list))
    save.RSS <- rep(NA,length(tau.list))
    list.psi <- list()
    
    
    for(old.tau in tau.list){
      
        
        if(g.iter ==1){
          
       #  old.beta <- 1
          old.rho <- 1
          old.sig.eps <- 1
          
        }else{
          
          old.rho <- tmp.rho
          old.sig.eps <- tmp.sig.eps
          
        }

      l <- l+1
      loglik <- rep(NA,Iter)
      
      tmp.psi <- matrix(NA,Iter,2)
      
      for(iter in 1:Iter){
        
        for(j in 1:p){
          if(j == 1){
            est.poly.K <- (old.rho^(-1))*square.z[[j]]
          }else{
            est.poly.K <- est.poly.K+(old.rho^(-1))*square.z[[j]]
          }
        }
        
        if(all(off.xi!=0)){
          
          for(j in 1:Off.index){
            est.poly.K <- est.poly.K+2*(off.xi[j])*off.square.z[[j]]
            
          }
        }
        est.poly.K <- as.matrix(est.poly.K)
        
        ### Estimating tau
        old.psi <- matrix(c(old.rho,old.sig.eps),2,1)
        
        diag.I <- diag(n)
        
        V <- old.tau*est.poly.K + old.sig.eps*diag.I
        
        if(any(V==Inf)){
          break;
        }
        
        inv.V <- try(solve(V),silent=TRUE)
        
        if ('try-error' %in% class(inv.V)){
          new.psi <- matrix(999,2,1)
          break
        }
        
        #ginv.V <- solve(t(X)%*%inv.V%*%X)
        ginv.V <- 0
        
        P <- inv.V-inv.V%*%X%*%ginv.V%*%t(X)%*%inv.V
        
        V.prime.rho <- -old.tau*(Reduce('+',square.z))/((old.rho)^2)
        V.prime2.rho <- 2*old.tau*(Reduce('+',square.z))/((old.rho)^3) 
        V.prime.tr <- -(Reduce('+',square.z))/((old.rho)^2)
        
        #### score function #####
        
        
        g.psi[1,] <- g_rho <- as.numeric(-0.5*tr(V.prime.rho%*%P)+0.5*t(scale.y)%*%P%*%V.prime.rho%*%P%*%scale.y )
        g.psi[2,] <- g_eps <-  as.numeric(-0.5*tr(P)+0.5*t(scale.y)%*%P%*%P%*%scale.y)
        
        
        H.psi[1,1] <- as.numeric(0.5*tr(V.prime.rho%*%P%*%V.prime.rho%*%P-V.prime2.rho%*%P)
                                 -0.5*t(scale.y)%*%P%*%(2*V.prime.rho%*%P%*%V.prime.rho - V.prime2.rho)%*%P%*%scale.y)
        
        H.psi[1,2] <- as.numeric(0.5*tr(P%*%V.prime.rho%*%P)-0.5*t(scale.y)%*%P%*%(P%*%V.prime.rho+V.prime.rho%*%P)%*%P%*%scale.y)
        
        H.psi[2,1] <- H.psi[1,2]
        
        H.psi[2,2] <- as.numeric(0.5*tr(P%*%P)-t(scale.y)%*%P%*%P%*%P%*%scale.y)
        
        if(!is.positive.definite(H.psi)){
          H.psi <- H.psi + 10^(-6)*diag(2)
        }
        
        inv.H.psi <- try(solve(H.psi), silent=TRUE)
        if ('try-error' %in% class(inv.H.psi)){
          new.psi <- matrix(999,2,1)
          break
        }else{ 
          
          new.psi <- old.psi - solve(H.psi)%*%g.psi
          new.psi <- matrix(new.psi,2,1)
        }
        
        if(any(is.na(new.psi))){
          new.psi <- matrix(999,2,1)
          cat(l," has NA! \n")
          break;
        }
        
        inv.lambda <- old.tau/old.sig.eps
         A <- solve(diag.I + inv.lambda*est.poly.K)%*%(inv.lambda*est.poly.K)
        inv.tmp <- solve(diag.I+inv.lambda*est.poly.K)
        #A <- inv.tmp%*%(inv.lambda*est.poly.K+X%*%solve(t(X)%*%inv.tmp%*%X)%*%t(X)%*%inv.tmp)
        est.h <- matrix(A%*%scale.y,ncol=1)
        
        new.beta <- ginv.V%*%t(X)%*%inv.V%*%scale.y
        est.h <- matrix(A%*%scale.y-X%*%new.beta,ncol=1)
        
        #new.psi[2,] <- as.numeric(t(scale.y-est.h)%*%(scale.y-est.h)/(n-tr(A)))
        new.psi[2,] <- as.numeric(t(scale.y-A%*%scale.y)%*%(scale.y-A%*%scale.y)/(n-tr(A)))
        
        #loglik[iter] <- -0.5*(log(det(V))+t(scale.y)%*%P%*%scale.y)
     #   if(det(V) < 0){
          
    #      loglik[iter] <- -0.5*(log(det(ginv.V))+t(scale.y)%*%P%*%scale.y)
     #   }else if(det(ginv.V)<0){ 
    #      loglik[iter] <- -0.5*(log(det(V))+t(scale.y)%*%P%*%scale.y)
          
    #    }else if(det(V) <0 & det(ginv.V)<0){
          
    #      loglik[iter] <- -0.5*(t(scale.y)%*%P%*%scale.y)
    #    }else{
          
    #      loglik[iter] <- -0.5*(log(det(V))+log(det(ginv.V))+t(scale.y)%*%P%*%scale.y)
          
    #    }
        
        
        if(max((new.psi-old.psi)^2) < 10^(-5)){
          
          opt.tau <- old.tau
          
          opt.rho <- new.psi[1,]
          opt.sig.eps <- new.psi[2,]
          if(opt.rho < 0){
            opt.rho = 999
            opt.sig.eps = 999
          }	
          
          if(opt.sig.eps <0){
            opt.sig.eps = 999
          }
          
          for(j in 1:p){
            if(j == 1){
              est.poly.K <- (opt.rho^(-1))*square.z[[j]]
            }else{
              est.poly.K <- est.poly.K+(opt.rho^(-1))*square.z[[j]]
            }
          }
          
          if(all(off.xi!=0)){
            
            for(j in 1:Off.index){
              est.poly.K <- est.poly.K+2*(off.xi[j])*off.square.z[[j]]
            }
          }
          est.poly.K <- as.matrix(est.poly.K)
          
          save.loglik[l] <- max(loglik[which(!is.na(loglik))])
          save.psi[l,] <-  matrix(c(old.tau,opt.rho,opt.sig.eps),1,3)
          
          V <- opt.sig.eps*diag.I+opt.tau*est.poly.K
          inv.V <- solve(V)
          
          est.h <-  opt.tau*est.poly.K%*%inv.V%*%(scale.y) 
          #est.h <-  opt.tau*est.poly.K%*%inv.V%*%(scale.y-X%*%new.beta) 
          
          lm.out <- lm(scale.y~est.h)
          summary.out <- summary(lm.out)
          save.r2[l] <- summary.out$r.squared
          
          #  ginv.V <- solve(t(X)%*%inv.V%*%X)
          
          #  beta.hat <- ginv.V%*%t(X)%*%inv.V%*%scale.y 
          
          cat(l," is successful! \n")
          
          break;
        }else{
          tmp.psi[iter,] <- matrix(new.psi,1,2)
          
          old.tau <- old.tau
          
          old.rho <- new.psi[1,]
          
          old.sig.eps <- new.psi[2,]
          
          cat("The number of iteration is ",iter,"\n")
          
        }
        
        list.psi[[l]] <- tmp.psi
      }
      
      if(iter==Iter){
        opt.rho = 999
        opt.sig.eps =999
      }
      
      
      if((abs(opt.sig.eps-true.sig.eps)<0.1)){
        
        psi.index <- l 
        stop.sign <- 1
        
        break;
      }else if(g.iter==g.Iter & (abs(save.psi[l,2]-true.xi^(-1))<0.2 | abs(save.psi[l,3]-true.sig.eps)<0.2)){
        
        psi.index <- l 
        stop.sign <- 1
        
        break;
        
      }else if(g.iter==g.Iter & l==100){
        
        l.index = 0 
        for(l in  1:100){
          if(all(save.psi[l,]!=999)){
            if(l.index[1] ==0){
              l.index <- l
            }else{
              l.index <- c(l.index,l)
            }
          }
        }   
        psi.index <- max(l.index)
        cat("Pick the last one. \n")        
        break;
        
      }
      
    }
    
    if(g.iter < g.Iter){
      if(stop.sign==1){
        cat("Stop the algorithm! \n")
        break;
      }else{
        
        conv.index <- which(save.psi[,3]!=999)
        sum.index <- sum(abs(diff(save.psi[conv.index,3]))<10^(-4))
        
        if(sum.index > 30){
          
          rho.diff <- abs(save.psi[conv.index,2]-true.xi^(-1))          
          rho.index <- which(rho.diff==min(rho.diff))
          start.pt <- conv.index[rho.index]+1
          end.pt <- start.pt+1
        }else{
          
          tmp.diff <- (save.psi[,3]-true.sig.eps)
          start.pt <- which(abs(tmp.diff)==min(abs(tmp.diff)))
          
          if(tmp.diff[start.pt] < 0){
            end.pt <- start.pt-1
          }else{
            end.pt <- start.pt+1
          }
          
          if(which(abs(tmp.diff)==min(abs(tmp.diff)))==1){
            start.pt <- 79; end.pt <- 80
          }
          
        }
        tau.list <- try(seq(tau.list[start.pt],tau.list[end.pt],length.out = 100),silent=TRUE)  ### error happens
        
        
        if ('try-error' %in% class(tau.list)){
          
          next.sign <- 1
          break;
          
        }else{
          tmp.rho <-save.psi[start.pt,2]     
          tmp.sig.eps <- save.psi[start.pt,3]
        }
        
      }
    }
    
  }
  
  if(next.sign==1){
    
    cat(Data.num,"th Data has an error! \n")
    result.summary <- NULL
    
    
  }else{
    
    ## Estimation of all parameters with REML
    
    if((psi.index)!=0){
      old.rho <- save.psi[psi.index,2]
     # old.beta <- new.beta
      old.tau <- save.psi[psi.index,1]
      old.sig.eps <- save.psi[psi.index,3]
    }else{
      old.rho <- opt.rho
      old.tau <- opt.tau
      old.sig.eps <- opt.sig.eps
    }
    
    
    if((opt.tau/opt.rho)< 100){
      
      for(j in 1:p){
        if(j == 1){
          est.poly.K <- (opt.rho^(-1))*square.z[[j]]
        }else{
          est.poly.K <- est.poly.K +(opt.rho^(-1))*square.z[[j]]
        }
      }
      if(all(off.xi!=0)){
        
        for(j in 1:Off.index){
          est.poly.K <- est.poly.K + 2*(off.xi[j])*off.square.z[[j]]
        }
      }    
      est.poly.K <- as.matrix(est.poly.K)
      
      V <- opt.sig.eps*diag.I+opt.tau*est.poly.K
      inv.V <- solve(V)
      
       est.h <-  opt.tau*est.poly.K%*%inv.V%*%(scale.y)
      #est.h <-  opt.tau*est.poly.K%*%inv.V%*%(scale.y-X%*%new.beta) 
      plot(scale.y,est.h)
      plot(scale.f,est.h)
      
      est.summary <- list(opt.rho = opt.rho, opt.tau=opt.tau, opt.sig.eps=opt.sig.eps)
      
      
    }else{	
      
      g.psi <- matrix(NA,3,1)
      H.psi <- diag(NA,3,3)
      
      Iter <- 100
      
      tmp.psi <- matrix(NA,Iter,3)
      
      for(iter in 1:Iter){
        
        for(j in 1:p){
          if(j == 1){
            est.poly.K <- old.rho^(-1)*square.z[[j]]
          }else{
            est.poly.K <- est.poly.K+old.rho^(-1)*square.z[[j]]
          }
        }
        if(all(off.xi!=0)){
          
          for(j in 1:Off.index){
            est.poly.K <- est.poly.K+2*(off.xi[j])*off.square.z[[j]]
            
          }
        }
        est.poly.K <- as.matrix(est.poly.K)
        
        ### Estimating tau
        old.psi <- matrix(c(old.rho,old.tau,old.sig.eps),3,1)
        
        diag.I <- diag(n)
        
        V <- old.tau*est.poly.K + old.sig.eps*diag.I
        inv.V <- solve(V)
         ginv.V <- 0
        #ginv.V <- solve(t(X)%*%inv.V%*%X)
        
        P <- inv.V-inv.V%*%X%*%ginv.V%*%t(X)%*%inv.V
        
        V.prime.rho <- -old.tau*(Reduce('+',square.z))/(old.rho)^2
        V.prime2.rho <- 2*old.tau*(Reduce('+',square.z))/(old.rho)^3
        V.prime.tr <- -(Reduce('+',square.z))/(old.rho)^2
        
        #### score function #####
        
        g.psi[1,] <- g_rho <- as.numeric(-0.5*tr(V.prime.rho%*%P)+0.5*t(scale.y)%*%P%*%V.prime.rho%*%P%*%scale.y )
        g.psi[2,] <- g_tau <- as.numeric(-0.5*tr(est.poly.K%*%P)+0.5*t(scale.y)%*%P%*%est.poly.K%*%P%*%scale.y)
        g.psi[3,] <- g_eps <-  as.numeric(-0.5*tr(P)+0.5*t(scale.y)%*%P%*%P%*%scale.y)
        
        
        H.psi[1,1] <- as.numeric(0.5*tr(V.prime.rho%*%P%*%V.prime.rho%*%P-V.prime2.rho%*%P)
                                 -0.5*t(scale.y)%*%P%*%(2*V.prime.rho%*%P%*%V.prime.rho - V.prime2.rho)%*%P%*%scale.y)
        
        H.psi[1,2] <- as.numeric(0.5*tr(est.poly.K%*%P%*%V.prime.rho%*%P-V.prime.tr%*%P)
                                 -0.5*t(scale.y)%*%P%*%(est.poly.K%*%P%*%V.prime.rho-V.prime.tr+V.prime.rho%*%P%*%est.poly.K)%*%P%*%scale.y)
        H.psi[2,1] <- H.psi[1,2]
        
        H.psi[1,3] <- as.numeric(0.5*tr(P%*%V.prime.rho%*%P)-0.5*t(scale.y)%*%P%*%(P%*%V.prime.rho+V.prime.rho%*%P)%*%P%*%scale.y)
        
        H.psi[3,1] <- H.psi[1,3]
        
        H.psi[2,2] <- as.numeric(0.5*tr(est.poly.K%*%P%*%est.poly.K%*%P)-t(scale.y)%*%P%*%est.poly.K%*%P%*%est.poly.K%*%P%*%scale.y )
        
        H.psi[2,3] <- as.numeric(0.5*tr(est.poly.K%*%P%*%P)-0.5*t(scale.y)%*%P%*%(est.poly.K%*%P+P%*%est.poly.K)%*%P%*%scale.y)
        
        H.psi[3,2] <- H.psi[2,3]
        
        H.psi[3,3] <- as.numeric(0.5*tr(P%*%P)-t(scale.y)%*%P%*%P%*%P%*%scale.y)
        
        if(!is.positive.definite(H.psi)){
          H.psi <- H.psi + 10^(-6)*diag(3)
        }
        
        
        new.psi <- matrix(old.psi - solve(H.psi)%*%g.psi,3,1)
        
        inv.lambda <- old.tau/old.sig.eps
        
         A <- solve(diag.I + inv.lambda*est.poly.K)%*%(inv.lambda*est.poly.K)
        # est.h <- matrix(A%*%scale.y,ncol=1)
        
        inv.tmp <- solve(diag.I+inv.lambda*est.poly.K)
        #A <- inv.tmp%*%(inv.lambda*est.poly.K+X%*%solve(t(X)%*%inv.tmp%*%X)%*%t(X)%*%inv.tmp)
        est.h <- matrix(A%*%scale.y,ncol=1)
        
        #new.beta <- ginv.V%*%t(X)%*%inv.V%*%scale.y
        #est.h <- matrix(A%*%scale.y-X%*%new.beta,ncol=1)
        
        new.psi[2,] <- as.numeric(t(scale.y-est.h)%*%(scale.y-est.h)/(n-tr(A)))
        #new.psi[2,] <- as.numeric(t(scale.y-A%*%scale.y)%*%(scale.y-A%*%scale.y)/(n-tr(A)))
        
        loglik[iter] <- -0.5*(log(det(V))+t(scale.y)%*%P%*%scale.y)
        #loglik[iter] <- -0.5*(log(det(V))+det(ginv.V)+t(scale.y)%*%P%*%scale.y)
        
        if(new.psi[2,] > 100){
          
          if((psi.index)!=0){
            opt.rho <- save.psi[psi.index,2]
           # opt.beta <- new.beta
            opt.tau <- save.psi[psi.index,1]
            opt.sig.eps <- save.psi[psi.index,3]
          }else{
            opt.rho <- opt.rho
            opt.tau <- opt.tau
            opt.sig.eps <- opt.sig.eps
          }
          cat("rho is too large! \n")
          break;
        }
        
        if(max((new.psi-old.psi)^2) < 10^(-5)){
          
          opt.rho <- new.psi[1,]
          opt.rho <- opt.rho*I(opt.rho>0)
          opt.tau <- new.psi[2,]
          
          
          opt.sig.eps <- new.psi[3,]
          
          for(j in 1:p){
            if(j == 1){
              est.poly.K <- opt.rho^(-1)*square.z[[j]]
            }else{
              est.poly.K <- est.poly.K+opt.rho^(-1)*square.z[[j]]
            }
          }
          if(all(off.xi!=0)){
            
            for(j in 1:Off.index){
              est.poly.K <- est.poly.K+2*(off.xi[j])*off.square.z[[j]]
              
            }
          }
          est.poly.K <- as.matrix(est.poly.K)
          
          V <- opt.sig.eps*diag.I+opt.tau*est.poly.K
          inv.V <- solve(V)
           ginv.V <- 0
          #ginv.V <- solve(t(X)%*%inv.V%*%X)
          
          #  beta.hat <- ginv.V%*%t(X)%*%inv.V%*%scale.y 
          
          cat("It is successful! \n")
          break;
        }else{
          tmp.psi[iter,] <- matrix(new.psi,1,3)
          old.tau <- new.psi[2,]
          old.rho <- new.psi[1,]
          
          old.sig.eps <- new.psi[3,]
          
          cat("The number of iteration is ",iter,"\n")
          
        }
        
      }
      V <- opt.sig.eps*diag.I+opt.tau*est.poly.K
      inv.V <- solve(V)
      est.h <-  opt.tau*est.poly.K%*%inv.V%*%(scale.y) 
     # est.h <-  opt.tau*est.poly.K%*%inv.V%*%(scale.y-X%*%new.beta) 
      par(mfrow=c(1,2))
      plot(scale.y,est.h)
      plot(scale.f,est.h)
      
      est.summary <- list(opt.rho = opt.rho, opt.tau=opt.tau, opt.sig.eps=opt.sig.eps)
    }
  }
  
  
  if(opt.tau/opt.rho > 100 | opt.rho==999){
    result.summary <- NULL
  } 
  ########################### NGK ########################################
  
  
  return(est.summary)
}


withGlobals <- function(FUN, ...){ 
  environment(FUN) <- list2env(list(...)) 
  FUN 
} 


result.summary <- parLapply(cl,Data.num,withGlobals(parallel.NGK,p.val=p.val,q.val=q.val,n.val=n.val,xi.val1=xi.val1,xi.val2=xi.val2,off.val1=off.val1,off.val2=off.val2))    

save.image(paste0("Semi_NGK_PK_PK_",Data.num[50],"_p40.RData"))

stopCluster(cl) 
