library(pbdDMAT, quiet=TRUE)
library(psych)
library(matrixcalc)


## read in the command line arguments
args = commandArgs(trailingOnly=TRUE)

Nnode = as.numeric(args[1])


init.grid()

bldim <- c(16,16)

# p=p.val
# q=q.val
# n=n.val
# xi1=xi.val1
# xi2=xi.val2
# off1=off.val1
# off2=off.val2

load("~/src/RData/Unif_Semi_NGK_PK_Data_p200.RData")

Data.num <- 1

  ############# Setting with Xi matrix ##############
  
  true.xi <- 1
  true.sig.eps <- 1
  true.tau <- 10
  
  
if(comm.rank()==0){  
  ############### Apply NGK with Gauss ################
  start_time <- proc.time()  
  result.summary <- list()
  est.summary <- list()
  est.prec.summary <- list()
  
  ######## Import dataset from the list
  
  n <- nrow(data.list[[Data.num]])
  y <- matrix(data.list[[Data.num]][,1],nrow=n)
  X <- matrix(data.list[[Data.num]][,2],nrow=n)
  z <- t(data.list[[Data.num]][,-c(1,2)])
  eps <- eps.list[[Data.num]]
  true.f <- matrix((unlist(y - eps)),nrow=n)
  Tot.p <- nrow(z)
  

  scale.y <- y-mean(y)
  X <- X -mean(X)
  scale.f <- true.f - mean(true.f)
  
  #######################################
  
  diag.I <- diag(1,n)

  p <- nrow(z)
  
  off.xi <- rep(0,c(p*(p-1)/2))
}else{
y <- NULL
X <- NULL
z <- NULL
diag.I <- NULL
scale.y <- NULL
}  
  
  
#  for(l in 2:p){
#    for(m in l-c((l-1):1)){
      
#      Off.index <- Off.index + 1
#      off.square.z[[Off.index]] <- matrix(z[m,],ncol=1)%*%matrix(z[l,],nrow=1)
      
#    }
#  }
  
  
dz <- as.ddmatrix(x=z, bldim=bldim)
dy <- as.ddmatrix(x=scale.y, bldim=bldim)
dX <- as.ddmatrix(x=X, bldim=bldim)
dI <- as.ddmatrix(x=diag.I, bldim=bldim)

d_D_sum <- crossprod(dz,dz)

  ### step 1
  
  
  #### initial step for tau
  
  scale.y <- as.matrix(dy)

  next.sign <- 0
  psi.index <- 0
  
  old.sig.eps <- 1
  opt.sig.eps <- 0
  opt.rho <- 0
  
  g.Iter <- 2
  
  
  tau.list <- exp(seq(-2,3,length.out = 150))
  
  
  g_psi <- matrix(NA,2,1)
  H_psi <- matrix(NA,2,2)
  
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
          
          old.beta <- 1
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
        
            d_poly_K <- (old.rho^(-1))*d_D_sum
        
       # if(all(off.xi!=0)){
          
       #   for(j in 1:Off.index){
       #     est.poly.K <- est.poly.K+2*(off.xi[j])*off.square.z[[j]]
            
       #   }
       # }
       # est.poly.K <- as.matrix(est.poly.K)
        
        ### Estimating tau
        old.psi <- matrix(c(old.rho,old.sig.eps),2,1)
        
        d_V <- old.tau*d_poly_K + old.sig.eps*dI
        V <- as.matrix(d_V)
       
  
    #    if(any(V==Inf)){
    #      break;
    #    }
        
    #    inv.V <- try(solve(V),silent=TRUE)
        
    #    if ('try-error' %in% class(inv.V)){
    #      new.psi <- matrix(999,2,1)
    #      break
    #    }
        
        
        d_inv_V <- solve(d_V)  
        d_inv_XVX <- solve(t(dX) %*% d_inv_V %*% dX)
        #ginv.V <- 0
        inv_XVX <- as.matrix(d_inv_XVX)        

        d_P <- d_inv_V - d_inv_V %*% dX %*% d_inv_XVX %*% t(dX) %*% d_inv_V
        
        d_dV_rho <- -(old.tau/(old.rho^2))*d_D_sum
        d_dV2_rho <- (2*old.tau/(old.rho^3))*d_D_sum 
        d_dV_tau <- -d_D_sum/((old.rho)^2)
        
        #### score function #####
      

	d_g_psi <- (-0.5*sum(diag(d_dV_rho %*% d_P))+0.5*t(dy) %*% d_P %*% d_dV_rho %*% d_P %*% dy)
	g_psi[1,] <- as.matrix(d_g_psi)      
        d_g_psi <- (-0.5*(sum(diag(d_P)))+0.5*t(dy) %*% d_P %*% d_P %*% dy)

        g_psi[2,] <- as.matrix(d_g_psi)        
        print(g_psi)
        
        d_H_psi <- (0.5*(sum(diag(d_dV_rho %*% d_P %*% d_dV_rho %*% d_P - d_dV2_rho %*% d_P)))
                                 -0.5*t(dy) %*% d_P %*% (2*d_dV_rho %*% d_P %*% d_dV_rho - d_dV2_rho) %*% d_P %*% dy)
        H_psi[1,1] <- as.matrix(d_H_psi) 

        d_H_psi <- (0.5*(sum(diag(d_P %*% d_dV_rho %*% d_P)))-0.5*t(dy) %*% d_P %*% (d_P %*% d_dV_rho+ d_dV_rho %*% d_P) %*% d_P%*% dy)
        
	H_psi[1,2] <- as.matrix(d_H_psi)
        H_psi[2,1] <- H_psi[1,2]
        
        d_H_psi <- (0.5* sum(diag(d_P %*% d_P))-t(dy) %*% d_P %*% d_P %*% d_P %*% dy)
	H_psi[2,2] <- as.matrix(d_H_psi)        
	
	print(H_psi)

        if(!is.positive.definite(H_psi)){
          H_psi <- H_psi + 10^(-6)*diag(1,2)
        }
        
        inv.H.psi <- try(solve(H_psi), silent=TRUE)
        if ('try-error' %in% class(inv.H.psi)){
          new.psi <- matrix(999,2,1)
          break
        }else{ 
          
          new.psi <- old.psi - inv.H.psi %*% g_psi
          new.psi <- matrix(new.psi,2,1)
        }
        
        if(any(is.na(new.psi))){
          new.psi <- matrix(999,2,1)
          cat(l," has NA! \n")
          break;
        }
        
        inv.lambda <- old.tau/old.sig.eps
        # A <- solve(diag.I + inv.lambda*est.poly.K)%*%(inv.lambda*est.poly.K)
        d_inv_tmp <- solve(dI + inv.lambda*d_poly_K)

        d_A <- d_inv_tmp %*% (inv.lambda*d_poly_K + dX %*% solve(t(dX) %*% d_inv_tmp %*% dX) %*% t(dX) %*% d_inv_tmp)
        #est.h <- matrix(A%*%scale.y,ncol=1)
        
        d_beta <- d_inv_XVX %*% t(dX) %*% d_inv_V %*% dy
        d_est_h <- d_A %*% dy - dX %*% d_beta
        
        
        #new.psi[2,] <- as.numeric(t(scale.y-est.h)%*%(scale.y-est.h)/(n-tr(A)))
        new.psi[2,] <- as.matrix(t(dy - d_A %*% dy) %*% (dy - d_A %*% dy)/(nrow(d_A)-sum(diag(d_A))))
        
        #loglik[iter] <- -0.5*(log(det(V))+t(scale.y)%*%P%*%scale.y)
        if(det(V) < 0){
          
          loglik[iter] <- as.matrix(-0.5*(log(det(inv_XVX))+t(dy) %*% d_P %*% dy))
       
        }else if(det(inv_XVX)<0){ 
       
          loglik[iter] <- as.matrix(-0.5*(log(det(V))+t(dy) %*% d_P %*% dy))
          
        }else if(det(V) <0 & det(inv_XVX)<0){
          
          loglik[iter] <- as.matrix(-0.5*(t(dy) %*% d_P %*% dy))
        }else{
          
          loglik[iter] <- as.matrix(-0.5*(log(det(V))+log(det(inv_XVX)) + t(dy) %*% d_P %*% dy))
          
        }
        
        
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
          
           d_poly_K <- (old.rho^(-1))*d_D_sum

          #if(all(off.xi!=0)){
            
           # for(j in 1:Off.index){
           #   est.poly.K <- est.poly.K+2*(off.xi[j])*off.square.z[[j]]
           # }
          #}
          
          save.loglik[l] <- max(loglik[which(!is.na(loglik))])
          save.psi[l,] <-  matrix(c(old.tau,opt.rho,opt.sig.eps),1,3)
          
          d_V <- opt.sig.eps*dI + opt.tau*d_poly_K
          d_inv_V <- solve(d_V)
          
          #est.h <-  opt.tau*est.poly.K%*%inv.V%*%(scale.y) 
          d_est_h <-  opt.tau*d_poly_K %*% d_inv_V %*% (dy-dX %*% d_beta) 
          

	  est.h <- as.matrix(d_est_h) 
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
 
if(comm.rank()==0){ 
  est.summary$opt.tau <- opt.tau
  est.summary$opt.rho <- opt.rho
  est.summary$opt.sig.eps <- opt.sig.eps

save.image(file="~/src/my_project/parallel_project/LSKM_MPI_Result.RData")

end_time <- proc.time()
ptm <- end_time - start_time

print(ptm)
}

finalize()
