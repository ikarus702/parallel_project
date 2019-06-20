
library(mvtnorm)

Data.num <- 50 ## number of datasets
p.val <- 20
q.val <- 180
n.val <- 200
xi.val1 <- 1
xi.val2 <- 1
off.val1 <- 0; off.val2 <- 0


data.list <- list()
eps.list <- list()



gen.NGK <- function(Data.num,p=p.val,q=q.val,n=n.val,xi1=xi.val1,xi2=xi.val2,off1=off.val1,off2=off.val2){
  
  ############# Setting with Xi matrix ##############
  
  Tot.p <- p + q 
  
  true.xi <- xi1
  true.sig.eps <- 1
  true.tau <- 10
  
  mat.I <- diag(n)
  
  zero.vec1 <- matrix(0,p,1)
  zero.vec2 <- matrix(0,q,1)
  zero.vec3 <- matrix(0,n,1)
  
  set.seed(123)
  
  data.list <- list()
  eps.list <- list()
  
  ######### correlation ###############
  
  model1=function(dim,xi.val,off.diag)#the function to generate a tridiagonal strcture
  {
    prec=diag(dim)*xi.val
    
    for(i in 1:dim){
      for(j in 1:dim){
        if(abs(i-j)==1){
          prec[i,j] <- off.diag
        }
      }
    }
    return(prec)
  }
  xi.val1 <- xi1
  xi.sigma1 <- model1(p,xi.val=xi1,off.diag=off1)
  xi.sigma2 <- model1(q,xi.val=xi2,off.diag =off2)
  
  
  load("NGK_PK_GP_Data_p200.RData")

  for(Data.index in 1:Data.num){
    
    y <- matrix(data.list[[Data.index]][,1],nrow=n)
    z <- t(data.list[[Data.index]][,-1])
    eps <- eps.list[[Data.index]]
    X <- runif(n,-1,1)
    X <- matrix(sort(X),ncol=1)
    beta <- 1
    
    
    Scale.z <- matrix(0,Tot.p,n)
    Square.z <- list()
   
    ###### z1 : sig, z2 : not sig
  
    
    for(j in 1:Tot.p){
      Scale.z[j,] <- scale(z[j,])
    }
    
    poly.K <- matrix(NA,n,n) 
    
    for(l in 1:n){
      for(m in 1:n){
        
        tmp.vec.l <- matrix(Scale.z[1:p,l],nrow=p)
        tmp.vec.m <- matrix(Scale.z[1:p,m],nrow=p)
        
        poly.K[l,m] <- as.numeric(t(tmp.vec.l)%*%xi.sigma1%*%(tmp.vec.m))
        
      }
    }
    
    y <- matrix(X%*%beta+y,nrow=n)
    eps.list[[Data.index]] <<- eps
    data.list[[Data.index]] <<- data.frame(cbind(y,X,t(z)))
  }
  
  save.image("Unif_Semi_NGK_PK_Data_p200.RData")
}


gen.NGK(Data.num,p.val,q.val,n.val,xi.val1,xi.val2,off.val1,off.val2)
