library(pbdDMAT, quiet=TRUE)
library(matrixcalc)

args = commandArgs(trailingOnly = TRUE)

Nnode = as.numeric(args[1])

init.grid()

nrows <- as.numeric(args[2])
ncols <- as.numeric(args[3])
Block.size <- c(16)

mn <- 10
sdd <- 100


###### AR(2): 2nd-order neighborhood matrix ##############

model2 <- function(dim, xi.val, off.diag){
	
	prec <- diag(1,dim)*xi.val

	for(i in 1:dim){
		for(j in 1:dim){
			if(abs(i-j)==1){
				prec[i,j] <- off.diag
			}else if(abs(i-j)==2){
				prec[i,j] <- off.diag^2
			}
		}
	}

return(prec)
}


################ full off-diagonal matrix ##########################

model3 <- function(dim,xi.val,off.min, off.max){

	prec <- diag(1, dim)*xi.val
	
	for(i in 1:dim){
		for(j in i:dim){
			if(i!=j){
				off.val <- runif(1,off.min,off.max)
				prec[i,j] <- off.val
				prec[j,i] <- off.val
			}
		}
	}
	return(prec)
}
if(comm.rank()==0){

Z <- matrix(rnorm(n=nrows*ncols,mean=mn,sd=sdd),nrow=nrows,ncol=ncols)



#Gamma <- diag(round(runif(n=ncols, 0.1, 10),2))
#Gamma[upper.tri(Gamma)] <- round(runif(n=(ncols*(ncols-1)/2),1,10),2)

Omega <- model2(dim=ncols,1,0.5)
is.positive.semi.definite(Omega)


}else{
Z <-NULL
Omega <- NULL

}



for(i in 1:100){



bldim <- c(Block.size,Block.size)

ptm <- proc.time()

dZ <-as.ddmatrix(x=Z,bldim=bldim)
dOmega <- as.ddmatrix(x=Omega, bldim=bldim)


#dGamma2 <- crossprod(dGamma,dGamma)
#Gamma2 <- as.matrix(dGamma2,proc.dest=0)

#if(comm.rank()==0){
#print(all.equal(Omega,Gamma2))

#}


cat("processing time for pbdDMAT is \n")
print(proc.time()-ptm)

cross_dZdOm <- crossprod(t(dZ),dOmega) 
dK <- crossprod(t(cross_dZdOm), t(dZ))


pbd_K <- as.matrix(dK, proc.dest=0)



if(i==1){
	mpi.time.mat <- reduce(proc.time()-ptm)/Nnode

}else{

        mpi.time.mat <- rbind(mpi.time.mat, reduce(proc.time()-ptm)/Nnode)
}

if(i==100){
 mpi.time <- apply(mpi.time.mat,2,mean)
}


}

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
 r_K <- Z%*%Omega%*%t(Z)
  
 print(pbd_K[1:5,1:5])
 print(r_K[1:5,1:5])
 single.time <- proc.time()-ptm1




 print(all.equal(pbd_K,r_K))

  save.image(paste0("~/src/my_project/parallel_project/Comp_Kernel_Prod_",ncols,"_Col.RData"))
}
finalize()
