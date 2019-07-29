library(pbdDMAT, quiet=TRUE)
library(matrixcalc)

args = commandArgs(trailingOnly = TRUE)

Nnode = as.numeric(args[1])

init.grid()
nrows <- 5000
ncols <- 5000

mn <- 10
sdd <- 100


bldim <- c(16,16)

#i <- 0
#for(ncols in ncol.cand){

#i <- i+1

if(comm.rank()==0){
##### AR(2): 2nd-order neighborhood matrix ##############
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
Z <- matrix(rnorm(n=nrows*ncols,mean=mn,sd=sdd),nrow=nrows,ncol=ncols)



#Gamma <- diag(round(runif(n=ncols, 0.1, 10),2))
#Gamma[upper.tri(Gamma)] <- round(runif(n=(ncols*(ncols-1)/2),1,10),2)

Omega <- model2(dim=ncols,1,0.5)
is.positive.semi.definite(Omega)


}else{
model2 <- NULL
model3 <- NULL
Z <-NULL
Omega <- NULL

}

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
single.time <- system.time({ r_K <- Z%*%Omega%*%t(Z)})

#single.time <- proc.time()-ptm1
print(single.time)



} 
barrier()



ptm <- comm.timer({

dZ <-as.ddmatrix(x=Z,bldim=bldim);
dOmega <- as.ddmatrix(x=Omega, bldim=bldim);


#dGamma2 <- crossprod(dGamma,dGamma)
#Gamma2 <- as.matrix(dGamma2,proc.dest=0)

#if(comm.rank()==0){
#print(all.equal(Omega,Gamma2))

#}


dK <- dZ %*% dOmega %*% t(dZ);


pbd_K <- as.matrix(dK, proc.dest=0);
})

comm.print(ptm)

#cat(paste0("processing time for pbdDMAT",comm.rank()," is \n"))
#if(i==1){
#	mpi.time <- reduce(proc.time()-ptm)/Nnode
#}else{

#        mpi.time <- rbind(mpi.time, reduce(proc.time()-ptm)/Nnode)
#}



if(comm.rank()==0){
 
 print(pbd_K[1:5,1:5])
 print(r_K[1:5,1:5])
 
print(all.equal(pbd_K,r_K))


  save.image("~/src/my_project/parallel_project/Comp_Kernel_Prod.RData")
}
finalize()
