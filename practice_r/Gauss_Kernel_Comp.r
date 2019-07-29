library(pbdDMAT, quiet=TRUE)
library(matrixcalc)

args = commandArgs(trailingOnly = TRUE)

Nnode = as.numeric(args[1])

init.grid()
nrows <- 100
ncols <- 100

mn <- 0
sdd <- 1


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

Omega <- diag(rep(c(1,1,1,1),each=c(ncols/4)),ncols)
is.positive.semi.definite(Omega)
print(diag(Omega)[1:5])

}else{
model2 <- NULL
model3 <- NULL
Z <-NULL
Omega <- NULL

}

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
single.time <- system.time({ 
new.Z <- Z %*% sqrt(Omega);
new.Dist.z <- as.matrix(dist(new.Z,diag=TRUE,upper=TRUE));
r_K <- matrix(exp(-new.Dist.z^2),ncol=ncols,byrow=TRUE);
})

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


d_new_Z <- dZ %*% sqrt(dOmega);
d_new_Z2 <- d_new_Z^2;
d_F_vec <- apply(d_new_Z2,1,sum);
d_D1 <- sweep(-2*crossprod(t(d_new_Z)),1,d_F_vec,FUN="+");
d_dist_Z2 <- sweep(d_D1,2,d_F_vec,FUN="+");
dK <- exp(-d_dist_Z2);
pbd_K <- as.matrix(dK, proc.dest=0);
diag(pbd_K) <- 1

})

comm.print(ptm)

#cat(paste0("processing time for pbdDMAT",comm.rank()," is \n"))
#if(i==1){
#	mpi.time <- reduce(proc.time()-ptm)/Nnode
#}else{

#        mpi.time <- rbind(mpi.time, reduce(proc.time()-ptm)/Nnode)
#}



if(comm.rank()==0){
 print(dim(pbd_K))
 print(pbd_K[1:5,1:5])
 print(dim(r_K))
 print(r_K[1:5,1:5])
 
print(all.equal(pbd_K,r_K))


  save.image("~/src/my_project/parallel_project/Gauss_Kernel.RData")
}
finalize()
