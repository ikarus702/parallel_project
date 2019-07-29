library(pbdDMAT, quiet=TRUE)
library(matrixcalc)
library(openblasctl)

openblas_set_num_threads(1)

args = commandArgs(trailingOnly = TRUE)

Nnode = as.numeric(args[1])

init.grid()



bldim <- c(4,4)

#i <- 0
#for(ncols in ncol.cand){

#i <- i+1

if(comm.rank()==0){
X <- matrix(runif(100,2,4),nrow=10,ncol=10)
K <- matrix(rep(1,100),10,10)

}else{
X<-NULL
K <- NULL
}

if(comm.rank()==0){

cat("processing time for single node is\n")
single.time <- system.time({ 
rows <- cbind(rep(1:nrow(X),each=nrow(X)), 1:nrow(X));
distances <- X[rows[,1],] - X[rows[,2],]; 
for(k in 1:10){
	distk <- matrix(distances[,k],10,10,byrow=TRUE)
	L <- distk*K
}
})

#single.time <- proc.time()-ptm1
print(single.time)



} 
barrier()



ptm <- comm.timer({

dX <-as.ddmatrix(x=X,bldim=bldim);
dK <- as.ddmatrix(x=K, bldim=bldim);

comm.print(dK[1,1],all.rank=FALSE)


})

comm.print(ptm)

#cat(paste0("processing time for pbdDMAT",comm.rank()," is \n"))
#if(i==1){
#	mpi.time <- reduce(proc.time()-ptm)/Nnode
#}else{

#        mpi.time <- rbind(mpi.time, reduce(proc.time()-ptm)/Nnode)
#}



if(comm.rank()==0){
 
# print(pbd_K[1:5,1:5])
# print(r_K[1:5,1:5])
 
#print(all.equal(pbd_K,r_K))


 # save.image("~/src/my_project/parallel_project/Comp_Kernel_Prod.RData")
}
finalize()
