library(pbdDMAT, quiet=TRUE)

args = commandArgs(trailingOnly=TRUE)
Nnode = as.numeric(args[1])

init.grid()

nrows <- 5e2
ncols <- 5e2

mn <- 10
sdd <- 100


if(comm.rank()==0){

x <- matrix(rnorm(n=nrows*ncols,mean=mn,sd=sdd),nrow=nrows,ncol=ncols)
b <- matrix(rnorm(n=ncols*2, mean=mn, sd=sdd),nrow=ncols,ncol=2)
}else{
x <-NULL
b <- NULL

}

Block.size <- c(2,4,8,16,32,64)
i <- 0

for(blm in Block.size){

i <- i+1

bldim <- c(blm, blm)
ptm <- proc.time()

dx <-as.ddmatrix(x=x,bldim=bldim)
db <- as.ddmatrix(x=b, bldim=bldim)

dx_inv <- solve(t(dx)%*%dx)
solnx <- solve(dx_inv,db)

pbd_dx_inv <- as.matrix(dx_inv, proc.dest=0)
pbd_solns <- as.matrix(solnx, proc.dest=0)

if(i==1){
 mpi.time <- reduce(proc.time()-ptm)/Nnode
}else{
  mpi.time <- rbind(mpi.time, reduce(proc.time()-ptm)/Nnode)
} 

}

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
 r_x_inv <- solve(t(x)%*%x)
 r_solns <- solve(r_x_inv,b)
 
 single.time <- proc.time()-ptm1
 print(all.equal(pbd_dx_inv,r_x_inv))
 print(all.equal(pbd_solns, r_solns))

 save.image("~/src/my_project/parallel_project/Comp_Matmult_Block.RData")
}


finalize()
