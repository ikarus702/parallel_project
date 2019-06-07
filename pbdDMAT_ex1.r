library(pbdDMAT, quiet=TRUE)
init.grid()

nrows <- 10e2
ncols <- 10e2

mn <- 10
sdd <- 100

bldim <- c(250,250)

if(comm.rank()==0){

x <- matrix(rnorm(n=nrows*ncols,mean=mn,sd=sdd),nrow=nrows,ncol=ncols)
b <- matrix(rnorm(n=ncols*2, mean=mn, sd=sdd),nrow=ncols,ncol=2)
}else{
x <-NULL
b <- NULL

}

dx <-as.ddmatrix(x=x,bldim=bldim)
db <- as.ddmatrix(x=b, bldim=bldim)

cat("matrix for ",comm.rank(),"\n")
print(dx)


print(dim(dx))

cat("processing time for pbdDMAT is \n")
ptm <- proc.time()

dx_inv <- solve(t(dx)%*%dx)
solnx <- solve(dx_inv,db)

pbd_dx_inv <- as.matrix(dx_inv, proc.dest=0)
pbd_solns <- as.matrix(solnx, proc.dest=0)

proc.time()-ptm

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
 r_x_inv <- solve(t(x)%*%x)
 r_solns <- solve(r_x_inv,b)
 
 print(proc.time()-ptm1)
 print(all.equal(pbd_dx_inv,r_x_inv))
 print(all.equal(pbd_solns, r_solns))
}

finalize()
