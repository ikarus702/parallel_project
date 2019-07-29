library(pbdDMAT, quiet=TRUE)
init.grid()

nrows <- 50e2
ncols <- 50e2

mn <- 10
sdd <- 100

bldim <- c(32,32)

if(comm.rank()==0){

x <- matrix(rnorm(n=nrows*ncols,mean=mn,sd=sdd),nrow=nrows,ncol=ncols)
b <- matrix(rnorm(n=ncols*nrows,mean=mn, sd=sdd),nrow=ncols,ncol=nrows)
}else{
x <-NULL
b <- NULL

}

y <- ddmatrix(rnorm(n=nrows*ncols, mean=mn, sd=sdd),nrow=nrows, ncol=ncols, bldim=bldim)
g <- ddmatrix(rnorm(n=ncols*nrows, mean=mn, sd=sdd), nrows=ncols, ncol=nrows, bldim=bldim)

ptm <- proc.time()


#dx <-as.ddmatrix(x=x,bldim=bldim)
#db <- as.ddmatrix(x=b, bldim=bldim)

#cat("matrix for ",comm.rank(),"\n")
#print(dx)


#print(dim(dx))

cat("processing time for pbdDMAT is \n")
dy_inv <- solve(t(y)%*%g)
solny <- solve(dy_inv,g)

pbd_dy_inv <- as.matrix(dy_inv, proc.dest=0)
pbd_solns <- as.matrix(solny, proc.dest=0)

proc.time()-ptm

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
 r_x_inv <- solve(t(x)%*%x)
 r_solns <- solve(r_x_inv,b)
 
 print(proc.time()-ptm1)
 print(all.equal(pbd_dy_inv,r_x_inv))
 print(all.equal(pbd_solns, r_solns))
}

finalize()
