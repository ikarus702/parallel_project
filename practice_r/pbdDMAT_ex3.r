library(pbdDMAT, quiet=TRUE)
init.grid()

nrows <- 35
ncols <- 1000e2

mn <- 10
sdd <- 100

bldim <- c(8,8)

if(comm.rank()==0){

x <- matrix(rnorm(n=nrows*ncols,mean=mn,sd=sdd),nrow=nrows,ncol=ncols)
b <- matrix(rnorm(n=ncols*ncols, mean=mn, sd=sdd),nrow=ncols,ncol=ncols)
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

cross_dxdb <- crossprod(t(dx),db) 

pbd_dxdb_prod <- as.matrix(cross_dxdb, proc.dest=0)

proc.time()-ptm

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
 r_xb_prod <- crossprod(t(x),b)
 
 print(proc.time()-ptm1)
 print(all.equal(pbd_dxdb_prod,r_xb_prod))
}

finalize()
