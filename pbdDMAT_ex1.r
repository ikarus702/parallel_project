library(pbdDMAT, quiet=TRUE)
suppressMessages(library(parallel))

args = commandArgs(trailingOnly = TRUE)

cores = as.numeric(args[1])
host = system("hostname", intern = TRUE)

mc.function = function(x) {
    ## Put code for mclapply cores here
    Sys.getpid() # returns process id
}

mycores = mclapply(1:cores, mc.function, mc.cores = cores)
mc = do.call(paste, mycores) # combines results from mclapply

## same cores available for OpenBLAS (see openblasctl package)
##            or for other OpenMP enabled codes

## Now report what happened and where
comm.cat("Hello World from rank", comm.rank(), "on host", host, "with", cores, "cores allocated.\n",
         "               Processes:", mc, "\n", quiet = TRUE, all.rank = TRUE)


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
ptm <- proc.time()


dx <-as.ddmatrix(x=x,bldim=bldim)
db <- as.ddmatrix(x=b, bldim=bldim)

cat("matrix for ",comm.rank(),"\n")
print(dx)


print(dim(dx))

cat("processing time for pbdDMAT is \n")
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
 print(pbd_solns[1:5,1:5])
 print(r_solns[1:5,1:5])

 print(all.equal(pbd_dx_inv,r_x_inv))
 print(all.equal(pbd_solns, r_solns))
}

finalize()
