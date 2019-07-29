library(pbdDMAT, quiet=TRUE)
library(matrixcalc)

args = commandArgs(trailingOnly = TRUE)

Nnode = as.numeric(args[1])

init.grid()



bldim <- c(4,4)

#i <- 0
#for(ncols in ncol.cand){

#i <- i+1

if(comm.rank()==0){
	
load("~/src/RData/Unif_Semi_NGK_PK_p70_Data.RData")

Data.num <- 1
 result.summary <- list()
  est.summary <- list()
  est.prec.summary <- list()

  ######## Import dataset from the list

  n <- nrow(data.list[[Data.num]])
  y <- matrix(data.list[[Data.num]][,1],nrow=n)
  X <- matrix(data.list[[Data.num]][,2],nrow=n)
  z <- t(data.list[[Data.num]][,-c(1,2)])
  eps <- eps.list[[Data.num]]
  true.f <- matrix((unlist(y - eps)),nrow=n)
  Tot.p <- nrow(z)


  scale.y <- y-mean(y)
  X <- X -mean(X)
  scale.f <- true.f - mean(true.f)

  diag.I <- diag(1,n)

}else{
z <- NULL
y <- NULL
X <- NULL
diag.I <- NULL

}

rho <- 1
tau <- 10
sigma.eps <- 1

if(comm.rank()==0){

cat("processing time for single node is\n")
 ptm1 <- proc.time()
single.time <- system.time({ 
r_D_sum <- t(z)%*%z;

print(dim(r_D_sum));
r_K <- rho*r_D_sum;
r_V <- sigma.eps*diag.I + tau*r_K;
r_inv_V <- solve(r_V);
r_inv_XVX <- solve(t(X) %*% r_inv_V %*% (X));
r_P <- r_inv_V - r_inv_V %*% X %*% r_inv_XVX %*% t(X) %*% r_inv_V
r_g <- (-0.5*sum(diag(r_V %*% r_P))+0.5*t(y) %*% r_P %*% r_V %*% r_P %*% y)

})


#single.time <- proc.time()-ptm1
print(single.time)

} 
barrier()



ptm <- comm.timer({

dz <-as.ddmatrix(x=z,bldim=bldim);
dy <- as.ddmatrix(x=y, bldim=bldim);
dX <- as.ddmatrix(x=X, bldim=bldim);
dI <- as.ddmatrix(x=diag.I, bldim=bldim);

d_D_sum <- crossprod(dz,dz);
d_K <- rho*d_D_sum;
d_V <- sigma.eps*dI + tau*d_K;
d_inv_V <- solve(d_V);
d_inv_XVX <- solve(t(dX) %*% d_inv_V %*% dX);
d_P <- d_inv_V - d_inv_V %*% dX %*% d_inv_XVX %*% t(dX) %*% d_inv_V;
d_g_psi <- (-0.5*sum(diag(d_V %*% d_P))+0.5*t(dy) %*% d_P %*% d_V %*% d_P %*% dy)

print(submatrix(d_g_psi))

pbd_D_sum <- as.matrix(d_D_sum, proc.dest=0);
pbd_K <- as.matrix(d_K, proc.dest=0);
pbd_P <- as.matrix(d_P, proc.dest=0);
pbd_g <- as.matrix(d_g_psi, proc.dest=0)

})

comm.print(ptm)




if(comm.rank()==0){
 
 print(pbd_K[1:5,1:5])
 print(r_K[1:5,1:5])
 
cat("Are the kernels same? \n")
print(all.equal(pbd_K,r_K))

cat("Are D_sums same? \n")
print(all.equal(pbd_D_sum, r_D_sum))

cat("Are P matrices same? \n")
print(all.equal(pbd_P, r_P))

cat("Are g same? \n")
print(all.equal(pbd_g, r_g))

  save.image("~/src/my_project/parallel_project/Prac_Kernel_Prod.RData")
}
finalize()
