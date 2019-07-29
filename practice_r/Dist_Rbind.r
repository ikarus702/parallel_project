library(pbdDMAT, quiet=TRUE)
library(matrixcalc)
library(openblasctl)

openblas_set_num_threads(1)
args = commandArgs(trailingOnly = TRUE)
init.grid()


n_num <- as.numeric(args[1])
p_num <- as.numeric(args[2])
blm <- as.numeric(args[3])

sigma <- p_num

if(comm.rank()==0){
load(paste0("~/src/RData/Unif_Semi_NGK_PK_Data_n",n_num,"_p",p_num,".RData"))

Data.num <- 1
Data <- data.list[[Data.num]]

y <- as.matrix(Data[,1])
X <- as.matrix(Data[,-c(1,2)])
n <- nrow(X)
d <- ncol(X)


}else{
y <-NULL
X <- NULL

}
if(comm.rank()==0){
r_I <- diag(1,500)
r_onevec <- matrix(rep(1,500),ncol=1)

r_kron_I <- kronecker(r_I,r_onevec)

 X.init <- X
      X.init.sd <- apply(X.init,2,sd)
      y.init <- y
      y.init.sd <- apply(y.init,2,sd)
      y.init.mean <- mean(y.init)
      X <- scale(X,center=TRUE,scale=X.init.sd)  

print(y.init.sd)
print(y.init.mean)  
y <- scale(y,center=y.init.mean,scale=y.init.sd)


K <- exp(-1*as.matrix(dist(X)^2)/sigma)

derivmat <- matrix(NA,nrow(X),ncol(X))
varavgderivmat <- matrix(NA,1,ncol(X))


 rows <- cbind(rep(1:nrow(X), each = nrow(X)), 1:nrow(X))
distances <- X[rows[,1],] - X[ rows[,2],] # d by n*n matrix of pairwise distances  
print(dim(distances))


 for(k in 1:d){       
       if(d==1){
             distk <-  matrix(distances,n,n,byrow=TRUE)
         } else {
             distk <-  matrix(distances[,k],n,n,byrow=TRUE) 
        }
L <- distk*K

# derivmat[,k] <- (-2/sigma)*L%*%out$coeff
         # variance for average derivative
# varavgderivmat[1,k] <- (1/n^2)*sum((-2/sigma)^2 * crossprod(L,vcovmatc%*%L))     
}
print(distances[1:5,1:5])
print(L[1:5,1:5])

}else{
r_kron_I <- NULL
r_kron_alpha <- NULL
}

barrier()




bldim <- c(blm,blm)

input.time <- comm.timer({

dX <-as.ddmatrix(x=X,bldim=bldim);
dy <- as.ddmatrix(x=y, bldim=bldim);

})
if(comm.rank()==0){
	cat("Input time is \n")
}

comm.print(input.time, all.rank=FALSE)

# scale vars

scale.time <- comm.timer({

dX.init <- dX;
dX.init.sd <- sd(dX.init);
dy.init <- dy;
dy.init.mean <- mean(dy.init);
dy.init.sd <- sd(dy.init)
dX <- scale(dX,center=TRUE, scale=dX.init.sd);
dy <- solve(dy.init.sd,dy.init - mean(dy.init));

})
if(comm.rank()==0){
	cat("scale time is \n")
}
comm.print(scale.time, all.rank=FALSE)

# kernel_matrix

kernel.time <- comm.timer({

dX2 <- dX^2;
d_F_vec <- apply(dX2,1,sum);
d_D1 <- sweep(-2*crossprod(t(dX)),1,d_F_vec,FUN="+");
d_D2 <- sweep(d_D1,2,d_F_vec,FUN="+");
dK <- exp(-d_D2/sigma);
#pbd_K <- as.matrix(dK, proc.dest=0);
#diag(pbd_K) <- 1;
})

if(comm.rank()==0){
	cat("Kernel time is \n")
}
comm.print(kernel.time, all.rank=FALSE)

	
### Rbind ###
#n <- nrow(dX)
#d <- ncol(dX)

#for(i in 1:n){
#	d_tmp_row1 <- dX[i,]
	

#	if(i==1){
		
#		d_Dist <- sweep(dX,MARGIN=2,STATS=t(d_tmp_row1))
		#comm.print(submatrix(d_Dist)[1:5,1:5],all.rank=FALSE)
#	}else{

#		d_Dist <- rbind(d_Dist, sweep(dX,2,t(d_tmp_row1)))
#	}


#}

#comm.print(dim(d_Dist),all.rank=FALSE)

# for(k in 1:d){       
#       if(d==1){
#             distk <-  matrix(d,n,n,byrow=TRUE)
#         } else {
#             distk <-  matrix(distances[,k],n,n,byrow=TRUE) 
#        }
#L <- distk*K


A <- as.ddmatrix(x=matrix(1:500^2,5*100,5*100),bldim=c(blm,blm))

dist.A <- as.ddmatrix(x=matrix(NA,nrow(A)^2,1),bldim=c(blm,blm))

#ptm <- comm.timer({

#for(j in 1:ncol(A)){
#tmp_A <- as.matrix(A[,j])

#d_A <- as.ddmatrix(kronecker(tmp_A,matrix(rep(1,nrow(A)),ncol=1)),bldim=c(blm,blm))


#for(i in 1:nrow(A)){
#	index <- c(25*(i-1)+1):c(25*i)
#
#	if(i==1){

#	dist.A <- t(A[,j] - d_A[1:25,])
	
#	}else{

#	dist.A <- rbind(dist.A, t(A[,j]-d_A[index,]))
#	}

#}


K <- as.ddmatrix(x=matrix(2,nrow(A),nrow(A)),bldim=c(blm,blm))
alpha <- as.ddmatrix(x=matrix(3,nrow(A),ncol=1),bldim=c(blm,blm))
#L <- K*dist.A


#}

#})

#comm.print(ptm, all.rank=FALSE)
#if(comm.rank()==0){
#print(dim(pbd_Dist))
#print(pbd_Dist[1:5,1:5])
#print(pbd_Dist[1001:1005,1:5])
#}

ptm1 <- comm.timer({

rows <- cbind(rep(1:nrow(A),each=nrow(A)),1:nrow(A));
d_dist <- A[rows[,1],]-A[rows[,2],];

})
comm.print(ptm1, all.rank=FALSE)


return.time <- comm.timer({


r_K <- as.matrix(K,proc.dest=0);

r_alpha <- as.matrix(alpha,proc.dest=0)

if(comm.rank()==0){
r_kron_alpha <- kronecker(r_onevec,r_alpha)
}
barrier()

vec_K <- as.ddmatrix(x=do.call(rbind,as.list(t(r_K))),bldim=c(blm,blm))
d_kron_I <- as.ddmatrix(x=r_kron_I,bldim=c(blm,blm))
d_kron_alpha <- as.ddmatrix(x=r_kron_alpha, bldim=c(blm,blm))

})

if(comm.rank()==0){
	print("return time is \n")
}
comm.print(return.time, all.rank=FALSE)

ptm2 <- comm.timer({


list_dist <- as.list(seq_len(ncol(A)))

#for(k in 1:ncol(A)){

distk <- sweep(d_dist,1,vec_K,FUN="*");
d_deriv_col <-sweep(distk,1,d_kron_alpha,FUN="*");
#deriv_col <- distk*d_kron_alpha
d_derivmat <- crossprod(d_kron_I,d_deriv_col)

#list_dist[[k]] <- deriv_vec

#if(k==1){
#	d_derivmat <- deriv_vec

#}else{
#	d_derivmat <- cbind(d_derivmat, deriv_vec)
#}


#}

comm.print(dim(distk),all.rank=FALSE)
comm.print(dim(d_deriv_col),all.rank=FALSE)
#comm.print(submatrix(list_dist[[2]]),all.rank=FALSE)

#derivmat <- do.call(cbind,list_dist)

comm.print(dim(d_derivmat),all.rank=FALSE)

})


#for(k in 1:2){

comm.print("ptm2",all.rank=FALSE)
# ptm2 <- comm.timer({
 
#	for(i in 1:nrow(A)){
#	index <- c(nrow(A)*(i-1)+1):c(nrow(A)*i)
	
#		list_dist[[i]] <- (d_dist[index,k])
#	if(i==1){
#			deriv.col <- d_K_prime %*% alpha
#	}else{
#		deriv.col <- deriv.col
#		comm.print(list_dist[[i]],all.rank=FALSE)
#	}			
#		}
#	distk <- do.call(rbind,list_dist)


#	comm.print(dim(distk),all.rank=FALSE)
#comm.print(dim(deriv.col),all.rank=FALSE)	
#})
comm.print(ptm2,all.rank=FALSE)


#ptm4 <- comm.timer({
#	if(k==1){
#		derivmat <- deriv.col
#	}else{
#	derivmat <-derivmat
#	}

#})
#comm.print("ptm4",all.rank=FALSE)
#comm.print(ptm4, all.rank=FALSE)

#}



comm.print(ptm2, all.rank=FALSE)

#comm.print(submatrix(derivmat), all.rank=FALSE)
#comm.print(dim(derivmat),all.rank=FALSE)

if(comm.rank()==0){

 #print(all.equal(scale.X,pbd_X))
 #save.image("~/src/my_project/parallel_project/Dist_Rbind.RData") 
}

finalize()
