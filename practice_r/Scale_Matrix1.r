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
y.init.sd <- as.numeric(as.matrix(dy.init.sd))

comm.print(submatrix(dy.init.sd),all.rank=FALSE)
comm.print(dy.init.mean,all.rank=FALSE)

dX <- scale(dX,center=TRUE, scale=dX.init.sd);
dy <- (dy.init - mean(dy.init))/y.init.sd;

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

# Eigendecomposition

eigen.time <- comm.timer({
d_Eigenobject <- eigen(dK,only.values=FALSE,symmetric=TRUE);
})

if(comm.rank()==0){
	cat("Eigendecomposition time is \n")
}
comm.print(eigen.time, all.rank=FALSE)
## lambda search

lambda.time <- comm.timer({

n <- nrow(dy);
tol <- 10^(-3)*n;

# upper_bound_lambda


	U <- n;

	while(sum(d_Eigenobject$values / (d_Eigenobject$values + U)) < 1){
		U <- U-1
	}
	comm.print(U,all.rank=FALSE)
	 
# lower_bound_lambda

	q <- which.min(abs(d_Eigenobject$values - (max(d_Eigenobject$values)/1000)));

	L = .Machine$double.eps;

	while(sum(d_Eigenobject$values / (d_Eigenobject$values + L)) > q){
		L <- L + 0.5
	}
	comm.print(L,all.rank=FALSE)

# create new search values #

X1 <- L + 0.381966*(U - L);
X2 <- U - 0.381966*(U - L);

	d_diag_eigen <- as.ddmatrix(x=diag(d_Eigenobject$values),bldim=c(blm,blm));
	d_I <- as.ddmatrix(x=diag(1,n),bldim=c(blm,blm));
	
	
	d_inv_mat <- solve(d_diag_eigen + X1*d_I);
	d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));
	
	d_diag_inv_G <- as.ddmatrix(x=diag(diag(d_inv_G)),bldim=c(blm,blm));

	d_coeffs <- crossprod(d_inv_G,dy);
	LOOE1 <- crossprod(solve(d_diag_inv_G,d_coeffs));
	
	comm.print(submatrix(dy)[1:5],all.rank=FALSE)
	comm.print(diag(submatrix(d_diag_inv_G))[1:5],all.rank=FALSE)
	comm.print(submatrix(d_coeffs)[1:5], all.rank=FALSE)

	comm.print(submatrix(LOOE1),all.rank=FALSE)

	d_inv_mat <- solve(d_diag_eigen + X2*d_I);
	d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));
	
	d_diag_inv_G <- as.ddmatrix(x=diag(diag(d_inv_G)),bldim=c(blm,blm));

	d_coeffs <- crossprod(d_inv_G,dy);
	LOOE2 <- crossprod(solve(d_diag_inv_G,d_coeffs));
	

	comm.print(submatrix(d_coeffs)[1:5], all.rank=FALSE)
	comm.print(submatrix(LOOE2),all.rank=FALSE)

	d_diff_LE <- (LOOE1-LOOE2);
	diff_LE <- as.matrix(d_diff_LE);

	while(abs(diff_LE) > tol){

	if(diff_LE < 0){		

		U <- X2
		X2 <- X1
		X1 <- L + 0.381966*(U-L)
		LOOE2 <- LOOE1	
		d_inv_mat <- solve(d_diag_eigen + X1*d_I);
		d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));
		d_diag_inv_G  <- diag(d_inv_G);


		d_coeffs <- crossprod(d_inv_G,dy)
		LOOE1 <- crossprod(d_coeffs/d_diag_inv_G)
	}else{
		L <- X1
		X1 <- X2
		X2 <- U - (0.381966)*(U-L)
		LOOE1 <- LOOE2
		d_inv_mat <- solve(d_diag_eigen + X2*d_I);
		d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));
		d_diag_inv_G  <- diag(d_inv_G);


		d_coeffs <- crossprod(d_inv_G,dy)
		LOOE2 <- crossprod(d_coeffs/d_diag_inv_G)

	}
	
	d_diff_LE <- (LOOE1-LOOE2)
	diff_LE <- as.matrix(d_diff_LE)
 }	
	comm.print(X1,all.rank=FALSE)
	comm.print(X2,all.rank=FALSE)
	lambda <- ifelse(diff_LE < 0, X1 ,X2);
	
	LOOE <- crossprod(d_coeffs/d_diag_inv_G)
})

comm.print(lambda)








if(comm.rank()==0){
 



 #print(all.equal(scale.X,pbd_X))
 save.image("~/src/my_project/parallel_project/Scale_Matrix.RData") 
}
finalize()
