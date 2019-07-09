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

	 
# lower_bound_lambda

	q <- which.min(abs(d_Eigenobject$values - (max(d_Eigenobject$values)/1000)));

	L = .Machine$double.eps;

	while(sum(d_Eigenobject$values / (d_Eigenobject$values + L)) > q){
		L <- L + 0.5
	}

# create new search values #

X1 <- L + 0.381966*(U - L);
X2 <- U - 0.381966*(U - L);

	d_diag_eigen <- as.ddmatrix(x=diag(d_Eigenobject$values),bldim=c(blm,blm));
	d_I <- as.ddmatrix(x=diag(1,n),bldim=c(blm,blm));
	
	
	d_inv_mat <- solve(d_diag_eigen + X1*d_I);
	d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));

	d_diag_inv_G <- as.ddmatrix(x=diag(diag(d_inv_G)),bldim=c(blm,blm));

	d_coeffs <- crossprod(d_inv_G,dy);
	LOOE1 <- crossprod(solve(d_diag_inv_G, d_coeffs));

	d_inv_mat <- solve(d_diag_eigen + X2*d_I);
	d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));
	d_diag_inv_G <- as.ddmatrix(x=diag(diag(d_inv_G)),bldim=c(blm,blm));

	d_coeffs <- crossprod(d_inv_G,dy);
	LOOE2 <- crossprod(solve(d_diag_inv_G, d_coeffs));

	
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
		d_diag_inv_G <- as.ddmatrix(x=diag(diag(d_inv_G)),bldim=c(blm,blm));

	d_coeffs <- crossprod(d_inv_G,dy);
	LOOE1 <- crossprod(solve(d_diag_inv_G, d_coeffs));



	}else{
		L <- X1
		X1 <- X2
		X2 <- U - (0.381966)*(U-L)
		LOOE1 <- LOOE2
		d_inv_mat <- solve(d_diag_eigen + X2*d_I);
		d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));
		d_diag_inv_G <- as.ddmatrix(x=diag(diag(d_inv_G)),bldim=c(blm,blm));

	d_coeffs <- crossprod(d_inv_G,dy);
	LOOE2 <- crossprod(solve(d_diag_inv_G, d_coeffs));


	}
	
	d_diff_LE <- (LOOE1-LOOE2)
	diff_LE <- as.matrix(d_diff_LE)
 }	
	lambda <- ifelse(diff_LE < 0, X1 ,X2);
	lambda <- as.numeric(lambda);
	LOOE <- crossprod(solve(d_diag_inv_G,d_coeffs))
})

if(comm.rank()==0){
	cat("search time for lambda is \n")
}

comm.print(lambda.time, all.rank=FALSE)





### fitted y ###

fit.time <- comm.timer({

	d_fit_y <- dK %*% d_coeffs

})

if(comm.rank()==0){
	cat("Fitting time is \n")
}

comm.print(fit.time, all.rank=FALSE)

### var-covariance matrix for coeffs ###

vcov.time <- comm.timer({

	d_sigmasq <- as.numeric(as.matrix(crossprod(dy-d_fit_y)/n));
	
	d_inv_mat <- solve(d_diag_eigen + lambda*d_I);
	d_covmatc <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_sigmasq*crossprod(d_inv_mat)),(d_Eigenobject$vectors));

	d_vcovmat_yhat <- crossprod(dK,crossprod(d_covmatc,dK))

})

if(comm.rank()==0){
	cat("Covariance time is \n")
}
comm.print(vcov.time,all.rank=FALSE)


### derivative #####

#deriv.time <- comm.timer({



n <- nrow(dX)
d <- ncol(dX)

row.time <- comm.timer({

rows <- cbind(rep(1:n,each=n),1:n);
d_dist <- dX[rows[,1],]-dX[rows[,2],];

})
if(comm.rank()==0){
	cat("row time is \n")
}
comm.print(row.time, all.rank=FALSE)
for(k in 1:10){

dd.time <- comm.timer({

        distk <- as.ddmatrix(x=matrix(as.matrix(d_dist[,k]),nrow=n,ncol=n,byrow=TRUE),bldim=c(blm,blm))


}) 


if(comm.rank()==0){
	cat("ddmatrix time is \n")
}

comm.print(dd.time, all.rank=FALSE)

deriv.time <- comm.timer({
d_K_prime <- distk * dK

        if(k==1){
                d_derivmat <- (-2/sigma)*d_K_prime %*% d_coeffs
        	d_varavgderivmat <- (1/n^2)*sum((-2/sigma)^2 * crossprod(d_K_prime,d_covmatc%*% d_K_prime)) 
	}else{
                d_derivmat <- cbind(d_derivmat,(-2/sigma)*d_K_prime %*% d_coeffs)
		d_varavgderivmat <- cbind(d_varavgderivmat, (1/n^2)*sum((-2/sigma)^2 * crossprod(d_K_prime,d_covmatc%*% d_K_prime))) 	

        }
})
if(comm.rank()==0){
	cat("derivative time is \n")
}

comm.print(deriv.time, all.rank=FALSE)



	if(comm.rank()==0){
		print(k)
	}

}


#})


#comm.print(deriv.time, all.rank=FALSE)
comm.print(dim(d_derivmat),all.rank=FALSE)
comm.print(class(d_derivmat),all.rank=FALSE)
comm.print(d_varavgderivmat,all.rank=FALSE)


## scale.back ##

scale.back.time <- comm.timer({

	d_fit_y <- d_fit_y %*% dy.init.sd + dy.init.mean
	d_vcov_c <- (y.init.sd^2)*d_covmatc
	d_vcov_fitted <- (y.init.sd^2)*d_vcovmat_yhat
	LOOE <- LOOE * y.init.sd
})

if(comm.rank()==0){
	cat("scaling-back time is \n")
}
comm.print(scale.back.time, all.rank=FALSE)

## R-squared ####

R2.time <- comm.timer({
	
	d_R2 <- 1-(crossprod(dy.init-d_fit_y)/(n*dy.init.sd^2))

})

if(comm.rank()==0){
	cat("R2 time is \n")
}
comm.print(R2.time,all.rank=FALSE)





## return ###

return.time <- comm.timer({

K <- as.matrix(dK, proc.dest=0);
coeffs <- as.matrix(d_coeffs, proc.dest=0);
fitted <- as.matrix(d_fit_y, proc.dest=0);
X <- as.matrix(dX.init, proc.dest=0);
y <- as.matrix(dy.init, proc.dest=0);
R2 <- as.matrix(d_R2, proc.dest=0);
vcov.c <- as.matrix(d_vcov_c, proc.dest=0);
vcov.fitted <- as.matrix(d_vcov_fitted, proc.dest=0);
LOOE <- as.matrix(LOOE, proc.dest=0)


})

 if(comm.rank()==0){
	cat("return time is \n")
}
comm.print(return.time, all.rank=FALSE)  




if(comm.rank()==0){
 
out <- list(K=K, 
	    coeffs=coeffs,
	    LOOE=LOOE,
	    fitted = fitted,
	    X=X,
	    y=y,
	    sigma=sigma,
	    lambda=lambda,
	    R2=R2,
	    vcov.c=vcov.c,
	    vcov.fitted=vcov.fitted
	)

#print(out)


 #print(all.equal(scale.X,pbd_X))
 save.image("~/src/my_project/parallel_project/pbd_KRLS.RData") 
}
finalize()
