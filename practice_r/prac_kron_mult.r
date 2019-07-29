library(pbdBASE, lib.loc="~/R/lib")
library(pbdDMAT, quiet=TRUE,lib.loc="~/R/lib")
library(openblasctl)
library(pryr)
library(pbdIO, quiet=TRUE)

openblas_set_num_threads(1)
args = commandArgs(trailingOnly = TRUE)
init.grid()


n_num <- as.numeric(args[1])
p_num <- as.numeric(args[2])
blm <- as.numeric(args[3])

sigma <- p_num

if(comm.rank()==0){
load(paste0("~/src/RData/NGK_PK_GP_Data_n",n_num,"_p",p_num,".RData"))

Data.num <- 1
rm(eps.list)
X <- as.matrix(data.list[[1]][,-1])
y <- as.matrix(data.list[[1]][,1])
n <- nrow(X)
d <- ncol(X)

r_I <- diag(1,n)
r_onevec <- matrix(rep(1,n),ncol=1)
#r_kron_I <- kronecker(r_I,r_onevec)
rm(data.list)
print(object.size(r_onevec))
	print("rank0 input rank") 
	print(mem_used())
	print("after input")
}else{
data.list <- NULL
Data <- NULL
y <-NULL
X <- NULL
r_I <- NULL
r_onevec <- NULL
r_kron_I <- NULL
}


barrier()


comm.print(mem_used(),all.rank=TRUE)





bldim <- c(blm,blm)

input.time <- comm.timer({

dX <-as.ddmatrix(x=X,bldim=bldim);
dy <- as.ddmatrix(x=y, bldim=bldim);

comm.print("mem of dX is ")
comm.print(object.size(dX))
comm.print("mem of dy is ")
comm.print(object.size(dy))


})
if(comm.rank()==0){
	print(mem_used())
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
X.init.sd <- as.numeric(as.matrix(dX.init.sd))

dX <- scale(dX,center=TRUE, scale=dX.init.sd);
dy <- (dy.init - mean(dy.init))/y.init.sd;

})
if(comm.rank()==0){
	
	print(mem_used())
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

comm.print("mem of dK is")
comm.print(object.size(dK))
#pbd_K <- as.matrix(dK, proc.dest=0);
#diag(pbd_K) <- 1;
})
rm(dX2)
rm(d_D1)
rm(d_D2)
rm(d_F_vec)
if(comm.rank()==0){
	
	print(mem_used())
	cat("Kernel time is \n")

}
comm.print(kernel.time, all.rank=FALSE)

# Eigendecomposition

eigen.time <- comm.timer({
d_Eigenobject <- eigen(dK,only.values=FALSE,symmetric=TRUE);
})

if(comm.rank()==0){
	print(mem_used())

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
                d_inv_mat <- solve(d_diag_eigen + lambda*d_I);
                d_inv_G <- tcrossprod(crossprod(t(d_Eigenobject$vectors),d_inv_mat),(d_Eigenobject$vectors));
                d_diag_inv_G <- as.ddmatrix(x=diag(diag(d_inv_G)),bldim=c(blm,blm));

        d_coeffs <- crossprod(d_inv_G,dy);

	LOOE <- crossprod(solve(d_diag_inv_G,d_coeffs))
})

if(comm.rank()==0){
	
	print(mem_used())
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
	comm.print(object.size(d_covmatc))

	d_vcovmat_yhat <- crossprod(dK,crossprod(d_covmatc,dK))
	comm.print(object.size(d_vcovmat_yhat))
	rm(d_inv_mat)
})

if(comm.rank()==0){
	
	print(mem_used())
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
	
	print(mem_used())
	cat("row time is \n")
}
rm(rows)
rm(dX)
comm.print(row.time, all.rank=FALSE)

return.time1 <- comm.timer({


r_K <- as.matrix(dK,proc.dest=0);

r_coeffs <- matrix(c(1:n),ncol=1)

if(comm.rank()==0){
r_kron_alpha <- kronecker(r_onevec,r_coeffs)

print("mem of r_kron_alpha is")
print(object.size(r_kron_alpha))
}else{
r_kron_alpha <- NULL
}
barrier()


vec_K <- as.ddmatrix(x=do.call(rbind,as.list(t(r_K))),bldim=c(blm,blm))
#d_kron_I <- as.ddmatrix(x=r_kron_I,bldim=c(blm,blm))

d_onevec <- as.ddmatrix(x=r_onevec, bldim=c(blm,blm));



#for(j in 1:c(d/100)){

j <- comm.rank()
tmp.vec <-  rep(0,n^2)

kron_col <- comm.chunk(d, type="equal")
comm.print(kron_col, all.rank=TRUE)
d_mat <- as.list(1:kron_col)



for(i in 1:kron_col){

index <- c((i-1)*n+1+kron_col*j*n):c(i*n+kron_col*j*n)
tmp.vec[index] <- 1
#d_mat <- ddmatrix(tmp.vec, ncol=1, bldim=c(blm,blm))
#d_kron_I[[i]] <- tmp.ve

d_mat[[i]] <- tmp.vec

tmp.vec <- rep(0,n^2)
}

tmp.mat <- do.call(cbind, d_mat)
rm(d_mat)
comm.print("mem of tmp.mat is")
comm.print(object.size(tmp.mat))
#d_kron_I <- as.list(1:kron_col)

#comm.print(tmp.mat[c(1:5)+kron_col*j*n,c(1:5)],all.rank=TRUE)
#comm.print(dim(tmp.mat),all.rank=TRUE)


d_kron_I <- new("ddmatrix",Data=tmp.mat,dim=c(dim(tmp.mat)[1],allreduce(dim(tmp.mat)[2])), bldim=c(dim(tmp.mat)[1], comm.chunk(d, type="equal",form="bldim")), ldim=dim(tmp.mat), ICTXT=1)
#comm.print((d_kron_I),all.rank=TRUE)
rm(tmp.mat)

comm.print("mem of d_mat is")
comm.print(object.size(d_kron_I))

d_kron_I <- as.blockcyclic((d_kron_I),bldim=c(blm,blm))

comm.print(d_kron_I)
comm.print("mem of blockcyclic is")
comm.print(object.size(d_kron_I))


#}

#d_kron_I <- do.call(cbind,d_tot_mat)


d_kron_alpha <- as.ddmatrix(x=r_kron_alpha, bldim=c(blm,blm))
comm.print("mem of kron_alpha is")
comm.print(object.size(d_kron_alpha))
#rm(d_tot_mat)

})
comm.print(mem_used(),all.rank=TRUE)

if(comm.rank()==0){
       	rm(r_kron_alpha)
#	rm(r_kron_I)
	print(mem_used())
	print("return time is \n")
}
comm.print(return.time1, all.rank=FALSE)


mult_result <- crossprod(d_kron_I,d_kron_alpha)
r_mult <- as.matrix(mult_result,proc.dest=0)



if(comm.rank()==0){
 print(r_mult)
 #print(all.equal(scale.X,pbd_X))
 save.image(file="~/src/my_project/parallel_project/Kron_Prac.RData")
}
finalize()
