library(pbdBASE, lib.loc="~/R/lib")
library(pbdDMAT, quiet=TRUE,lib.loc="~/R/lib")
library(openblasctl)
library(pryr)
library(pbdIO, quiet=TRUE)
library(dmx)

openblas_set_num_threads(1)
args = commandArgs(trailingOnly = TRUE)
init.grid()


n_num <- as.numeric(args[1])
p_num <- as.numeric(args[2])
blm <- as.numeric(args[3])

sigma <- p_num

load(paste0("~/src/RData/NGK_PK_GP_Data_n",n_num,"_p",p_num,".RData"))

Data.num <- 1
rm(eps.list)
X <- as.matrix(data.list[[Data.num]][,-1])
y <- as.matrix(data.list[[Data.num]][,1])
n <- nrow(X)
d <- ncol(X)

scale.X <- scale(X)

r_I <- diag(1,n)
rm(data.list)

	comm.print("rank0 input rank") 
	comm.print(mem_used())
	comm.print("after input")


comm.print(mem_used())


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
dX.init.sd <- apply(dX.init,2,sd)
dy.init <- dy;
dy.init.mean <- mean(dy.init);
dy.init.sd <- sd(dy.init)
y.init.sd <- as.numeric(as.matrix(dy.init.sd))
X.init.sd <- as.numeric(as.matrix(dX.init.sd))

dX <- scale(dX.init)
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
rm(dX)
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


return.time1 <- comm.timer({

r_K <- as.matrix(dK)

kron_row <- comm.chunk(n, type="equal",form="vector")
j <- comm.rank()

vec_K <- matrix(as.vector(r_K[,kron_row]),ncol=1)

comm.print("vec_K dimension")
comm.print(dim(vec_K))



n_kron_row <- length(kron_row)

r_onevec1 <- matrix(rep(1,n_kron_row),ncol=1)
#r_coeffs <- as.matrix(d_coeffs,proc.dest=0);

r_coeffs <- as.matrix(d_coeffs)
r_kron_alpha <- kronecker(r_onevec1,r_coeffs)
rm(d_coeffs)
#if(comm.rank()==0){
#print("mem of r_kron_alpha is")
#print(object.size(r_kron_alpha))
#}
#barrier()


r_onevec2 <- matrix(rep(1,n),ncol=1)
tmp_I <- diag(1,n_kron_row)
kron_I <- kronecker(tmp_I,r_onevec2)

#kron_I <-  matrix(0,nrow=n_kron_row*n, ncol=n_kron_row)
#comm.print(kron_col[1], all.rank=TRUE)
#index <- 1:(n_kron_row*n) + c((kron_row[1]-1)*n)

#comm.print(index[1], all.rank=TRUE)
#kron_I[index,] <- tmp_kron_I
#rm(tmp_kron_I)


#d_mat <- ddmatrix(tmp.vec, ncol=1, bldim=c(blm,blm))
#d_kron_I[[i]] <- tmp.vec


comm.print("mem of tmp.mat is")
comm.print(object.size(kron_I))
#d_kron_I <- as.list(1:kron_col)

#comm.print(tmp.mat[c(index[1]+0:4),c(1:5)],all.rank=TRUE)
#comm.print(dim(tmp.mat),all.rank=TRUE)



#d_kron_I <- new("ddmatrix",Data=tmp.mat,dim=c(dim(tmp.mat)[1],allreduce(dim(tmp.mat)[2])), bldim=c(dim(tmp.mat)[1], comm.chunk(d, type="equal",form="bldim")), ldim=dim(tmp.mat), ICTXT=1)
#comm.print((d_kron_I),all.rank`=TRUE)

#comm.print("mem of d_mat is")
#comm.print(object.size(d_kron_I))


})
comm.print(mem_used())

if(comm.rank()==0){

#       	rm(r_kron_alpha)
#	rm(r_kron_I)
	print(mem_used())
	print("return time is \n")
}
comm.print(return.time1, all.rank=FALSE)



deriv.time <- comm.timer({

row.num <- comm.chunk(n, type="equal", form="vector")
comm.print(row.num, all.rank=TRUE)
comm.print(mem_used())
if(length(row.num)==0){

	r_derivmat <- NULL

}else{

rows <- cbind(rep(row.num,each=n),1:n);
r_dist <- scale.X[rows[,1],]-scale.X[rows[,2],];
comm.print(mem_used())
rm(rows); rm(scale.X)
comm.print(object.size(r_dist))

distk <- sweep(r_dist,1,vec_K,FUN="*");
rm(r_dist); rm(vec_K)
comm.print(mem_used())
comm.print(object.size(distk))


deriv_col <-sweep(distk,1,r_kron_alpha,FUN="*");
rm(distk);rm(r_kron_alpha)
comm.print(mem_used())
comm.print(object.size(deriv_col))

r_derivmat <- -(2/sigma)*(crossprod(kron_I, deriv_col))
rm(kron_I); rm(deriv_col)
comm.print(mem_used())
comm.print(object.size(r_derivmat))
}
barrier()

derivmat <- gather(r_derivmat)
rm(r_derivmat)
comm.print(length(derivmat))
comm.print(dim(derivmat[[1]]))

if(comm.rank()==0){

derivmat <- do.call(rbind, derivmat)

print(dim(derivmat))

avgderivmat <- apply(derivmat,2,mean)

}


})



if(comm.rank()==0){
	print(mem_used())
	cat("derivative time is \n")
}

comm.print(deriv.time, all.rank=FALSE)


## scale.back ##

scale.back.time <- comm.timer({

	d_fit_y <- d_fit_y %*% dy.init.sd + dy.init.mean
	d_vcov_c <- (y.init.sd^2)*d_covmatc
	d_vcov_fitted <- (y.init.sd^2)*d_vcovmat_yhat
	LOOE <- LOOE * y.init.sd
	if(comm.rank()==0){
	derivmat <- y.init.sd*derivmat/X.init.sd
	avgderivmat <- y.init.sd*avgderivmat/X.init.sd
	}
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

return.time2 <- comm.timer({

#r_K <- as.matrix(dK, proc.dest=0);
fitted <- as.matrix(d_fit_y, proc.dest=0);
X <- as.matrix(dX.init, proc.dest=0);
y <- as.matrix(dy.init, proc.dest=0);
R2 <- as.matrix(d_R2, proc.dest=0);
vcov.c <- as.matrix(d_vcov_c, proc.dest=0);
vcov.fitted <- as.matrix(d_vcov_fitted, proc.dest=0);
LOOE <- as.matrix(LOOE, proc.dest=0);
#derivmat <- as.matrix(d_derivmat, proc.dest=0)
#avgderiv <- as.matrix(d_avgderivmat, proc.dest=0)
})

 if(comm.rank()==0){
	
	print(mem_used())
	cat("return time is \n")
}
comm.print(return.time2, all.rank=FALSE)  




if(comm.rank()==0){
 
out <- list(K=r_K, 
	    coeffs=r_coeffs,
	    LOOE=LOOE,
	    fitted = fitted,
	    X=X,
	    y=y,
	    sigma=sigma,
	    lambda=lambda,
	    R2=R2,
	    vcov.c=vcov.c,
	    vcov.fitted=vcov.fitted,
	    derivatives=derivmat,
	    avgderivatives=avgderivmat	
	)

total.time <- sum(input.time[3],scale.time[3],kernel.time[3],eigen.time[3],lambda.time[3],
		fit.time[3],vcov.time[3],return.time1[3],deriv.time[3],
			R2.time[3],scale.back.time[3],return.time2[3])
print(total.time)


 #print(all.equal(scale.X,pbd_X))
 save(out, file=paste0("~/src/my_project/parallel_project/pbd_KRLS_p",p_num,".RData"))
}
finalize()
