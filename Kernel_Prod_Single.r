library(parallel)

.Internal(setMaxNumMathThreads(1)); .Internal(setNumMathThreads(1))

nrow.cand <- 1000*c(5)
ncols <- 5000

mn <- 10
sdd <- 100

###### AR(2): 2nd-order neighborhood matrix ##############

model2 <- function(dim, xi.val, off.diag){
	
	prec <- diag(1,dim)*xi.val

	for(i in 1:dim){
		for(j in 1:dim){
			if(abs(i-j)==1){
				prec[i,j] <- off.diag
			}else if(abs(i-j)==2){
				prec[i,j] <- off.diag^2
			}
		}
	}

return(prec)
}


################ full off-diagonal matrix ##########################

model3 <- function(dim,xi.val,off.min, off.max){

	prec <- diag(1, dim)*xi.val
	
	for(i in 1:dim){
		for(j in i:dim){
			if(i!=j){
				off.val <- runif(1,off.min,off.max)
				prec[i,j] <- off.val
				prec[j,i] <- off.val
			}
		}
	}
	return(prec)
}
i <- 0
for(nrows in nrow.cand){

i <- i+1


Z <- matrix(rnorm(n=nrows*ncols,mean=mn,sd=sdd),nrow=nrows,ncol=ncols)

Omega <- model2(dim=ncols,1,0.5)






cat("processing time for single node is\n")
 ptm1 <- proc.time()
 r_K <- Z%*%Omega%*%t(Z)
  
 print(r_K[1:5,1:5])
 
if(i==1){
	single.time <- proc.time()-ptm1

}else{

        single.time <- rbind(single.time, proc.time()-ptm1)
}
}
  save.image("~/src/my_project/parallel_project/Comp_Kernel_Prod_Single.RData")
