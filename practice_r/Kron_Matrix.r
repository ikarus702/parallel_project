library(pbdDMAT, quiet=TRUE)
library(pryr)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
Nnode = as.numeric(args[1])

init.grid()

n <- 1e3
p <- 1e3



blm <- c(32)


if(comm.rank()==0){
r_I <- Diagonal(n)
r_onevec <- matrix(rep(1,n),ncol=1)
r_kron_I <- kronecker(r_I,r_onevec)

}else{
r_I <- NULL
r_onevec <- NULL
r_kron_I <- NULL
}
barrier()

if(comm.rank()==0){

print("memory usage before ddmat")

}
comm.print(mem_used())

############# option 1 ##################################

kron_time <- comm.timer({ 

dA <- as.ddmatrix(x=matrix(c(1:1e4),100,100),bldim=c(blm,blm))
d_onevec <- as.ddmatrix(x=r_onevec, bldim=c(blm,blm))
d_kron_I <- as.ddmatrix(x=matrix(submatrix(dA[,1]),sqrt(nrow(submatrix(dA[,1]))),sqrt(nrow(submatrix(dA[,1])))),bldim=c(blm,blm))
#tmp.vec <- rep(0,n^2)


#d_tot_mat <- as.list(1:10)
#d_kron_I <- as.list(1:100)

#for(j in 1:10){

#for(i in 1:100){

#index <- c((i-1)*n+1):c(i*n)
#tmp.vec[index] <- 1
#d_mat <- ddmatrix(tmp.vec, ncol=1, bldim=c(blm,blm))
#d_kron_I[[i]] <- d_mat

#tmp.vec <- rep(0,n^2)
#}


#d_tot_mat[[j]] <- do.call(cbind,d_kron_I)
#d_kron_I <- as.list(1:100)



#}

#d_kron_I <- do.call(cbind,d_tot_mat)

})

comm.print(class(d_kron_I))
comm.print(submatrix(d_kron_I),all.rank=TRUE)
#comm.print(submatrix(d_kron_I))

if(comm.rank()==0){
print("memory usage after ddmat")

}
comm.print(mem_used(),all.rank=TRUE)

if(comm.rank()==0){
print("kronecker time is \n")
}
comm.print(kron_time, all.rank=FALSE)






if(comm.rank()==0){


 save.image("~/src/my_project/parallel_project/Comp_Kronecker_Mat.RData")
}


finalize()
