args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

library(bigKRLS)
library(openblasctl)



openblas_set_num_threads(1)

load(paste0("~/src/RData/NGK_PK_GP_Data_n",n_num,"_p",p_num,".RData"))

Data <- data.list[[1]]
y <- as.matrix(Data[,1])
X <- as.matrix(Data[,-1])

ptm <- as.list(1:20)

for(i in 1:20){

ptm[[i]] <- system.time(krls.summary <- bigKRLS(X=X,y=y))


}

ptm1 <- do.call(rbind,ptm)

ptm.summary <- apply(ptm1,2,mean)

print(ptm.summary)

