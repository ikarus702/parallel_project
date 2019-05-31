suppressMessages(library(pbdMPI))
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
finalize()
