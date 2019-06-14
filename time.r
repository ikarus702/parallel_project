suppressMessages(library(pbdMPI))
comm.set.seed(1234, diff=TRUE)

t.local = runif(1, min=.5, max=2)
comm.print(t.local, all.rank=TRUE)
barrier()

t = comm.timer(Sys.sleep(t.local))
comm.print(t)

finalize()
