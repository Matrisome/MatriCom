# input and output data folders
dataraw <- system.file("extdata", package = "MatriCom")
datarda <- system.file("data", package = "MatriCom")

intlist <- readRDS(paste0(dataraw, "/", "lr2.RDS"))
save(intlist, file = paste0(datarda, "/", "intlist.rda"))

mlist <- readRDS(paste0(dataraw, "/", "matrisome.list.RDS"))[[1]]
save(mlist, file = paste0(datarda, "/", "mlist.rda"))

excl <- as.data.frame(fread(paste0(dataraw, "/", "multimers.txt"), sep="\t"))
save(excl, file = paste0(datarda, "/", "excl.rda"))

signs <- as.data.frame(fread(paste0(dataraw, "/", "naba-1.csv")))
save(signs, file = paste0(datarda, "/", "signs.rda"))

CCgenes2 <- readRDS(paste0(dataraw, "/", "CCgenes2.RDS"))
save(CCgenes2, file = paste0(datarda, "/", "CCgenes2.rda"))
