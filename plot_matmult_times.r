

blas <- read.csv("matmult_times_blas.csv",header=FALSE)
fast <- read.csv("matmult_times_fast.csv",header=FALSE)
slow <- read.csv("matmult_times_naive.csv",header=FALSE)
range <- 1:(length(fast))

dif = slow - fast

dif = dif / slow
dif = dif * 100
print(mean(dif))
print(mean(slow - blas))

jpeg(file="plot.jpeg",width = 1000, height = 1000)

plot(range, fast, type='l', col='red', xlab='size', ylab='t(ns)', lwd=2)
#add second line to plot
lines(range, slow, col='green',lwd=2)
lines(range, blas, col='blue',lwd=2)

dev.off()

jpeg(file="difference.jpeg",width = 1000, height = 1000)

plot(range, dif, type='l', col='blue', xlab='size', ylab='%', lwd=2)
dev.off()

print(paste0(paste0("total slow: ", sum(slow)),"ns"))
print(paste0(paste0("total fast: ", sum(fast)),"ns"))
print(paste0(paste0("total BLAS: ", sum(blas)),"ns"))
print(paste0(paste0("difference: ", ((sum(slow)-sum(fast))/sum(slow))*100), "%"))
print(paste0(paste0("difference BLAS ", ((sum(slow)-sum(blas))/sum(slow))*100), "%"))
