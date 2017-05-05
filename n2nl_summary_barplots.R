#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

s <- read.table(args[1])
outf <- args[2]

png(outf, height=1000, width=3000)

par(mfrow=c(3, 6))
med <- s[seq(1,dim(s)[1],3),]
mean <- s[seq(2,dim(s)[1],3),]
var <- s[seq(3,dim(s)[1],3),]

for (i in 1:6) {
    if (i == 1) {
        barplot(as.matrix(med[i,][3:dim(med)[2]]), ylim=c(0,1), xaxt="n", ylab="median", cex.axis=2)
    }
    else {
        barplot(as.matrix(med[i,][3:dim(med)[2]]), ylim=c(0,1), axes=F, axisnames=F)
    }
}

for (i in 1:6) {
    if (i == 1) {
        barplot(as.matrix(mean[i,][3:dim(med)[2]]), ylim=c(0,1), xaxt="n", ylab="mean", cex.axis=2)
    }
    else {
        barplot(as.matrix(mean[i,][3:dim(med)[2]]), ylim=c(0,1), axes=F, axisnames=F)
    }
}

for (i in 1:6) {
    if (i == 1) {
        barplot(as.matrix(var[i,][3:dim(med)[2]]), ylim=c(0,0.1), xaxt="n", ylab="variance", cex.axis=2)
    }
    else {
        barplot(as.matrix(var[i,][3:dim(med)[2]]), ylim=c(0,0.1), axes=F, axisnames=F)
    }
}

graphics.off()
