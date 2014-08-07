#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

data <- read.table(args[1], row.names=1)
outf <- args[2]

png(outf, height=6000, width=6000)

par(mfrow=c(dim(data)[1]/6, 6))

for (i in 1:dim(data)[1]) {
	if (i %% 6 == 1) {
		barplot(as.matrix(data[i,]), ylim=c(0,1), xaxt="n", ylab=strsplit(rownames(data)[1], split="_")[[1]][1])
	}
	else {
		barplot(as.matrix(data[i,]), ylim=c(0,1), axes=F, axisnames=F)
	}
}

dev.off()
