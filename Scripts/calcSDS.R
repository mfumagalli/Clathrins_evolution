
args <- commandArgs(T)
fin <- args[1]
fout <- args[2]
rm(args)

res <- read.table(fin, stringsAsFact=F, head=T, sep="\t")
res <- res[which(res$DAF>0.05),]

quants <- quantile(res$DAF, seq(0,1,0.05))
cat("These are the quantiles:\n", quants, ".\n")

bins <- findInterval(res$DAF, quants, rightmost.closed=T)

stdSDS <- rep(NA, nrow(res))
for (i in 1:max(bins)) {
	ind <- which(bins==i)
	stdSDS[ind] <- scale(res$rSDS[ind])[,1]
}

cat("This is the range of stdSDS:", range(stdSDS), ".\n")

write.table(cbind(res, stdSDS=stdSDS), sep="\t", file=fout, col.names=T, row.names=F, quote=F)



