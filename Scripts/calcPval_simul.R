
# from the simulated distruibution, compute percentile ranks for all tests for CLTC CLTCL1

library(gplots)

gene <- commandArgs(T)

pops <- c("YRI","CEU","CHB")

# stats for the gene of interest
stats <- read.table(paste("Results/Genes_allsites/", gene, ".stats.txt", sep="", collapse=""), head=T)[,-2]
# remove stuff (do the same for all genes, see below)
stats <- stats[which(stats[,1]!="GI"),]
stats <- stats[,which(colnames(stats)!="FST_GI")]

# read simulations
fsim <- "msms.CLTCL1.stats.txt"
sims <- read.table(fsim, stringsAs=F, sep=" ", head=F)[,-2]
colnames(sims) <- colnames(stats)

# for each pop
for (j in 1:length(pops)) {

	ind1 <- which(stats[,1]==pops[j])
	ind2 <- which(sims[,1]==pops[j])

	data <- as.matrix(rbind(stats[ind1,-1], sims[ind2,-1]))

	apply(FUN=rank, MAR=2, X=data, ties.method="random")/nrow(data)

}



# initialise results
res <- stat[,-1]
pvs <- res
nvals <- res

# for each pop
for (i in 1:nrow(res)) {
	# row index for pop i
	ind <- which(as.character(ctr[,1])==as.character(stat[i,1]))
	# but remove controls where S=0
	ind <- ind[which(ctr[ind,2]>0)]
	# for each stat
	for (j in 1:ncol(res)) {
		val <- c(stat[i,j+1], ctr[ind,j+1])
		if (is.na(val[1])) stop("stat is NA.")
		val <- val[which(!is.na(val))]
		pr <- rank(val, ties.method="random")/length(val)
		res[i,j] <- pr[1]
		nvals[i,j] <- length(val)
		pvs[i,j] <- ""
		if (pr[1]<0.05 | pr[1]>0.95) pvs[i,j] <- "*"
		if (pr[1]<0.01 | pr[1]>0.99) pvs[i,j] <- "**"
		if (pr[1]<0.005 | pr[1]>0.995) pvs[i,j] <- "***"
	}
}

# select statistics
tokeep <- c(1,2,3,4,5,7,8,11)
sres <- res[,tokeep]
pvs <- pvs[,tokeep]

# pop names
pops <- as.character(stat[,1])
rownames(res) <- pops

labs <- read.table("Data/pop_labels.txt", sep="\t", stringsAsFact=F)
rownames(sres) <- paste(pops, " (", trimws(labs$V3[match(pops, labs[,1])]), ")", sep="" )

cat("\nProducing plot...")

pdf(file=paste("Plots/", gene, ".revision.stats.pdf", sep="", collapse=""))

#heatmap.2(as.matrix(sres), scale="none", main="CLTCL1", dendrogram="both", cellnote=pvs, notecol="black", notecex=1.5, key.xlab="empirical rank", key.ylab="Frequency", trace="non")
heatmap.2(as.matrix(sres), scale="none", main="CLTCL1", dendrogram="none", cellnote=pvs, notecol="black", notecex=1.5, key.xlab="empirical rank", key.ylab="Frequency", trace="non")

dev.off()

cat("\nMedian nr of S>0:", median(nvals$S))

# fst

# YRI, CEU, CHB, GI

fst <- as.numeric(c(stat[which(pops=="YRI"),12:14], stat[which(pops=="CEU"),13:14], stat[which(pops=="CHB"),14]))
fst_ctr <- cbind(ctr[which(ctr[,1]=="YRI"),12:14], ctr[which(ctr[,1]=="CEU"),13:14], ctr[which(ctr[,1]=="CHB"),14])

pp <- c("YRI-CEU","YRI-CHB","YRI-GI","CEU-CHB","CEU-GI","CHB-GI")

#xl <- c(0,max(as.numeric(c(unlist(fst_ctr), fst)), na.rm=T)+0.05 )

#par(mfrow=c(2,3))
cat("\nFST empirical ranks:")
for (i in 1:length(fst)) {

	all <- fst_ctr[,i]
	all <- all[which(!is.na(all))]
	all[which(all<0)] <- 0

	cat("\n",pp[i],":",length(which(all>fst[i]))/length(which(!is.na(all))))

	#hist(all, xlab="FST", ylab="Counts", main=pp[i], xlim=xl)
	#abline(v=fst[i], lty=2)
}













