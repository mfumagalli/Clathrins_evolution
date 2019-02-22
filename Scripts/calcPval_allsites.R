
# from the null distruibution, compute percentile ranks for all tests for CLTC CLTCL1

library(gplots)

gene <- commandArgs(T)

# stats for the gene of interest
stat <- read.table(paste("Results/Genes_allsites/", gene, ".stats.txt", sep="", collapse=""), head=T)[,-2]
stat$TW <- stat$TW/stat$LEN
stat$PI <- stat$PI/stat$LEN

ll <- system("ls Results/Genes_allsites/*.stats.txt", intern=T)
# remove clathrin
ll <- setdiff(ll, c("Results/Genes/CLTCL1.stats.txt","Results/Genes/CLTC.stats.txt", "Results/Genes/CLTA.stats.txt", "Results/Genes/CLTB.stats.txt", "Results/Genes/PCDHA2.stats.txt","Results/Genes/PCDHA3.stats.txt", "Results/Genes/PCDHA4.stats.txt", "Results/Genes/PCDHA5.stats.txt", "Results/Genes/PCDHA6.stats.txt", "Results/Genes/PCDHA7.stats.txt", "Results/Genes/PCDHA8.stats.txt"))
cat("\nFound",length(ll),"control genes.")

# read all files, remove unwanted pops or statistics
ctr <- c()
for (i in ll) {
	tmp <- read.table(i, head=T)
	# remove sample size
	ctr <- rbind(ctr, tmp[,-2])
}
ctr$TW <- ctr$TW/ctr$LEN
ctr$PI <- ctr$PI/ctr$LEN

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
		if (!is.na(val[1])) {
			val <- val[which(!is.na(val))]
			pr <- rank(val, ties.method="random")/length(val)
			res[i,j] <- pr[1]
			nvals[i,j] <- length(val)
			pvs[i,j] <- ""
			if (pr[1]<0.05 | pr[1]>0.95) pvs[i,j] <- "*"
			if (pr[1]<0.025 | pr[1]>0.975) pvs[i,j] <- "**"
			if (pr[1]<0.005 | pr[1]>0.995) pvs[i,j] <- "***"
		}
	}
}

# select statistics
#tokeep <- c(1,2,3,4,5,7,8,11)
tokeep <- c(2:13)
sres <- res[,tokeep]
pvs <- pvs[,tokeep]

# pop names
pops <- as.character(stat[,1])
rownames(res) <- pops

labs <- read.table("Data/pop_labels.txt", sep="\t", stringsAsFact=F)
rownames(sres) <- paste(pops, " (", trimws(labs$V3[match(pops, labs[,1])]), ")", sep="" )

cat("\nProducing plot...")

pdf(file=paste("Plots/", gene, ".all_sites.stats.pdf", sep="", collapse=""))

#heatmap.2(as.matrix(sres), scale="none", main="CLTCL1", dendrogram="both", cellnote=pvs, notecol="black", notecex=1.5, key.xlab="empirical rank", key.ylab="Frequency", trace="non")
heatmap.2(as.matrix(sres), scale="none", main="CLTCL1", dendrogram="none", cellnote=pvs, notecol="black", notecex=1.5, key.xlab="empirical rank", key.ylab="Frequency", trace="non")

dev.off()

cat("\nMedian nr of S>0:", median(nvals$S))

if(0) {

# fst

# YRI, CEU, CHB, GI

fst <- as.numeric(c(stat[which(pops=="YRI"),12:14], stat[which(pops=="CEU"),13:14], stat[which(pops=="CHB"),14]))
fst_ctr <- cbind(ctr[which(ctr[,1]=="YRI"),12:14], ctr[which(ctr[,1]=="CEU"),13:14], ctr[which(ctr[,1]=="CHB"),14])

pp <- c("YRI-CEU","YRI-CHB","CEU-CHB")

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

}











