
# from the null distruibution, compute percentile ranks for all tests for CLTC CLTCL1

library(gplots)

gene <- commandArgs(T)

# stats for the gene of interest
stat <- read.table(paste("Results/Genes/", gene, ".stats.txt", sep="", collapse=""), head=T)[,-2]
stat$TW <- stat$TW/stat$LEN
stat$PI <- stat$PI/stat$LEN

ll <- system("ls Results/Genes/*.stats.txt", intern=T)
# remove clathrin
ll <- setdiff(ll, c("Results/Genes/CLTCL1.stats.txt","Results/Genes/CLTC.stats.txt", "Results/Genes/CLTA.stats.txt", "Results/Genes/CLTB.stats.txt"))
# remove this gene family which seems unusually polymorphic
#ll <- setdiff(ll, c("Results/Genes/PCDHA2.stats.txt","Results/Genes/PCDHA3.stats.txt", "Results/Genes/PCDHA4.stats.txt", "Results/Genes/PCDHA5.stats.txt", "Results/Genes/PCDHA6.stats.txt", "Results/Genes/PCDHA7.stats.txt", "Results/Genes/PCDHA8.stats.txt"))

cat("\nFound",length(ll),"control genes.")

# read all files, remove unwanted pops or statistics
genenames <- c()
ctr <- c()
for (i in ll) {
	tmp <- read.table(i, head=T)
	genenames <- c(genenames, strsplit(strsplit(i, split="Results/Genes/")[[1]][2], split=".stats.txt")[[1]][1])
	ctr <- rbind(ctr, tmp[,-2]) # remove sample size
}
ctr$TW <- ctr$TW/ctr$LEN
ctr$PI <- ctr$PI/ctr$LEN

# check if all different (because of the cutoff on length done later on)
if(length(unique(ctr$LEN))!=length(ll)) stop("ll diff than len")
minmaxlen <- range(unique(ctr$LEN)[sort(abs(unique(ctr$LEN)-stat$LEN[1]), ind=T)$ix[1:500]])
# these are 500 closest genes even in genomic length

# initialise results
res <- stat[,-1] # no popname
pvs <- res
nvals <- res

# for each pop
for (i in 1:nrow(res)) {
	# row index for pop i
	ind <- which(as.character(ctr[,1])==as.character(stat[i,1]))
	if (length(ind)!=length(ll)) stop("missing pops?")
	# but remove controls where S<2 and length is outside 
	ind <- ind[which(ctr$S[ind]>1 & ctr$LEN[ind]>=minmaxlen[1] & ctr$LEN[ind]<=minmaxlen[2])]
	# for each stat
	for (j in 1:ncol(res)) {
		val <- c(stat[i,j+1], ctr[ind,j+1])
		if (!is.na(val[1])) {
			val <- val[which(!is.na(val))]
			pr <- rank(val, ties.method="average", na.last=NA)/length(val)
			res[i,j] <- pr[1]
			nvals[i,j] <- length(val)
			pvs[i,j] <- ""
			toput <- round(pr[1]*100)/100
			if (toput==1) toput <- 0.99
			if (toput==0) toput <- 0.01
			if (pr[1]<=0.10 | pr[1]>=0.90) pvs[i,j] <- toput
			#if (pr[1]<0.025 | pr[1]>0.975) pvs[i,j] <- "*"
			#if (pr[1]<0.005 | pr[1]>0.995) pvs[i,j] <- "**"
			#if (pr[1]<0.0025 | pr[1]>0.9975) pvs[i,j] <- "***"
		}
	}
}

# select statistics
tokeep <- c(4,5,6:10, 13)
sres <- res[,tokeep]
colnames(sres)[ncol(sres)] <- "H2H1"
spvs <- pvs[,tokeep]

# pop names
pops <- as.character(stat[,1])
rownames(res) <- pops

labs <- read.table("Data/pop_labels.txt", sep="\t", stringsAsFact=F)
rownames(sres) <- paste(pops, " (", trimws(labs$V3[match(pops, labs[,1])]), ")", sep="" )

pdf(file=paste("Plots/", gene, ".revision.stats.pdf", sep="", collapse=""))

heatmap.2(as.matrix(sres), scale="none", main="Summary statistics", dendrogram="none", cellnote=spvs, notecol="black", notecex=1, key.xlab="Empirical rank", key.ylab="Frequency", trace="non")

dev.off()

cat("\nMedian nr of S>0:", median(nvals$S))

# fst

# YRI, CEU, CHB

# YRI-CEU
fst <- stat[which(stat[,1]=="YRI"),16]
fst_ctr <- ctr[which(ctr[,1]=="YRI"),16]
(rank(c(fst, fst_ctr), ties.method="average", na.last=NA)[1])/(length(fst_ctr)+1)

# YRI-CHB
fst <- stat[which(stat[,1]=="YRI"),17]
fst_ctr <- ctr[which(ctr[,1]=="YRI"),17]
(rank(c(fst, fst_ctr), ties.method="average", na.last=NA)[1])/(length(fst_ctr)+1)

# CEU-CHB
fst <- stat[which(stat[,1]=="CEU"),17]
fst_ctr <- ctr[which(ctr[,1]=="CEU"),17]
(rank(c(fst, fst_ctr), ties.method="average", na.last=NA)[1])/(length(fst_ctr)+1)












