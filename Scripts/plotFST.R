
# plot heterozygosity vs fst

# look at control genes
ll <- system("ls Results/Genes/*.stats.txt", intern=T)
# remove other clathrins and strange gene family ?
ll <- setdiff(ll, c("Results/Genes/CLTC.stats.txt", "Results/Genes/CLTA.stats.txt", "Results/Genes/CLTB.stats.txt"))# "Results/Genes/PCDHA2.stats.txt", "Results/Genes_allsites/PCDHA3.stats.txt", "Results/Genes_allsites/PCDHA4.stats.txt", "Results/Genes_allsites/PCDHA5.stats.txt", "Results/Genes_allsites/PCDHA6.stats.txt", "Results/Genes_allsites/PCDHA7.stats.txt", "Results/Genes_allsites/PCDHA8.stats.txt"))

# read all files, keep only pop, length, pi, fst
ctr <- c()
for (i in ll) {
	name <- rep(strsplit(strsplit(i, split="/")[[1]][3], split=".s")[[1]][1],3)
        tmp <- read.table(i, head=T)
	ind <- which(tmp$POP %in% c("YRI","CEU","CHB"))
        ctr <- rbind(ctr, cbind(name,tmp[ind,]))
}

# calculate average het
het <- c()
for (i in unique(ctr$name)) {
	ind <- which(ctr$name==i)
	het <- c(het, rep(weighted.mean(ctr$PI[ind],w=ctr$NSAM[ind]),3)/ctr$LEN[ind[1]]) 
}
ctr <- cbind(ctr, het)

# for each comparison

par(mfrow=c(1,3)) 

ind <- which(ctr$POP=="YRI")
it <- ind[which(ctr$name[ind]=="CLTCL1")]
plot(ctr$het[ind], ctr$FST_CEU[ind]); points(ctr$het[it], ctr$FST_CEU[it], col="red", pch=16)

ind <- which(ctr$POP=="YRI")
it <- ind[which(ctr$name[ind]=="CLTCL1")]
plot(ctr$het[ind], ctr$FST_CHB[ind]); points(ctr$het[it], ctr$FST_CHB[it], col="red", pch=16)

ind <- which(ctr$POP=="CEU")
it <- ind[which(ctr$name[ind]=="CLTCL1")]
plot(ctr$het[ind], ctr$FST_CHB[ind]); points(ctr$het[it], ctr$FST_CHB[it], col="red", pch=16)

# countour plot

#ind <- which(ctr$POP=="YRI")
#it <- ind[which(ctr$name[ind]=="CLTCL1")]

#quants <- quantile(ctr$het[ind], seq(0,1,0.10))
#ints <- findInterval(ctr$het[ind], quants, all.inside=T)
#zmat <- matrix(NA, nrows=length(ind), ncols=length(ind))
#xval <- rep(NA, length(ind))
#for (i in 1:max(ints)) {
#	ii <- which(ints==i)
#	xval[ii] <- rank(ctr$FST_CHB[ind][ii], ties="random")/length(ii)
#}







