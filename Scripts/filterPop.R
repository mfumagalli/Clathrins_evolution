
args=commandArgs(T)
fin=args[1]
fint=args[2]
pop=args[3]
min_freq=as.numeric(args[4])
rm(args)

fout=paste(fin,".",pop,".haps.fasta",sep="",collapse="")
foutt=paste(fint,".",pop,".txt",sep="",collapse="")

pops=read.table("/data/data/CLTCL1/integrated_call_samples_v3.20130502.ALL.panel", stringsAsFactors=F, head=T)
tokeep=c(unique(pops$pop[which(pops$super_pop==pop)]), "CHIMP")

freqs=read.csv(fint)
freqs2=freqs[, c(1, match(tokeep, colnames(freqs)))]
colnames(freqs2)[1]=""
rm(freqs)

ind=which(apply(FUN=sum, MAR=1, X=freqs2[,2:ncol(freqs2)])>=min_freq)
cat("Processed", length(ind), "haplotypes.\n")

write.csv(freqs2[ind,], file=foutt, quote=F, row.names=F)

ss=readLines(fin)
pops=substring(trimws(ss[seq(1,length(ss),2)]),2)
haps=trimws(ss[seq(2,length(ss),2)])

cat("", file=fout)
for (i in ind) {
	cat(paste(">",pops[i],sep="",collapse=""), "\n", file=fout, append=T)
        cat(haps[i], "\n", file=fout, append=T)
}



cat("Written these files:", fout, foutt, "\n")




