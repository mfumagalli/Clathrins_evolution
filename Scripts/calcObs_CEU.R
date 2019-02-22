
source("Scripts/popgen.rtx")

# extract haplotypes for CEU and CHB and compute stats

haplist=list()

# CEU

f_freqs="Results/out_pop_CLTCL1.fasta_0_0.traits.txt.EUR.txt"
f_seqs="Results/out_pop_CLTCL1.fasta_0_0.haps.fasta.EUR.haps.fasta"

haps=trimws(readLines(f_seqs)); haps=haps[seq(2,length(haps),2)]
freqs=read.csv(f_freqs, head=T)$CEU

haplist[[1]]=rep(haps, freqs)

# CHB

f_seqs="Results/out_pop_CLTCL1.fasta_0_0.haps.fasta.EAS.haps.fasta"
f_freqs="Results/out_pop_CLTCL1.fasta_0_0.traits.txt.EAS.txt"

haps=trimws(readLines(f_seqs)); haps=haps[seq(2,length(haps),2)]
freqs=read.csv(f_freqs, head=T)$CHB

haplist[[2]]=rep(haps, freqs)


## CHECK WHY THESE SAMPLES SIZES ARE DIFFERENT FROM COUNTS FROM HETERO !!!!!!!!!!!!!!!

## SUMMARY STATS on CEU
taj=tajima(haplist[[1]])[5] #TD
fu=fuli(haplist[[1]]) # 1 and 2 are S and SS, 3 and 4 are Ds and Fs
hs=homohapl(haplist[[1]]) # H1, H2, H2/H1
fsts=reynolds(haplist)
outres=c(fu[1:2],taj,fu[3:4],hs,fsts)

#cat(outres,"\n", file=fout, append=T)

sims=read.table("summary_stats_NEU.txt.gz", head=T, stringsAsFactors=F)
rbind(outres,apply(FUN=quantile, X=sims, MAR=2, probs=c(0.05,0.50,0.95)))
names(outres)=colnames(sims[1:9])

par(mfrow=c(1,2))
hist(sims$S_20, breaks=20, freq=T, main="No. SNPs", xlab=""); abline(v=outres[1], lty=2)
hist(sims$H2_20, breaks=20, freq=T, main="H2", xlab="", sub="H1=sum(hfreqs^2);H2=H1-(hfreqs[1]^2)"); abline(v=outres[7], lty=2)

## temp

sims=read.table("summary_stats_BAL.txt", head=T, stringsAsFactors=F)
pars=read.table("params_BAL.txt", head=T, stringsAsFactors=F)[1:nrow(sims),5]
#rbind(outres,apply(FUN=quantile, X=sims, MAR=2, probs=c(0.05,0.50,0.95)))
names(outres)=colnames(sims[1:9])

rr=sort(sqrt((sims$H2_5-outres[7])^2), ind=T)


par(mfrow=c(1,2))
hist(sims$S_20, breaks=20, freq=T, main="No. SNPs", xlab=""); abline(v=outres[1], lty=2)
hist(sims$H2_20, breaks=20, freq=T, main="H2", xlab="", sub="H1=sum(hfreqs^2);H2=H1-(hfreqs[1]^2)"); abline(v=outres[7], lty=2)


## MAKE SOME PREDICTIONS ON LIKELY SEL STRENGTH?


