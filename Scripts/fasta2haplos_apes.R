
# read, create a fasta with unique haplotypes and record their frequency (in traits file)

args=commandArgs(T)

fin=args[1]
max_hap=as.numeric(args[2])
min_freq=as.numeric(args[3])
binary=as.numeric(args[4])

if (binary) {
	fout_fa=paste(fin, "_", max_hap, "_", min_freq, ".bin.haps.fasta", sep="", collapse="")
	fout_trait=paste(fin, "_", max_hap, "_", min_freq, ".bin.traits.txt", sep="", collapse="")
} else {
	fout_fa=paste(fin, "_", max_hap, "_", min_freq, ".haps.fasta", sep="", collapse="")
        fout_trait=paste(fin, "_", max_hap, "_", min_freq, ".traits.txt", sep="", collapse="")
}

cat("", file=fout_fa)
cat("", file=fout_trait)

# read haplo
ss=readLines(fin)
pops=substring(trimws(ss[seq(1,length(ss),2)]),2)
haps=trimws(ss[seq(2,length(ss),2)])
rm(ss)

uhap=names(rev(sort(table(haps)))) # sorted unique
upop=unique(pops)
cat(paste(c("",upop),sep="",collapse=","), "\n", file=fout_trait, append=T)

if (max_hap==0) max_hap=length(uhap);
if (max_hap>length(uhap)) max_hap=length(uhap);

uhap=uhap[1:max_hap]

# binary format 0/1, fake A/C
if (binary) {

	buhap=c(); for (i in 1:length(uhap)) buhap[i]=""

	for (i in 1:nchar(uhap[1])) {

		major=substring(uhap[1],i,i)
	
		if (length(unique(substring(uhap,i,i)))>2) {
			cat("\n",i); print(table(substring(uhap,i,i)))
		}

		alls=rep("A", length(uhap))
		alls[which(substring(uhap,i,i)!=major)]="C"

		for (j in 1:length(uhap)) buhap[j]=paste(buhap[j],alls[j],sep="",collapse="")
	}
}

#444
# C  H  R 
# 3  1 33 

#  1316
#  M  T  V 
#  33  2  2

# get frequencies and print
nh=0 # incremental
for (i in 1:length(uhap)) {
	freqs=c()
	ii=which(!is.na(match(haps, uhap[i])))
	for (j in 1:length(upop)) {
		freqs[j]=length(which(pops[ii]==upop[j]))
	}
	if ( (sum(freqs)>=min_freq) ) { 
		nh=nh+1
		nn=paste("hap_",i,sep="",collapse="")
		cat(paste(">",nn,sep="",collapse=""), "\n", file=fout_fa, append=T)
		if (binary) cat(buhap[i], "\n", file=fout_fa, append=T) else cat(uhap[i], "\n", file=fout_fa, append=T)
		cat(paste(c(nn,freqs),sep="",collapse=","), "\n", file=fout_trait, append=T)
	}
}

cat("Processed", nh, "haplotypes.\n")
cat("Written these files:", fout_fa, fout_trait, "\n")






