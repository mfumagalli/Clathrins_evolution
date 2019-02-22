
# Results/out_sup_CLTCL1.fasta_20_0.haps.fasta Results/out_sup_CLTCL1.fasta_20_0.traits.txt

args <- commandArgs(T)
fin <- args[1]
ffreq <- args[2]
rm(args)

fout <- paste(fin,".table.csv",sep="",collapse="")
fout2 <- paste(fin,".diff.txt",sep="",collapse="")

cat("", file=fout)
cat("", file=fout2)

ss <- readLines(fin)
haps <- trimws(ss[seq(2,length(ss),2)])
freqs <- trimws(readLines(ffreq))

pops <- substring(freqs[1], 2, nchar(freqs[1]))
cat("I found these labels:", pops, "\n")

cat(paste(c("hap_ID",pops, seq(1,nchar(haps[1])), ""), sep="", collapse=","), "\n", file=fout)

if (length(haps)>1) {

	for (i in 1:length(haps)) cat(paste(c(freqs[i+1], ",", paste(strsplit(haps[i],split="")[[1]],sep="",collapse=","), ","), sep="", collapse=""), "\n", file=fout, append=T)

	# get differences
	cat("", file=fout2)
	for (i in 1:(length(haps)-1)) {
		for (j in (i+1):length(haps)) {
			h1 <- haps[i]; h2 <- haps[j]
			diffs <- which(strsplit(h1, split="")[[1]]!=strsplit(h2, split="")[[1]]) 
			cat(i,j,paste(diffs, sep="",collapse=","), sep="\t", "\n", file=fout2, append=T)
		}
	}
}

cat("Written these files:", fout, fout2, "\n")







