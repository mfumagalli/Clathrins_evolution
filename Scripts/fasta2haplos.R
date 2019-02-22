
# read, create a fasta with unique haplotypes and record their frequency (in traits file)

args <- commandArgs(T)
fin <- args[1]
max_hap <- as.numeric(args[2])
min_freq <- as.numeric(args[3])
rm(args)

# however the last 3 sequences are altai, deni and chimp which I want to keep

fout_fa <- paste(fin, "_", max_hap, "_", min_freq, ".haps.fasta", sep="", collapse="")
fout_trait <- paste(fin, "_", max_hap, "_", min_freq, ".traits.txt", sep="", collapse="")

cat("", file=fout_fa)
cat("", file=fout_trait)

# read haplotypes
ss <- readLines(fin)
pops <- substring(trimws(ss[seq(1,length(ss),2)]),2)
haps <- trimws(ss[seq(2,length(ss),2)])
rm(ss)

chimp <- haps[which(pops=="CHIMP")]
altai <- haps[which(pops=="ALTAI")]
deni <- haps[which(pops=="DENISOVA")]

uhap <- names(rev(sort(table(haps)))) # sorted unique
upop <- unique(pops)
cat(paste(c("",upop),sep="",collapse=","), "\n", file=fout_trait, append=T)

if (max_hap==0) max_hap <- length(uhap);
if (max_hap>length(uhap)) max_hap <- length(uhap);

uhap <- uhap[1:max_hap]
if (length(altai)>0) {
	if (!(altai %in% uhap)) uhap <- c(uhap, altai)
}
if (length(deni)>0) {
	if (!(deni %in% uhap)) uhap <- c(uhap, deni)
}
if (length(chimp)>0) {
	if (!(chimp %in% uhap)) uhap <- c(uhap, chimp)
}

# get frequencies and print
nh <- 0 # incremental
for (i in 1:length(uhap)) {
	freqs <- c()
	ii <- which(!is.na(match(haps, uhap[i])))
	for (j in 1:length(upop)) {
		freqs[j] <- length(which(pops[ii]==upop[j]))
	}
	if ( (sum(freqs)>=min_freq) | (sum(freqs[6:8])>0) ) { # the second condition forces print if at least on ancient/chimp has the haplotype
		nh <- nh+1
		nn <- paste("hap_",i,sep="",collapse="")
		cat(paste(">",nn,sep="",collapse=""), "\n", file=fout_fa, append=T)
		cat(uhap[i], "\n", file=fout_fa, append=T)
		cat(paste(c(nn,freqs),sep="",collapse=","), "\n", file=fout_trait, append=T)
	}
}

cat("Processed", nh, "haplotypes.\n")
cat("Written these files:", fout_fa, fout_trait, "\n")



