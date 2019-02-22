
# extract haplotypes for all populations and compute summary statistics
source("Scripts/popgen.rtx")

# extract haplotypes for CEU and CHB and compute stats

args <- commandArgs(T)
gene <- toupper(args[1])
length <- abs(as.numeric(args[3]) - as.numeric(args[2])) + 1
rm(args)

fout <- paste("Results/Genes/",gene,".stats.txt",sep="",collapse="")

# read haplos for all populations

f_freqs <- paste("Results/Genes/out_pop_",gene,".fasta_0_0.traits.txt",sep="",collapse="")
f_seqs <- paste("Results/Genes/out_pop_", gene,".fasta_0_0.haps.fasta",sep="",collapse="")

haps <- trimws(readLines(f_seqs)); haps=haps[seq(2,length(haps),2)]
freqs <- read.csv(f_freqs, head=T)
pops <- colnames(freqs)[-1]

# here I assume the last three are altai, deni, chimp!!!

# fill in for each pop
haplist <- list()
for (i in 1:(length(pops)-3)) {
	haplist[[i]] <- rep(haps, freqs[,(i+1)])
}

## summary stats
cat(c("POP","NSAM","LEN","S","SS","TW","PI","TD","FLDs","FLFs","H1","H2","H12","H2/H1","H2/H1''","FST_YRI","FST_CEU","FST_CHB"), sep="\t", "\n", file=fout)

yri <- which(pops=="YRI"); ceu <- which(pops=="CEU"); chb <- which(pops=="CHB");
for (i in 1:length(haplist)) {

	nsam <- length(haplist[[i]])
        taj <- tajima(haplist[[i]])[3:5] # TW PI TD
        fu <- fuli(haplist[[i]]) # 1 and 2 are S and SS, 3 and 4 are Ds and Fs
        hs <- homohapl(haplist[[i]]) # H1, H2, H2/H1
        yri <- which(pops=="YRI"); ceu <- which(pops=="CEU"); chb <- which(pops=="CHB"); gi <- which(pops=="GI")
        fst1 <- reynolds(list(haplist[[i]], haplist[[yri]]))
        fst2 <- reynolds(list(haplist[[i]], haplist[[ceu]]))
        fst3 <- reynolds(list(haplist[[i]], haplist[[chb]]))

        res <- c(pops[i], nsam, length, fu[1], fu[2], taj, fu[3:4], hs, fst1, fst2, fst3)
        cat(res, sep="\t", "\n", file=fout, append=T)

}

cat("Written this file:", fout, "\n")

