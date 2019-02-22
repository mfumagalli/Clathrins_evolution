
# take the 5 most common haplotypes in humans and convert it into full protein sequence
# then take the merge (?) with chimp protein (or just take the most common one?)
# finally do a network

# 1) get reference human protein 
# (this first part is taken from makeProteinChimp)

source("Scripts/functions.R")

genename <- "CLTCL1"

refgene <- read.table("Data/refGene.txt")
colnames(refgene) <- strsplit(readLines("Data/refGene.txt")[1], split="\t")[[1]]

ii <- which(refgene$name2==genename)
if (length(ii)==1) refgene <- refgene[ii,] else stop("Gene name not found or multiple hits.");
rm(ii)

strand <- as.character(refgene$strand)
chrom <- as.numeric(strsplit(as.character(refgene$chrom), split="chr")[[1]][2])
starts <- as.numeric(strsplit(as.character(refgene$exonStarts), split=",")[[1]])+1 # because they are 0-based
ends <- as.numeric(strsplit(as.character(refgene$exonEnds), split=",")[[1]])

# retrieve nucleotide of coding (bases)
bases <- c()
for (i in 1:length(starts)) {
        primo <- ultimo <- 0;
        if (refgene$cdsStart<=ends[i] & refgene$cdsEnd>=starts[i]) {
                if (refgene$cdsStart>=starts[i] & refgene$cdsStart<=ends[i]) primo <- 1 else primo <- 0;
                if (refgene$cdsEnd>=starts[i] & refgene$cdsEnd<=ends[i]) ultimo <- 1 else ultimo <- 0;
                if (primo==1) starts[i] <- refgene$cdsStart+1 # 0-based
                if (ultimo==1) ends[i] <- refgene$cdsEnd
                bases <- c(bases, seq(starts[i],ends[i]))
                #cat(i,":",length(seq(starts[i],ends[i])),"=",length(bases),"\n")
        } else cat(i,"non-coding exon skipped\n")
}
cat(length(bases)," ", length(bases)/3) # 4923   1641

# read human haplotypes
haps <- read.csv("Results/Genes/out_sup_CLTCL1.fasta_0_0.haps.fasta.table.csv", stringsAsFact=F)
# they are sorted by frequency: apply(X=haps[,2:6],MAR=1,FUN=sum); first 8 have more than 50 frequency
anno <- read.csv("Results/Genes/anno_CLTCL1.csv", head=F, stringsAsFact=F)

# get reference sequence for coding
mine <- list(); for (j in 1:8) mine[[j]] <- ""
for (i in 1:length(bases)) {
        cmd <- paste(paste("samtools faidx ~/Data/hg19/chr", chrom, ".fa.gz", sep="", collapse=""), paste(as.character(refgene$chrom),":", paste(rep(bases[i],2),sep="",collapse="-"), sep="", collapse=""))
	nucl <- system(cmd, intern=T)[2]
	ind <- which(anno$V2==bases[i])
	if (length(ind)==0) {
		for (j in 1:8) mine[[j]] <- c(mine[[j]], nucl)
	} else {
		if (length(ind)==1) {
			cat("\t",ind)
			if (nucl %in% c(anno[ind,4:5])) {
				for (j in 1:8) mine[[j]] <- c(mine[[j]], as.character(haps[j,ind+9]))
			} else stop("diff alleles in annotation?")
		} else stop("multiple lines?!")
	}	
}
for (j in 1:8) mine[[j]] <- paste(mine[[j]], sep="", collapse="")

# put reverse in case
if (strand=="-") {
        for (j in 1:8) mine[[j]] <- reverse(mine[[j]])
        bases <- rev(bases)
}

# translate
prots <- rep(NA, 8)
for (j in 1:8) prots[j] <- translate(mine[[j]])

# 2) get chimp data

gene <- "CLTCL1"
fin <- paste("Results/Chimp/",gene,".fasta",sep="",collapse="")

fa <- trimws(readLines(fin))

haps_out <- fa[c(seq(2,length(fa),3),seq(3,length(fa),3))]
df <- as.data.frame(table(haps_out))
prots_out <- as.character(df$haps[which(df$Freq>10)])

# 3) write fasta

# with AA yo do phylogenetics (does it work?)
#fout <- "Results/CLTCL1.hg-chimp.protein.fasta"
#cat("", file=fout)
#for (i in 1:length(prots)) cat("\n>Homo_sapiens\n", prots[i], file=fout, append=T)
#for (i in 1:length(prots_out)) cat("\n>Ancestral\n", prots_out[i], file=fout, append=T)

fout <- "Results/CLTCL1.hg-chimp.protein.fasta"
cat("", file=fout)

uhap <- c(prots, prots_out)
buhap <- c(); for (i in 1:length(uhap)) buhap[i] <- ""

for (i in 1:nchar(uhap[1])) {

	major <- substring(uhap[1],i,i)

	nrall <- length(unique(substring(uhap,i,i)))

	if (nrall>2) cat("\n!!!",i,":",substring(uhap,i,i))

	if (nrall==2) {
		cat("\n",i,  nrall)
        	alls <- rep("A", length(uhap))
        	alls[which(substring(uhap,i,i)!=major)] <- "G"
        	for (j in 1:length(uhap)) buhap[j] <- paste(buhap[j], alls[j], sep="", collapse="")
	}
}

for (i in 1:length(prots)) cat("\n>Homo sapiens\n", buhap[i], file=fout, append=T)
for (i in (length(prots)+1):(length(prots)+length(prots_out))) cat("\n>Ancestral\n", buhap[i], file=fout, append=T)

# 4) plot net

library(pegas)

fin <- fout

d <- ape::read.dna(fin, format="fasta")
e <- dist.dna(d)

h <- pegas::haplotype(d)
h <- sort(h, what = "label")

net <- pegas::haploNet(h)
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
        table(hap=ind, pop=rownames(d)[values])
)

human_size <- round(apply(X=haps[,2:6],MAR=1,FUN=sum)[1:8] / sum(apply(X=haps[,2:6],MAR=1,FUN=sum)[1:8]) * 100)
chimp_size <- round(df$Freq[which(df$Freq>10)]/sum(df$Freq[which(df$Freq>10)]) * 100)

pdf(file="Plots/network.ancestral.pdf")

plot(net, size=c(human_size,chimp_size), scale.ratio=40, pie=ind.hap, labels=F, legend=F, show.mutation=1, th=0, main="CLTCL1")
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)

dev.off()











