
# http://stackoverflow.com/questions/25755930/how-to-plot-pie-charts-in-haplonet-haplotype-networks-pegas
# but be careful there is a bug, need to sort haplotypes
# this correct: http://stackoverflow.com/questions/31220586/r-how-to-plot-correct-pie-charts-in-haplonet-haplotyp-networks-pegas-ape-a

library(pegas)

# CLTC
fin <- "Results/Genes_allsites/out_sup_CLTC.fasta"

d <- ape::read.dna(fin, format="fasta")
e <- dist.dna(d)

h <- pegas::haplotype(d)
h <- sort(h, what = "label")

net <- pegas::haploNet(h)
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
	table(hap=ind, pop=rownames(d)[values])
)

plot(net, size=10*log2(attr(net, "freq")), scale.ratio=10, pie=ind.hap, labels=F, legend=F, show.mutation=0, th=0, main="A")
mtext("CLTC",3,0)
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)


# CLTCL1
fin <- "Results/Genes_allsites/out_sup_CLTCL1.fasta"

res <- readLines(fin)
df <- as.data.frame(table(res[seq(2,length(res),2)]))

# select haplotype with freq>10
ind <- which(res[seq(2,length(res),2)] %in% as.character(df[df$Freq>10,1]))
# but add archaic
ind <- unique(c(ind, (length(res)/2)-1, (length(res)/2)))

cat("", file="tmp.fasta")
heads <- res[seq(1,length(res),2)]
haps <- res[seq(2,length(res),2)]
for (i in ind) cat(heads[i], "\n", haps[i], "\n", sep="", collapse="", file="tmp.fasta", append=T)

d <- ape::read.dna("tmp.fasta", format="fasta")
e <- dist.dna(d)

h <- pegas::haplotype(d)
h <- sort(h, what = "label")

net <- pegas::haploNet(h)
ind.hap<-with(
	stack(setNames(attr(h, "index"), rownames(h))),
        table(hap=ind, pop=rownames(d)[values])
)

ss <- (attr(net, "freq"))
#ss[1:2] <- ss[1:2]/3

plot(net, size=ss, scale.ratio=1, pie=ind.hap, labels=F, legend=F, show.mutation=0, th=0, fast=F, main="B")

mtext("CLTCL1",3,0)
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)




# chimp


	fin <- "Results/Chimp/CLTCL1.pop.fasta_0_0.bin.haps.fasta"
	res <- readLines(fin)
	heads <- res[seq(1,length(res),2)]
        haps <- res[seq(2,length(res),2)]
	counts <- read.csv("Results/Chimp/CLTCL1.pop.fasta_0_0.bin.traits.txt")

#	counts <- counts[which(apply(counts[,-1], 1, sum)>2),]

	# THIS WORKS ONLY IF FREQUENCIES ARE SORTED AND DECREASE!!!

	cat("", file="tmp.chimp.fasta")
	for (i in 1:nrow(counts)) {
		for (j in 2:ncol(counts)) {
			volte <- counts[i,j]
			if (volte>0) {
				for (t in 1:volte) cat(">",as.character(colnames(counts)[j]),"\n",haps[i],"\n",sep="",collapse="",append=T, file="tmp.chimp.fasta")
			}
		}
	}

	fin = "tmp.chimp.fasta"

	d <- ape::read.dna(fin, format="fasta")

        e <- dist.dna(d)

        h <- pegas::haplotype(d)
        h <- sort(h, what = "label")

        net <- pegas::haploNet(h)
        ind.hap<-with(
                stack(setNames(attr(h, "index"), rownames(h))),
                table(hap=ind, pop=rownames(d)[values])
        )

	ss <- attr(net, "freq")
	ss[ss<3] <- 6

        plot(net, size=ss, scale.ratio=50, pie=ind.hap, labels=F, legend=F, show.mutation=0, th=0)
        mtext("CLTCL1",3,0)
        legend("bottomleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)









