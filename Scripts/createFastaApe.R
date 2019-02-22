
# reads the raw fasta (id) and create fasta based on pop (pop) species (spec)

gene <- toupper(commandArgs(T))

fin <- paste("Results/Apes/",gene,".fasta",sep="",collapse="")
fout1 <- paste("Results/Apes/", gene, ".id.fasta", sep="", collapse="")
fout2 <- paste("Results/Apes/", gene, ".pop.fasta", sep="", collapse="")

fa <- trimws(readLines(fin))
# remove white spaces and rewrite
# cat(fa, sep="\n", file=fin)

# id and pops
ids <- fa[seq(1,length(fa),3)]
tmp <- unlist(strsplit(ids, split="-"))
pop <- tmp[seq(1, length(tmp),2)]
rm(tmp)

fa1 <- c()
fa2 <- c()
for (i in 1:length(ids)) { 
	fa1 <- c(fa1, c(ids[i], fa[(i*3)-1], ids[i], fa[(i*3)]))
	fa2 <- c(fa2, c(pop[i], fa[(i*3)-1], pop[i], fa[(i*3)]))
}

cat(fa1, sep="\n", file=fout1)
cat(fa2, sep="\n", file=fout2)






