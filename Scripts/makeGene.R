
## this is to reconstruct the gene from vcf2fasta results and the sequence given by Marine

source("Scripts/functions.R")

# read Marine's sequence, remember that this gene is supposed to be reverse-stranded

marine=paste(readLines("CHC22-FL-Val.txt")[-1], sep="", collapse="")
cat("Length gene:", nchar(marine), "bp\n")
cat("Length protein:", nchar(marine)/3,"aa\n")

# read gene annotation file

refgene=read.table("refGene.txt")
colnames(refgene)=strsplit(readLines("refGene.txt")[1], split="\t")[[1]]

# retrieve table with coding exons (from cds start to end, in reversed order)

starts=as.numeric(strsplit(as.character(refgene$exonStarts), split=",")[[1]])+1 # because they are 0-based
ends=as.numeric(strsplit(as.character(refgene$exonEnds), split=",")[[1]])

bases=c()
for (i in 1:length(starts)) {
	primo=ultimo=0;
	if (refgene$cdsStart<=ends[i] & refgene$cdsEnd>=starts[i]) {
		if (refgene$cdsStart>=starts[i] & refgene$cdsStart<=ends[i]) primo=1 else primo=0;
		if (refgene$cdsEnd>=starts[i] & refgene$cdsEnd<=ends[i]) ultimo=1 else ultimo=0;
		if (primo==1) starts[i]=refgene$cdsStart+1 # 0-based
		if (ultimo==1) ends[i]=refgene$cdsEnd
		bases=c(bases, seq(starts[i],ends[i]))
		#cat(i,":",length(seq(starts[i],ends[i])),"=",length(bases),"\n")
	} else cat(i,"skipped\n")
}
cat(length(bases)," ", length(bases)/3)

# check ATG
cmd=paste("samtools faidx /data/data/CLTCL1/hs37d5.fa.gz", paste(as.numeric(strsplit(as.character(refgene$chrom), split="chr")[[1]][2]),":",paste(bases[c(1,3)],sep="",collapse="-"),sep="",collapse=""))
atg=system(cmd, intern=T)[2]
# retrieve sequence
mine=c()
for (i in 1:length(bases)) {
	cmd=paste("samtools faidx /data/data/CLTCL1/hs37d5.fa.gz", paste(as.numeric(strsplit(as.character(refgene$chrom), split="chr")[[1]][2]),":",paste(rep(bases[i],2),sep="",collapse="-"),sep="",collapse=""))
	mine=c(mine, system(cmd, intern=T)[2])
}
mine=paste(mine,sep="",collapse="")

if (atg!="ATG") mine=reverse(mine)

diff_genomic=seqdiff(marine, mine)
diff_protein=protdiff(translate(marine), translate(mine))

#diff_genomic
#   pos   aa base1 base2
#1 3189 1063     G     A
#2 3457 1153     T     C
#3 3459 1153     G     A
#4 3462 1154     G     T
#5 3946 1316     G     A
#6 4523 1508     C     G
#7 4922 1641     A     G
#diff_protein
#   pos aa1 aa2
#1 1316   V   M
#2 1508   T   R

# extract the SNP position for each haplotype and check marine's protein
anno=read.csv("Results/anno_CLTCL1.csv", head=F)
ind_gen=match(anno[,2], rev(bases))

haps=readLines("Results/out_sup_CLTCL1.fasta_20_0.haps.fasta")
haps=trimws(haps[seq(2,length(haps),2)])

# extract marine sequence from haps/anno
iif (atg!="ATG") marine_target=reverse(paste(rev(strsplit(marine,split="")[[1]][ind_gen]),sep="",collapse="")) else cat("check! case non-reverse!\nstop!")


# random...

#1316
#rev(bases)[ind_gen[9]]
#[1] 19184095

#1508?
#1508-1316+1
#[1] 193

#rev(bases)[ind_gen[9]-193]
#[1] 19188852























