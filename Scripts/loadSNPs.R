
for (chrom in c(17,22)) {

cat("\n",chrom)

snp144=paste("Data/snp144_chr", chrom, ".txt", sep="", collapse="")
snp144_coding=paste("Data/snp144_coding_chr", chrom, ".txt", sep="", collapse="")
snp144_ortho=paste("Data/snp144_ortho_chr", chrom, ".txt", sep="", collapse="")

## select only rs for coding snps, and annotated function, and select chimp allele
snps=      read.table(snp144, stringsAsFactors=F, head=F, sep="\t")
snps_cod=  read.table(snp144_coding, stringsAsFactors=F, head=F, sep="\t")
snps_ortho=read.table(snp144_ortho, stringsAsFactors=F, head=F, sep="\t")

save(snps, snps_cod, snps_ortho, file=paste("Data/snps_",chrom,".Rdata",sep="",collapse=""))

}

