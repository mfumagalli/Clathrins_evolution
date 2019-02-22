
# these files are from snp150 versions

datadir <- "~/Data/UCSC/SNP/"

for (chrom in 1:22) {

	cat(chrom,"\n")

	snps <- read.table(paste(datadir, "snp_chr", chrom, ".gz", sep="", collapse=""), stringsAsFactors=F, header=F, sep="\t")
	cat(".")
	snps_cod <- read.table(paste(datadir, "snp_cod_chr", chrom, ".gz", sep="", collapse=""), stringsAsFactors=F, header=F, sep="\t")
	cat(".")
	snps_ortho <- read.table(paste(datadir, "snp_ortho_chr", chrom, ".gz", sep="", collapse=""), stringsAsFactors=F, header=F, sep="\t")
	cat(".")
	save(snps, snps_cod, snps_ortho, file=paste(datadir,"snps_chr", chrom, ".Rdata", sep="", collapse=""))

	rm(snps)
	rm(snps_cod)
	rm(snps_ortho)

}


