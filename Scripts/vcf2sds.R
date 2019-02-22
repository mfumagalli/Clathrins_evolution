
# this script reads a VCF file and produces files for SDS

args <- commandArgs(T)
fin <- args[1]
fout_testsnps <- args[2]
fout_singletons <- args[3]
rm(args)

# fin="/data/data/CLTCL1/Genes/CHROM22GBR.vcf"
#fout_testsnps="chrom22.testsnp"
#fout_singletons="chrom22.singletons"

source("Scripts/functions.R")

# read VCF file
vcf <- read.table(fin, stringsAsFactors=F, sep="\t")
cat("I found", nrow(vcf), "sites.\n")
ninds <- ncol(vcf)-9
cat("I found", ninds, "samples.\n")
chrom <- vcf[1,1]
cat("I am using only", chrom, "chromosome out of", unique(vcf[,1]), "total.\n")

# read orthologous data
ortho <- read.table(paste("/data/data/hg19/SNP/snp147OrthoPt4Pa2Rm3.chr",chrom,".gz", sep="", collapse=""), stringsAsFactors=F, head=F, sep="\t")

# prepare output files
cat("", file=fout_testsnps)
cat("", file=fout_singletons)
singletons <- rep("", ninds)

vcfOffset <- 10

# cycle over sites
for (i in 1:nrow(vcf)) {

	if ((i %% 1e4)==0) cat(i,"/",nrow(vcf),"\n")

	id <- vcf[i,3]
	pos <- vcf[i,2]
	ref <- vcf[i,4]
	alt <- vcf[i,5]

	# if biallelic
	if (nchar(ref)==1 & nchar(alt)==1) {

		bases <- paste(sort(c(ref,alt)),sep="",collapse="/")

		# chimp allele
		chimp <- ortho[which(ortho$V5==id),]
		# forward or reverse?
		anc <- der <- "N"
		if (nrow(chimp)>0) {
			if (bases==chimp$V6) {
				anc <- chimp$V12
			}
			if (bases==reverse(chimp$V6)) {
				anc <- reverse(chimp$V12)
			}
		}
		if (!(anc %in% c(ref,alt))) anc <- ref
		if (anc==ref) {
			der <- alt
		} else {
			der <- ref
		}

		# assign genotypes and fill in missing data
		genos <- vcf[i,(vcfOffset:(vcfOffset+ninds-1))]
		genos[which(genos=="0|0")] <- 0
		genos[which(genos=="0|1" | genos=="1|0")] <- 1
		genos[which(genos=="1|1")] <- 2
		freqs <- c(length(which(genos==0)), length(which(genos==1)), length(which(genos==1)))
		genos[which(genos==".")] <- which.max(freqs)-1
	
		# fix ancestral state
		if (anc==alt) {
			tmp <- genos
			tmp[which(genos==0)] <- 2
			tmp[which(genos==2)] <- 0
			genos <- tmp
			rm(tmp)
			der <- ref
		} 

		# if polymorphic then it is a tested snp
		if (length(unique(as.numeric(genos)))>1) {

			# write tested snp
			cat(c(chrom, anc, der, pos, as.numeric(genos)), sep=" ", file=fout_testsnps, append=T)
			cat("\n", file=fout_testsnps, append=T)

			# if singleton
			if (length(which(genos==1))==1) {
				ind <- which(genos==1)
				singletons[ind] <- paste(singletons[ind], pos)
			}

		}

	} # end if biallelic
} # 

cat(singletons, sep="\n", file=fout_singletons)




