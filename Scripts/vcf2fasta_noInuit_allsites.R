
# merge Altai, Denisova, Inuit to 1000G
# select only missense
# create fasta

# prepare file names
args <- commandArgs(T)
gene <- toupper(args[1])
chrom <- as.numeric(args[2])
chrom_start <- as.numeric(args[3])
chrom_end <- as.numeric(args[4])
rm(args)

# input files
ins <- c()
ins[1] <- paste("~/Data/CLTCL1/Genes_allsites/", gene, ".vcf", sep="", collapse="")
ins[2] <- "~/Data/CLTCL1/names.txt"
ins[3] <- ""
#ins[4] <- paste("~/Data/CLTCL1/Genes/", gene, ".Altai.vcf", sep="", collapse="")
#ins[5] <- paste("~/Data/CLTCL1/Genes/", gene, ".Denisova.vcf", sep="", collapse="")

# output files (in Results folder)
outs <- c()
outs[1] <- paste("Results/Genes_allsites/out_sup_", gene, ".fasta", sep="", collapse="")
outs[2] <- paste("Results/Genes_allsites/out_id_", gene, ".fasta", sep="", collapse="")
outs[3] <- paste("Results/Genes_allsites/out_pop_", gene, ".fasta", sep="", collapse="")

# other files needed
fixeds <- c()
fixeds[1] <- "~/Data/1000G/related_samples_panel.20140910.ALL.panel"
fixeds[2] <- "~/Data/1000G/integrated_call_samples_v3.20130502.ALL.panel"

#cat("\nUsing these fixed files:", "\n1000G ped file:", fixeds[1], "\n1000G panel file:", fixeds[2])

#cat("\nAssuming these are input files:", "\nVCF 1000G:", ins[1], "\nNames ID 1000G:", ins[2], "\nsnp file:", paste("~/Data/UCSC/SNP/snps_", chrom, ".Rdata", sep="", collapse=""), "\nVCF Altai:", ins[4], "\nVCF Denisova:", ins[5])

#cat("\nLooking at this region:", chrom,":",chrom_start,"-",chrom_end)

#cat("\nProducing these output files:", outs)

## read 1000G vcf file
vcf <- read.table(ins[1], sep="\t", stringsAsFact=F)
nind <- ncol(vcf)-9
nsites <- nrow(vcf)
nam <- scan(ins[2], what="char", quiet=T) #IDs samples
#cat("\nHow many individuals?", nind, "and are they correct?", (length(nam)==nind))

## preliminary: check unrelated individuals ID
doUnrel <- 1
if (doUnrel) {
        ids <- as.character(read.table(file=fixeds[1], stringsAsFact=F, sep="\t", head=F)[,1])
        related <- length(which(nam %in% ids))
        if (related) stop("There are related samples. Remove them and restart.")
	rm(ids)
}

vcf <- vcf[,c(1:5,10:ncol(vcf))]
cat("How many sites?", nsites)

## read pop and superpop information from panel file
pops <- read.table(fixeds[2], stringsAsFactors=F, head=T)

#cat("\nReading VCFs...")
## read Altai
#altai_vcf <- read.table(ins[4], sep="\t", stringsAsFact=F)
## read Denisova
#deni_vcf <- read.table(ins[5], sep="\t", stringsAsFact=F)

## creata fasta for all individuals, plus archaic (no chimp)
# initialise
fas0 <- fas1 <- c()
for (i in 1:length(nam)) fas0[i] <- fas1[i] <- ""
#altai <- deni <- ""

# cycle across all sites and fill in merging archaic

# for each site in the 1000G dataset
for (i in 1:nsites) {

	# alternate and reference alleles
	ref <- vcf[i,4]; alt <- vcf[i,5];

	# alleles must be unique
	if ( (nchar(ref)==1) & (nchar(alt)==1)) {

		# get rs and function
		rs <- vcf$V3[i]; pos <- vcf$V2[i]
		cat(".")

		# 1000G 
		for (j in 1:length(nam)) {
			s <- j+5
			all <- strsplit(vcf[i,s],split="|")[[1]]
			if (as.numeric(all[1])==0) fas0[j] <- paste(fas0[j], ref, sep="", collapse="") else fas0[j] <- paste(fas0[j], alt, sep="", collapse="")
			if (as.numeric(all[3])==0) fas1[j] <- paste(fas1[j], ref, sep="", collapse="") else fas1[j] <- paste(fas1[j], alt, sep="", collapse="")
		}

	}
}


## print to files
for (i in 1:3) cat("", file=outs[i])

# 1000g
for (i in 1:length(nam)) {

	pp <- pops$pop[pops$sample==nam[i]]
	sp <- pops$super_pop[pops$sample==nam[i]]

	# id
	cat(paste(">",nam[i], "_",pp,"_",sp, ".0",sep="",collapse=""),"\n", file=outs[2], append=T)
	cat(fas0[i],"\n", file=outs[2], append=T)
	cat(paste(">",nam[i], "_",pp,"_",sp, ".1",sep="",collapse=""),"\n", file=outs[2], append=T)
        cat(fas1[i],"\n", file=outs[2], append=T)

	# superpop
	cat(paste(">",sp,sep="",collapse=""),"\n", file=outs[1], append=T)
	cat(fas0[i],"\n", file=outs[1], append=T)
	cat(paste(">",sp,sep="",collapse=""),"\n", file=outs[1], append=T)
	cat(fas1[i],"\n", file=outs[1], append=T)

	# pop
	cat(paste(">",pp,sep="",collapse=""),"\n", file=outs[3], append=T)
	cat(fas0[i],"\n", file=outs[3], append=T)
	cat(paste(">",pp,sep="",collapse=""),"\n", file=outs[3], append=T)
	cat(fas1[i],"\n", file=outs[3], append=T)

}


cat("\nWritten these files:", outs, sep="\t", "\n")




