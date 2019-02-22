
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
ins[1] <- paste("~/Data/CLTCL1/Genes/", gene, ".vcf", sep="", collapse="")
ins[2] <- "~/Data/CLTCL1/names.txt"
ins[3] <- paste("~/Data/Greenland/Exome/phased/Chr",chrom,".trioChr",chrom,".filt.phased",sep="",collapse="")
ins[4] <- paste("~/Data/CLTCL1/Genes/", gene, ".Altai.vcf", sep="", collapse="")
ins[5] <- paste("~/Data/CLTCL1/Genes/", gene, ".Denisova.vcf", sep="", collapse="")

# output files (in Results folder)
outs <- c()
outs[1] <- paste("Results/Genes/out_sup_", gene, ".fasta", sep="", collapse="")
outs[2] <- paste("Results/Genes/out_id_", gene, ".fasta", sep="", collapse="")
outs[3] <- paste("Results/Genes/out_pop_", gene, ".fasta", sep="", collapse="")

fanno <- paste("Results/Genes/anno_", gene,".csv", sep="", collapse="")
fcod <- paste("Results/Genes/cods_", gene, ".txt", sep="", collapse="")

# other files needed
fixeds <- c()
fixeds[1] <- "~/Data/1000G/related_samples_panel.20140910.ALL.panel"
fixeds[2] <- "~/Data/1000G/integrated_call_samples_v3.20130502.ALL.panel"

cat("\nUsing these fixed files:", "\n1000G ped file:", fixeds[1], "\n1000G panel file:", fixeds[2])

cat("\nAssuming these are input files:", "\nVCF 1000G:", ins[1], "\nNames ID 1000G:", ins[2], "\nsnp file:", paste("~/Data/UCSC/SNP/snps_", chrom, ".Rdata", sep="", collapse=""), "\nPhased data Inuit:", ins[3], "\nVCF Altai:", ins[4], "\nVCF Denisova:", ins[5])

cat("\nLooking at this region:", chrom,":",chrom_start,"-",chrom_end)

# cat("\nProducing these output files:", outs, fanno, fcod)
cat("\nProducing these output files:", outs, fanno, fcod)

## read 1000G vcf file
vcf <- read.table(ins[1], sep="\t", stringsAsFact=F)
nind <- ncol(vcf)-9
nsites <- nrow(vcf)
nam <- scan(ins[2], what="char", quiet=T) #IDs samples
cat("\nHow many individuals?", nind, "and are they correct?", (length(nam)==nind))

## preliminary: check unrelated individuals ID
doUnrel <- 1
if (doUnrel) {
        ids <- as.character(read.table(file=fixeds[1], stringsAsFact=F, sep="\t", head=F)[,1])
        related <- length(which(nam %in% ids))
        if (related) stop("There are related samples. Remove them and restart.")
	rm(ids)
}

vcf <- vcf[,c(1:5,10:ncol(vcf))]
cat("\nHow many sites?", nsites)

## read pop and superpop information from panel file
pops <- read.table(fixeds[2], stringsAsFactors=F, head=T)

cat("\nReading SNPs...")
load(paste("~/Data/UCSC/SNP/snps_chr", chrom, ".Rdata", sep="", collapse=""))
## select only rs for coding snps, and annotated function, and select chimp allele

## read Inuit
inu <- read.table(file=ins[3], stringsAsFactors=F, head=T)
num_inu <- 18 # I assume we have 18 chromosomes (9 individuals)
inu_pos <- unique(sort(as.numeric(unlist(strsplit(inu[,2], split=paste(chrom,"_",sep="",collapse=""))))))
# get only snps in the gene
# e.g. 22:19167732-19279164
inu_ind_pos <- which(inu_pos>=chrom_start & inu_pos<=chrom_end)
# inu_pos[inu_ind_pos] %in% vcf[,2] # are all snps there?
# so it seems that are not extra snps in the Inuit
# there is one missing SNP rs2793062 with ambiguous position but it is not coding anyway

cat("\nReading VCFs...\n")
## read Altai
altai_vcf <- read.table(ins[4], sep="\t", stringsAsFact=F)
## read Denisova
deni_vcf <- read.table(ins[5], sep="\t", stringsAsFact=F)

## creata fasta for all individuals, plus chimp and archaic
# initialise
fas0 <- fas1 <- c()
for (i in 1:(length(nam)+num_inu)) fas0[i] <- fas1[i] <- ""
chimp <- altai <- deni <- ""

## cycle across all sites and fill in merging inuit and chimp and archaic
# keep annotation of coding sites  
anno <- cods <- c()

# for each site in the 1000G dataset
for (i in 1:nsites) {

	# alternate and reference alleles
	ref <- vcf[i,4]; alt <- vcf[i,5];

	# matching by rs and also the alternate must be unique
	if ( (vcf$V3[i] %in% snps$V5) & (nchar(alt)==1) ) {

		# get rs and function
		rs <- vcf$V3[i]
		pos <- vcf$V2[i]
		func <- snps$V16[which(snps$V5==rs)]

		# if non-synonymous
		 if ( length(grep("missense",func)) | length(grep("nonsense",func)) | length(grep("stop-loss",func)) | length(grep("frameshift",func)) ) {
		
			cat(".")

			# extract annotation
			anno <- rbind(anno, cbind(vcf[i,1:5],func))
			cods <- rbind(cods, snps_cod[which(snps_cod$V5==rs),])

			# chimp allele
			ch_all <- snps_ortho[which(snps_ortho$V5==rs),12]
			if (length(ch_all)==0) ch_all <- ref
			chimp <- paste(chimp, ch_all, sep="", collapse="")

			# altai allele
			ind_altai <- which(altai_vcf$V2==pos)
			if (length(ind_altai)==1) {
				geno_altai <- sum(as.numeric(strsplit(strsplit(altai_vcf[ind_altai,10],split=":")[[1]][1],split="/")[[1]]))
				if (geno_altai<2) altai <- paste(altai, ref, sep="",collapse="") else altai <- paste(altai, alt, sep="",collapse="")
			} else {
				altai <- paste(altai, chimp, sep="",collapse="")
			}

			# denisova allele
			ind_deni <- which(deni_vcf$V2==pos)
			if (length(ind_deni)==1) {
				geno_deni <- sum(as.numeric(strsplit(strsplit(deni_vcf[ind_deni,10],split=":")[[1]][1],split="/")[[1]]))
				if (geno_deni<2) deni <- paste(deni, ref, sep="",collapse="") else deni <- paste(deni, alt, sep="",collapse="")
			} else {
				deni <- paste(deni, chimp, sep="",collapse="")
			}

			# 1000G 
			for (j in 1:length(nam)) {
				s <- j+5
				all <- strsplit(vcf[i,s],split="|")[[1]]
				if (as.numeric(all[1])==0) fas0[j] <- paste(fas0[j], ref, sep="", collapse="") else fas0[j] <- paste(fas0[j], alt, sep="", collapse="")
				if (as.numeric(all[3])==0) fas1[j] <- paste(fas1[j], ref, sep="", collapse="") else fas1[j] <- paste(fas1[j], alt, sep="", collapse="")
			}

			# Inuit
			ii <- which(inu[,2]==paste(chrom,"_",vcf[i,2],sep="",collapse=""))
			# if (length(ii)>0) cat("\n",i,":",ii) # these are multiallelic snps in Inuit
			#print(sort(unique(as.character(inu[ii,3:38])))==sort(c(ref,alt)))
			tt <- 3
			for (j in (length(nam)+1):length(fas0)) {
				if (length(ii)==0) {
					fas0[j] <- paste(fas0[j], ref, sep="", collapse="")
					fas1[j] <- paste(fas1[j], ref, sep="", collapse="")
				} else {
					fas0[j] <- paste(fas0[j], as.character(inu[ii,tt]), sep="", collapse="")
					fas1[j] <- paste(fas1[j], as.character(inu[ii,(tt+1)]), sep="", collapse="")
				}
				tt <- tt+2	
			}
		}
	} else { # end if not matched by rs
		if ( (vcf$V2[i] %in% inu_pos[inu_ind_pos]) ) cat("\nWarning! Missed an Inuit SNP:", i, ";", unlist(vcf[i,1:5]),".To be checked.")
	}
}

rm(snps)
rm(snps_cod)
rm(snps_ortho)

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

# inuit
for (i in (length(nam)+1):length(fas0)) {

	# id
	cat(paste(">GI_0_",i,sep="",collapse=""),"\n", file=outs[2], append=T)
	cat(fas0[i],"\n", file=outs[2], append=T)
	cat(paste(">GI_1_",i,sep="",collapse=""),"\n", file=outs[2], append=T)
	cat(fas1[i],"\n", file=outs[2], append=T)

	# superpop=pop
	cat(paste(">GI",sep="",collapse=""),"\n", file=outs[1], append=T)
        cat(fas0[i],"\n", file=outs[1], append=T)
        cat(paste(">GI",sep="",collapse=""),"\n", file=outs[1], append=T)
	cat(fas1[i],"\n", file=outs[1], append=T)

	# pop=superpop
	cat(paste(">GI",sep="",collapse=""),"\n", file=outs[3], append=T)
	cat(fas0[i],"\n", file=outs[3], append=T)
	cat(paste(">GI",sep="",collapse=""),"\n", file=outs[3], append=T)
	cat(fas1[i],"\n", file=outs[3], append=T)

}


for (i in 1:3) {

	cat(">ALTAI\n", file=outs[i], append=T)
	cat(altai, "\n", file=outs[i], append=T)

      	cat(">DENISOVA\n", file=outs[i], append=T)
    	cat(deni, "\n", file=outs[i], append=T)

	cat(">CHIMP\n", file=outs[i], append=T)
	cat(chimp, "\n", file=outs[i], append=T)

}


# coordinates
write.table(anno, file=fanno, sep=",", quote=F, row.names=F, col.names=F)

# coordinates
write.table(cods, file=fcod, sep="\t", quote=F, row.names=F, col.names=T)

# samples
# write.table(pops[match(nam, pops$sample),], sep="\t", quote=F, row.names=F, col.names=T, file=fnam)

# cat("\nWritten these files:", outs, fanno, fcod, fnam, sep="\t", "\n")
cat("\nWritten these files:", outs, fanno, fcodsep="\t", "\n")




