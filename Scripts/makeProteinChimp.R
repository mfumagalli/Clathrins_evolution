
source("Scripts/functions.R")

genename <- toupper(commandArgs(T))

# check whether chimp snps in vcf are nonsyno or not
# also extract 5' promoter region
# extract reference gene sequence hg19 vcf file

# retrieve table with coding exons (from cds start to end, and strand)
refgene <- read.table("Data/refGene.txt")
colnames(refgene) <- strsplit(readLines("Data/refGene.txt")[1], split="\t")[[1]]

ii <- which(refgene$name2==genename)
if (length(ii)==1) refgene <- refgene[ii,] else stop("Gene name not found or multiple hits.");
rm(ii)

strand <- as.character(refgene$strand)
chrom <- as.numeric(strsplit(as.character(refgene$chrom), split="chr")[[1]][2])
starts <- as.numeric(strsplit(as.character(refgene$exonStarts), split=",")[[1]])+1 # because they are 0-based
ends <- as.numeric(strsplit(as.character(refgene$exonEnds), split=",")[[1]])

cat("Strand:", strand, "\n")

# retrieve nucleotide of coding (bases) + promotorial 5' (proms)
bases <- c()
for (i in 1:length(starts)) {
	primo=ultimo=0;
	if (refgene$cdsStart<=ends[i] & refgene$cdsEnd>=starts[i]) {
		if (refgene$cdsStart>=starts[i] & refgene$cdsStart<=ends[i]) primo=1 else primo=0;
	        if (refgene$cdsEnd>=starts[i] & refgene$cdsEnd<=ends[i]) ultimo=1 else ultimo=0;
		if (primo==1) starts[i]=refgene$cdsStart+1 # 0-based
		if (ultimo==1) ends[i]=refgene$cdsEnd
		bases=c(bases, seq(starts[i],ends[i]))
		#cat(i,":",length(seq(starts[i],ends[i])),"=",length(bases),"\n")
	} else cat(i,"non-coding exon skipped\n")
}
cat(length(bases)," ", length(bases)/3)
# promotorial 5k upstream txStart/end
if (strand=="+") proms=seq(refgene$txStart-10000, refgene$txStart-1) else proms=seq(refgene$txEnd+1, refgene$txEnd+10000)

# get reference sequence for coding
mine=c()
for (i in 1:length(bases)) {
	cmd=paste(paste("~/Software/samtools-1.9/samtools faidx ~/Data/hg19/chr", chrom, ".fa.gz", sep="", collapse=""), paste(as.character(refgene$chrom),":",paste(rep(bases[i],2),sep="",collapse="-"),sep="",collapse=""))
	mine=c(mine, system(cmd, intern=T)[2])
}
mine=paste(mine,sep="",collapse="")
# and for promotorial
promot=c()
for (i in 1:length(proms)) {
	        cmd=paste(paste("~/Software/samtools-1.9/samtools faidx ~/Data/hg19/chr", chrom, ".fa.gz", sep="", collapse=""), paste(as.character(refgene$chrom),":",paste(rep(proms[i],2),sep="",collapse="-"),sep="",collapse=""))
        promot=c(promot, system(cmd, intern=T)[2])
}
promot=toupper(paste(promot,sep="",collapse=""))

# put reverse in case
if (strand=="-") {
	mine=reverse(mine)
	promot=reverse(promot)
	bases=rev(bases)
	proms=rev(proms)
}

# translate
prot_ref=translate(mine)

# now read VCF file and for each snp check whether the translated proteins with ref/alt are different
# modify these options according to the VCF specific features:
# n_line_header
n_line_header=6
fin=paste("~/Data/CLTCL1/Chimp/", genename, ".phased.vcf", sep="", collapse="")
vcf=read.table(fin, sep="\t", stringsAsFact=F)
nind=ncol(vcf)-9
heads=strsplit(readLines(fin,n=50)[6], split="\t")[[1]]
colnames(vcf)=heads
names=heads[10:length(heads)]

## CODING

# cycle across all sites that are coding
# bases are the coding sites
ind=match(vcf$POS, bases) # coding
ind_prom=match(vcf$POS, proms)
# vcf=vcf[which(!is.na(ind)),]
nsites=nrow(vcf)
#ind_cod=ind_cod[which(!is.na(ind_cod))]
#ind_prom=ind_prom[which(!is.na(ind_prom))]

# this is to get annotation sites that are nonsynonymous/promotirla and haplotypes for nucleotides

# if they are nonsyno or promotor, retain the nucleotide in a "haplotype"-like format
fhaps=c(); for (i in 1:(nind*2)) fhaps[i]=""

# if they are nonsyno, retain the amino acid in a "haplotype"-like format
phaps=matrix("?",nrow=nind*2,ncol=nchar(prot_ref))
for (i in 1:nrow(phaps)) phaps[i,]=strsplit(prot_ref, split="")[[1]]

anno=genomic=c()
for (i in 1:nsites) {

	cat("\ni")

	# if CODING
	if (!is.na(ind[i])) {

		if (strand=="-") {
			if (!(as.character(vcf[i,4])==reverse(substring(mine,ind[i],ind[i])))) {
				cat("\nreferences do not match!")
			  	cat(i, as.character(vcf[i,4]), reverse(substring(mine,ind[i],ind[i])))	
			} #stop("references do not match!");
		} else {
			if (!(as.character(vcf[i,4])==(substring(mine,ind[i],ind[i])))) {
				#stop("references do not match!");
				cat("\nreferences do not match!")
				cat(i, as.character(vcf[i,4]), (substring(mine,ind[i],ind[i])))
			}
		}

		# put them in reverse in case
		ref=as.character(vcf[i,4])
		alt=as.character(vcf[i,5])
		if (strand=="-") {
			ref=reverse(ref)
			alt=reverse(alt)
		}

		tmp=strsplit(mine,split="")[[1]]; tmp[ind[i]]=alt;
		tmp=paste(tmp,sep="",collapse="")
		prot_alt=translate(tmp)

		dd=protdiff(prot_ref, prot_alt)

		# if nonsyno
		if (nrow(dd)==1) {

			cat("nonsyno")

			# extract nucleotides, aminoacids and positions
			sd=seqdiff(mine,tmp)
			anno=rbind(anno, cbind(vcf[i,1:2],sd, dd[-1]))
			pref=as.character(unlist(dd[2]))
			palt=as.character(unlist(dd[3]))
			aa_pos=as.numeric(as.character(unlist(dd[1])))
			cat("\n",i,"(",aa_pos,"), ref:", ref, "(",pref,"), alt:", alt, "(", palt, ")")

			# for each individual fill in the protein if alternate
			for (j in 1:nind) {

				jj=j+9;j2=(j*2);j1=j2-1; js=c(j1,j2)

				alls=as.numeric(strsplit(vcf[i,jj],split="\\|")[[1]])

				if (is.na(sum(alls))) stop("alleles not defined in VCF!")

				for (z in 1:length(js)) if (alls[z]==1) phaps[js[z], aa_pos]=palt;

			}
		}
	} # end if coding

	if (!is.na(ind_prom[i])) {

		if (strand=="-") {
	       		if (!(as.character(vcf[i,4])==reverse(substring(promot,ind_prom[i],ind_prom[i])))) stop("references do not match!");
                } else {
		       	if (!(as.character(vcf[i,4])==(substring(promot,ind_prom[i],ind_prom[i])))) stop("references do not match!");
	        }

                # put them in reverse in case
                ref=as.character(vcf[i,4])
	        alt=as.character(vcf[i,5])
	        if (strand=="-") {
			ref=reverse(ref)
			alt=reverse(alt)
		}

		# for each individual fill in the nucleotide
		for (j in 1:nind) {
			jj=j+9;j2=(j*2);j1=j2-1; js=c(j1,j2)

		        alls=as.numeric(strsplit(vcf[i,jj],split="\\|")[[1]])
			if (is.na(sum(alls))) stop("alleles not defined in VCF!")

                        for (z in 1:length(js)) {
				if (alls[z]==1) fhaps[js[z]]=paste(fhaps[js[z]], alt, sep="", collapse="") else fhaps[js[z]]=paste(fhaps[js[z]], ref, sep="", collapse="")
			}

		}
		genomic=c(genomic, proms[ind_prom[i]])

	} # end if prom

}
cat("\nLength promoter:", length(genomic))
write.table(anno, file=paste("Results/Chimp/", genename,".anno.txt", sep="", collapse=""), col.names=T, row.names=F, quote=F, sep="\t")
cat(genomic, file=paste("Results/Chimp/", genename,".promoter.txt", sep="", collapse=""), sep="\n")

# write fasta but for the whole protein

fout=paste("Results/Chimp/", genename,".fasta",sep="",collapse="")
cat("", file=fout)
for (i in 1:nind) {
	cat(paste(">",names[i],sep="",collapse=""), "\n", file=fout, append=T)
	cat(paste(phaps[((i*2)-1),], sep="", collapse=""), "\n", file=fout, append=T)
	cat(paste(phaps[(i*2),], sep="", collapse=""), "\n", file=fout, append=T)
}

# write fasta for promoter + the whole protein

fout=paste("Results/Chimp/", genename,".promoter.fasta",sep="",collapse="")
cat("", file=fout)
for (i in 1:nind) {
        cat(paste(">",names[i],sep="",collapse=""), "\n", file=fout, append=T)
        cat(paste(fhaps[((i*2)-1)], paste(phaps[((i*2)-1),], sep="", collapse=""), sep="", collapse=""), "\n", file=fout, append=T)
	cat(paste(fhaps[(i*2)],     paste(phaps[(i*2),], sep="", collapse=""), sep="", collapse=""), "\n", file=fout, append=T)
}
rm(phaps)
rm(fhaps)




