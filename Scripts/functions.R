
reverse<-function(seq) {
	ss=strsplit(seq,split="")[[1]]
        ss2=ss
	ss2[which(ss=="A")]="T"
	ss2[which(ss=="T")]="A"
	ss2[which(ss=="C")]="G"
	ss2[which(ss=="G")]="C"
	rm(ss)
	ss2=paste(rev(ss2),sep="",collapse="")
	ss2
}


seqdiff<-function(seq1, seq2) {
	seq1=strsplit(seq1, split="")[[1]]
	seq2=strsplit(seq2, split="")[[1]]
	dd=which(seq1!=seq2)
	res=cbind(pos=dd, aa=ceiling(dd/3), base1=seq1[dd], base2=seq2[dd])
	as.data.frame(res)
}

protdiff<-function(seq1, seq2) {
	seq1=strsplit(seq1, split="")[[1]]
        seq2=strsplit(seq2, split="")[[1]]
	dd=which(seq1!=seq2)
	res=cbind(pos=dd, aa1=seq1[dd], aa2=seq2[dd])
	as.data.frame(res)
}

translate<-function(seq) {
	nn=nchar(seq)/3
        aa=rep("?", nn)
	for (i in 1:nn) {
		ind=c(((i*3)-2),(i*3))
	        sub=substring(seq,((i*3)-2),(i*3))
		if (sub %in% c("TTT","TTC")) aa[i]="F"
		if (sub %in% c("TTA","TTG","CTT","CTC","CTA","CTG")) aa[i]="L"
		if (sub %in% c("ATT","ATC","ATA")) aa[i]="I"
		if (sub %in% c("ATG")) aa[i]="M"
		if (sub %in% c("GTT","GTC","GTA","GTG")) aa[i]="V"
		if (sub %in% c("TCT","TCC","TCA","TCG")) aa[i]="S"
		if (sub %in% c("CCT","CCC","CCA","CCG")) aa[i]="P"
		if (sub %in% c("ACT","ACC","ACA","ACG")) aa[i]="T"
		if (sub %in% c("GCT","GCC","GCA","GCG")) aa[i]="A"
		if (sub %in% c("TAT","TAC")) aa[i]="Y"
		if (sub %in% c("TAA","TAG","TGA")) aa[i]="-"
		if (sub %in% c("CAT","CAC")) aa[i]="H"
		if (sub %in% c("CAA","CAG")) aa[i]="Q"
		if (sub %in% c("AAT","AAC")) aa[i]="N"
		if (sub %in% c("AAA","AAG")) aa[i]="K"
		if (sub %in% c("GAT","GAC")) aa[i]="D"
		if (sub %in% c("GAA","GAG")) aa[i]="E"
		if (sub %in% c("TGT","TGC")) aa[i]="C"
		if (sub %in% c("TGG")) aa[i]="W"
		if (sub %in% c("CGT","CGC","CGA","CGG","AGA","AGG")) aa[i]="R"
		if (sub %in% c("AGT","AGC")) aa[i]="S"
		if (sub %in% c("GGT","GGC","GGA","GGG")) aa[i]="G"
	}
	aa=paste(aa,sep="",collapse="")
	aa
}


Dpr2<-function(haplos, fout) {

	cat("", file=fout)

	len<-nchar(haplos[1])
	nsam<-length(haplos)

	for (i in 1:(len-1)) {

		#cat("\n",i,"/",len-1)
	
		Dprime=rsquare=pv=NA;

		# alleles
		subi<-substring(haplos,i,i); subj<-substring(haplos,len,len)

		if (prod(subi==subj)) { # if exactly the same haplotypes

			Dprime=rsquare=1;pv=0;

		} else {

			if (length(unique(subi))==2) { # if biallelic SNP		

				# unique alleles
				ali<-unique(subi); alj<-unique(subj)
                        	if (length(ali)==1) ali<-c(ali, "N") # se monomorfico, soluzione temporanea
                        	if (length(alj)==1) alj<-c(alj, "N")

				# allele frequencies
				pp<-rep(NA,4)
				pp[1]<-length(which(subi==ali[1]))/nsam
				pp[2]<-length(which(subi==ali[2]))/nsam
				pp[3]<-length(which(subj==alj[1]))/nsam
				pp[4]<-length(which(subj==alj[2]))/nsam

				# haplotypes
				haps<-paste(subi,subj,sep="")

				# possible haplotypes
				posshap<-c()
				for (ai in 1:2) {
					for (aj in 1:2) {
						posshap<-c(posshap, paste(ali[ai],alj[aj],sep="",collapse=""))
					}
				}

				# haplotype frequencies
				P<-rep(NA,4)
				for (c in 1:4) P[c]<-length(which(haps==posshap[c]))/nsam

				# D, Dprime
				D<-P[1]*P[4]-P[2]*P[3]
				Dmax<-min(c(pp[1]*pp[4]),c(pp[2]*pp[3]))
				Dmin<-max(c(-pp[1]*pp[3]),c(-pp[2]*pp[4]))
				if (D<0) {
					Dprime<-abs(D/Dmin)
				} else {
					Dprime<-abs(D/Dmax)
				}

				# rsquare
				rsquare<-(D^2)/prod(pp)
				chiquadro<-rsquare*nsam
				pv<-pchisq(q=chiquadro, df=1, lower.tail=F)

			} # end if SNP biallelic

		} # fine if the same

		if (!is.na(Dprime)) {

			cat(paste(c(i,Dprime,rsquare),sep="", collapse="\t"), file=fout,append=T)
			cat("\n", file=fout,append=T)

			if ((Dprime>0.9) & (rsquare>0.1)) cat(i,Dprime,rsquare,pv,"\n")
		}

        } # fine for in i

}





