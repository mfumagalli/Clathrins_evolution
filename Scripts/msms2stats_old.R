
## FROM MS FILE TO SUMMARY STATS

args=commandArgs(T)
simfile=args[1]
len=as.numeric(args[2]) # tot length
reg_100=args[3] # file with regions to consider
reg_200=args[4]
reg_500=args[5]
fout=args[6]
nsam1=as.numeric(args[7])
nsam2=as.numeric(args[8])
#print(args)
rm(args)

source("Scripts/popgen.rtx")

#cat(c("S_100","Ss_100","TD_100","FLDs_100","FLFs_100","H1_100","H2_100","H2H1_100","nSL_100","FST_100",
#"S_200","Ss_200","TD_200","FLDs_200","FLFs_200","H1_200","H2_200","H2H1_200","nSL_200","FST_200",
#"S_500","Ss_500","TD_500","FLDs_500","FLFs_500","H1_500","H2_500","H2H1_500","nSL_500","FST_500"), "\n", file=fout, append=T)

# regions
regs=c(reg_100, reg_200, reg_500)

ftemp=simfile;
outres=c()

for (x in 1:length(regs)) {
	
		haplos=readMs_sub( ftemp, nsam1+nsam2, len, freg=regs[x])
#		haplos=readMs_pos( ftemp, 230, len) # debug

		# which sites are SNPs in GR?
		#dd=countDaf(haplos$hap[195:230],haplos$hap[195], F)
		#gr_snps=which(dd[[1]]>0)
		# old way, sample the exact nr of obs SNPs
		# if (length(gr_snps)>nsub[x]) iran=sort(sample(gr_snps,nsub[x]))
		# new way, sample a fraction based on simulated-observed ration of polymorphisms, this ratio is 3.25
		#correcting_factor=3.25
#		correcting_factor=1
		#iran=sort(sample(gr_snps, round(length(gr_snps)/correcting_factor) ))
		#rm(dd); rm(gr_snps)

		pos=haplos$pos
		haplos=haplos$hap

		#if (  substring(haplos,which(pos==250000),which(pos==250000))

#		for (l in 1:length(haplos)) haplos[l]=paste(strsplit(haplos[l],split="")[[1]][iran],sep="",collapse="")
#		pos=pos[iran]

		haplist=list()
		haplist[[1]]=haplos[1:nsam1] # CEU
		haplist[[2]]=haplos[(nsam1+1): (nsam1+nsam2)] # CHB
		rm(haplos)
		# anc=paste(rep("0",nchar(haplist[[1]][1])),sep="",collapse="")

		## SUMMARY STATS on CEU
		taj=tajima(haplist[[1]])[5] #TD
		fu=fuli(haplist[[1]]) # 1 and 2 are S and SS, 3 and 4 are Ds and Fs
		hs=homohapl(haplist[[1]]) # H1, H2, H2/H1

		fsts=reynolds(haplist)

		subres=c(fu[1:2],taj,fu[3:4],hs,fsts)
		outres=c(outres, subres)

}

## print
cat(outres,"\n", file=fout, append=T)




