
## H1 vs H2/H1 plot for all genes

fgene <- "/home/mfumagal/Documents/Projects/ARCHIVED/Greenland_Selection/Exome_GW/refGene_unique.txt"
genes <- read.table(fgene, heade=F, stringsAsFac=F)

for (g in 1:nrow(genes)) {

	cat((g/nrow(genes))*100, "%\n")

        chrom <- as.numeric(as.character(strsplit(genes[g,3],split="chr")[[1]][2]))
        start <- genes[g,5]; end <- genes[g,5]

	if (chrom %in% seq(1,22,1)) {

	filein <- paste("~/Data/Greenland/Exome/phased/Chr",chrom,".trioChr",chrom,".filt.phased", sep="", collapse="")
	dat <- readLines(filein)

        tmp <- strsplit(dat[1],split=" ")[[1]]; ind <- tmp[3:length(tmp)]
        nsam <- length(ind)
        pos <- c(); haplos <- rep("", nsam)

	for (z in 2:length(dat)) {

     		tmp=strsplit(dat[z], split=" ")[[1]]
            	tpos=strsplit(tmp[2], split="_")[[1]][2]

                if (tpos>=start & tpos<=end & prod(nchar(tmp[3:length(tmp)])==1) & (length(which(tmp[3:length(tmp)]=="."))==0) ) {
             		pos=c(pos, tpos)
                     	for (q in 1:nsam) haplos[q]=paste(haplos[q], tmp[2+q], sep="", collapse="") 
		}
   	}

	pos <- as.numeric(pos)

	H1=H12=H2=H2H1=NA
        hfreqs=sort(as.numeric(table(unlist(haplos)))/length(haplos), dec=T)
        if (length(hfreqs)>0) H1=sum(hfreqs^2)
        if (length(hfreqs)>1) {
        	H12=H1+(2*hfreqs[1]*hfreqs[2])
                H2=H1-(hfreqs[1]^2)
                H2H1=H2/H1
     	}

	if (!is.na(H2H1)) {
		cat(genes[g,13], "\t", H1, "\t", H2H1, "\n", file="tmp.txt", append=T)
		cat((g/nrow(genes))*100, "%\n")
	}

	} # if chrom

}	


#genes <- matrix(c(9,36190852, 36212059, 5, 175819455, 175843570, 17, 57697049, 57774317, 22, 19166986, 19279239), nrow=4, ncol=3, byrow=T)
#rownames(genes)=c("CLTA","CLTB","CLTC","CLTCL1")
#colnames(genes)=c("CHROM" , "START", "END")

source("Scripts/popgen.rtx")

#fgene <- "/home/mfumagal/Documents/Projects/ARCHIVED/Greenland_Selection/Exome_GW/refGene_unique.txt"

#genes_1kg <- read.table(fgene, heade=F, stringsAsFac=F)

#lint <- 5e5

#fout <- "Plots/"

#g <- 4 # CLTCL1

#fout_here <- paste(fout,rownames(genes)[g],"_InuitSweep.pdf",sep="",collapse="")


#chr22:19,167,712-19,279,239

chrom <- 22
start <- 19167712
end <- 19279239

filein <- paste("~/Data/Greenland/Exome/phased/Chr",chrom,".trioChr",chrom,".filt.phased", sep="", collapse="")
dat <- readLines(filein)

tmp <- strsplit(dat[1],split=" ")[[1]]
ind <- tmp[3:length(tmp)]
nsam <- length(ind)

pos_tot <- pos <- c()
haplos_tot <- haplos <- rep("", nsam)

for (z in 2:length(dat)) {

	tmp=strsplit(dat[z], split=" ")[[1]]
        tpos=strsplit(tmp[2], split="_")[[1]][2]

       	if (tpos>=start & tpos<=end & prod(nchar(tmp[3:length(tmp)])==1) & (length(which(tmp[3:length(tmp)]=="."))==0) ) {
        		pos_tot=c(pos_tot, tpos)
                     	for (q in 1:nsam) haplos_tot[q]=paste(haplos_tot[q], tmp[2+q], sep="", collapse="") 
#			cat(c(chrom,"rs",0,tpos,unique(tmp[3:length(tmp)])), sep="\t", "\n", file="mytmp.bim", append=T)
		}
}
pos_tot=as.numeric(pos_tot)

# get ancestral
#system("bash Scripts/getAnc.sh")
#anc=paste(unlist(readLines("mytmp.anc")),sep="",collapse="")

#CLTCL1
#rs2793062
#chr22:19167801-19167801
which(pos_tot==19167801)
#[1] 61

haveT=which(substring(haplos_tot, 61, 61)=="T")
#2  3  4  8 10 11 17 21 27





#CLTC
#rs76405328
#chr17:57721536-57721536
#which(pos_tot==57721536)
#[1] 17
#which(substring(haplos_tot, 17, 17)=="C")
#[1] 13 17 20 23 24 27 29 33 35






















# plot color sweep for Inuit data

source("Scripts/popgen.rtx")

fin <- "Results/Genes/out_sup_CLTCL1.fasta"

haps <- readLines(fin)

chimp <- trims(haps[length(haps)])

hs <- trimws(haps[which(haps==">GI ")+1])

colorSweep(haplos=hs, outgroup=chimp, mut=9, show.labels=TRUE, show.grid=TRUE, show.ticks=TRUE, main="")

################



