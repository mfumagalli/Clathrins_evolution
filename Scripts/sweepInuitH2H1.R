
## H1 vs H2/H1 plot for all genes

fgene <- "/home/mfumagal/Documents/Projects/ARCHIVED/Greenland_Selection/Exome_GW/refGene_unique.txt"
genes <- read.table(fgene, heade=F, stringsAsFac=F)

cat("", file="tmp.txt")

for (g in 1:nrow(genes)) {

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

        if(!is.na(H2H1)) {
                cat(genes[g,13], "\t", H1, "\t", H2H1, "\n", file="tmp.txt", append=T)
                cat((g/nrow(genes))*100, "%\n")
        }

        } # if chrom

}



