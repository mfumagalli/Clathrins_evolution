
## extract control genes with similar coding length

load("~/Data/UCSC/refGene.unique.Rdata")
# "genes"

args <- commandArgs(T)
th <- as.numeric(args[1])
fout <- args[2]
rm(args)

# calc exonic lengths for each gene

# which(genes[,13]=="CLTCL1")

cat("initial set:" , nrow(genes))

gg <- c()
for (i in 1:nrow(genes)) {

	if (!( genes[i,3] %in% c("chr6_dbb_hap3","chr6_ssto_hap7","chr6_cox_hap2","chr17_ctg5_hap1","chr6_mcf_hap5","chr6_mann_hap4","chr6_qbl_hap6","chrUn_gl000228","chr6_apd_hap1","chr19_gl000209_random","chrUn_gl000213","chr7_gl000195_random","chrX","chrY"))) {
		starts <- as.numeric(strsplit(genes[i,10], split=",")[[1]])
		ends <- as.numeric(strsplit(genes[i,11], split=",")[[1]])
		len <- sum(ends-starts)
		gg <- rbind(gg, cbind(i,genes[i,c(3,5,6,13)],len))
	}
}

ind <- which(gg[,5]=="CLTCL1")
t_len <- gg[ind,6]
cat("target len:", t_len)

ind_c <- which( (gg[,6]<=(t_len+(t_len*th))) & (gg[,6]>=(t_len-(t_len*th))) & !(gg[,2] %in% c("chr6")) & !(gg[,5] %in% c("CLTC","CLTCL1","CLTA","CLTB")) )

controls <- gg[ind_c,]
controls <- controls[sort(controls[,6], ind=T)$ix,]

cat("nr controls:", nrow(controls))

# now further filter so that CLTLC1 is the exact median

more <- length(which(controls[,6]>t_len))
less <-  length(which(controls[,6]<t_len))

if (more<less) { 
	controls <- controls[-(1:(less-more)),]
} else {
	stop("to implement!")
}

cat("nr controls:", nrow(controls))

cat("", file=fout)
for (i in 1:nrow(controls)) {
	chr <- as.numeric(strsplit(controls[i,2],split="chr")[[1]][2])
	cat(c(controls[i,5],chr,controls[i,3],controls[i,4]), sep="\t", "\n", file=fout, append=T)
}




