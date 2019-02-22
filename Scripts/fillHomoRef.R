
args <- commandArgs(T)
fin <- args[1]
nheader <- as.numeric(args[2])
rm(args)

vcf <- read.table(fin, stringsAsFact=F)

fillHomo <- function(arr) {
	arr[which(arr==".")]="0/0"
	arr[3]="." # ID
	arr
}

vcf2 <- t(apply(X=vcf, MAR=1, FUN=fillHomo))

cat(readLines(fin, n=(nheader-1)), sep="\n")
cat(strsplit(readLines(fin, n=nheader)[nheader], split="\t")[[1]], sep="\t")
cat("\n")
for (i in 1:nrow(vcf2)) {
	cat(vcf2[i,], sep="\t")
	cat("\n")
}



