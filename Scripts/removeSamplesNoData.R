
args <- commandArgs(T)
fin_ss <- args[1]
fin_te <- args[2]
fout_ss <- args[3]
fout_te <- args[4]
rm(args)

# chrom22.mod.singletons chrom22.testsnp chrom22.filt.singletons chrom22.filt.testsnp

ss <- readLines(fin_ss)
nodata <- which(nchar(ss)==0)
cat(ss[-nodata], sep="\n", file=fout_ss)

maxSS <- 0
sss <- strsplit(ss, split=" ")
for (i in 1:length(ss)) maxSS <- max(c(maxSS, length(sss[[i]])  ))
cat("Highest number of singletons:", maxSS, ".\n")

nodata <- nodata+4
cat("Number of samples with data:", length(ss)-length(nodata), ".\n")
te <- read.table(fin_te, sep=" ", stringsAsFact=F, head=F)
write.table(te[,-nodata], sep=" ", quote=F, row.names=F, col.names=F, file=fout_te)

cat("Range of positions:", range(te$V4), ".\n")




