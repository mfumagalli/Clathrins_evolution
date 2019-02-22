
# read all ref genes
fgene <- "~/Data/UCSC/refGene.gz"
agenes <- read.table(fgene, head=F, stringsAsFact=F)

cat("All genes:", nrow(agenes))

ug <- unique(agenes[,13])
cat("\nUnique genes:", length(ug))

longest <-c()
for (i in 1:length(ug)) {
        ind <- which(agenes[,13]==ug[i])
        longest <- c(longest, ind[which.max(abs(agenes[ind,6]-agenes[ind,5]))])
}

# take only the longest isoform
genes <- agenes[longest,]
rm(agenes)

cat("\nUnique data:", nrow(genes))

save(genes, file="~/Data/UCSC/refGene.unique.Rdata")




