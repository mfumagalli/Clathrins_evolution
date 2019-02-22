
args <- commandArgs(T)
fin <- args[1]
flab <- args[2]
fout <- args[3]
rm(args)

# fin="Results/out_pop_CLTCL1.hetero.allelic.csv"
# flab="pops.geo.txt"
# fout="Plots/hetero_map.pdf"

library(ggmap)
library(maptools)
library(maps)

pops <- read.table(flab, sep="\t", head=T, stringsAsFact=F)

res <- read.csv(fin)

# already done
#visited <- c("SFO", "Chennai", "London", "Melbourne", "Johannesbury, SA")
#ll.visited <- geocode(visited)
#visit.x <- ll.visited$lon
#visit.y <- ll.visited$lat

freqH <- log(res[,3]/apply(X=res[,2:4],FUN=sum,MAR=1))
freqH <- freqH+1-min(freqH)


# Create some breaks and use colorRampPalette to transform the breaks into a color code
#x <- res$RatioHetero
#gr <- .bincode(x, seq(min(x), max(x), len=length(x)), include.lowest = T)
#col <- colorRampPalette(c("red", "white", "yellow"))(length(x))[gr]

#tmp <- 1-res$RatioHetero
#cols <- grey(tmp-min(tmp), alpha = NULL)
     
# pdf(file="Plots/hetero_map.pdf")

pchs <- rep(15, nrow(pops))
qq <- quantile(res$RatioHetero, seq(0,1,1/3))
pchs[which(res$RatioHetero>=qq[3])] <- 17
pchs[which(res$RatioHetero<qq[2])] <- 16

cols <- rep("gray30", nrow(pops))
cols[which(res$RatioHetero>=qq[3])] <- "red"
cols[which(res$RatioHetero<qq[2])] <- "blue"


#> table(pchs)
#pchs
#15 16 17 
#10  8  9

if (prod(pops[,1]==res[,1])) {

	pdf(file=fout)

	map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
	points(pops$lon, pops$lat, col=cols, pch=pchs, cex=freqH)

	dev.off()

}





