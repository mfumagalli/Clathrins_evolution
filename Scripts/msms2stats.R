
## calculate summary stats from msms file

fin <- commandArgs(T)[1]

source("Scripts/popgen.rtx")

nsams <- c(216,198,206)
pops <- c("YRI","CEU","CHB")

haplos <- readMs(fin, sum(nsams))

# for each simulation calculate stats for each pop
for (i in 1:length(haplos)) {

	# create list for pops
	haplist <- list()
        haplist[[1]] <- haplos[[i]][1:nsams[1]] # YRI
	haplist[[2]] <- haplos[[i]][(nsams[1]+1): (nsams[1]+nsams[2])] # CEU
	haplist[[3]] <- haplos[[i]][(nsams[1]+nsams[2]+1): (nsams[1]+nsams[2]+nsams[3])] # CHB

	# for each pop
	fsts <- NA
	for (j in 1:length(nsams)) {
		taj <- tajima(haplist[[j]])[3:5] #TD
                fu <- fuli(haplist[[j]]) # 1 and 2 are S and SS, 3 and 4 are Ds and Fs
                hs <- homohapl(haplist[[j]]) # H1, H2, H12, H2/H1, H2/H1''
		fwh <- NA
		if (j==1) fsts <- reynolds(haplist)

		cat(pops[j], nsams[j], fu[1], fu[2], taj, fu[3:4], fwh, hs, fsts)
		cat("\n")
	}

}


