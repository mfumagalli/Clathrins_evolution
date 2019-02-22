
# input: haplotypes
# marker: index

args=commandArgs(T)
fin=args[1]
rm(args)

source("Scripts/functions.R")

fin="Results/Apes/CLTCL1.promoter.fasta"

haps=readLines(fin)

ind_hs=grep(">Homo_sapiens-", haps)
ind_pt=grep(">Pan_troglodytes_ellioti-", haps)

hs=trimws(haps[c(ind_hs+1, ind_hs+2)], "both")
pt=trimws(haps[c(ind_pt+1, ind_pt+2)], "both")

len_prom=500
mark=1316 # M/V in humans, M/T in pt

hs=paste(substring(hs, 1, len_prom), substring(hs, len_prom+mark, len_prom+mark), sep="")
pt=paste(substring(pt, 1, len_prom), substring(pt, len_prom+mark, len_prom+mark), sep="")

# human 1316

Dpr2(hs, "Results/Apes/CLTCL1.ld.hs.txt")
Dpr2(pt, "Results/Apes/CLTCL1.ld.pt.txt")

#Dpr2(hs, "Results/Apes/CLTCL1.ld.hs.txt")
#15 1 0.16 0.08968602 
#329 1 0.16 0.08968602 
#Dpr2(pt, "Results/Apes/CLTCL1.ld.pt.txt")
#376 1 0.2647059 0.02139757 
#416 1 0.1666667 0.06788915 
#447 1 0.6428571 0.0003361935 



