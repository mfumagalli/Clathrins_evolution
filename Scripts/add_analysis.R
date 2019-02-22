
# check where the stop-codon occurs in great apes
# it is the position 64 C/T and the stop-codon is caused by T
# but it is reversed so it is at beginning of the protein

aa=readLines("Results/Apes/CLTCL1.pop.fasta")

hh=aa[seq(2,length(aa),2)]
pops=aa[seq(1,length(aa),2)]

table(substring(hh, 64, 64))

C   T 
171   5 

5/(171+5)

table(pops[which(substring(hh, 64, 64)=="T")])
>Gorilla_gorilla_gorilla 
                       5 

ind=which(substring(hh, 64, 64)=="T")

source("Scripts/popgen.R")

 tt=tajima(hh); tt
0.000000 65.000000 11.314470 13.506429  0.593029

tt=tajima(hh[-ind]); tt
0.0000000 61.0000000 10.6718873 13.2189886  0.7295555


## check frequency of 2 main haplotypes in humans (but based on ape data set)

hu=hh[which(pops==">Homo_sapiens")]
nn=c("Dai", "French", "Han", "Karitiana", "Madenka", "Mbuti", "Papuan", "San", "Sardinian")

cbind(rep(nn,each=2), match(hu, names(table(hu))[2:3]))
# the papuan can be considered as 1/2

# the mutation that separates the haplotypes is 15
# check if excess of heterozygosity

substring(hu, 15, 15)

table(substring(hu, 15, 15))

 A  G 
10  8 

f2=0.44
f1=1-f2

obs=c(2,6,1)  # AA AG GG
exp=c((f1/(f1+f2))^2*9, 2*(f1/(f1+f2))*(f2/(f1+f2))*9,(f2/(f1+f2))^2*9)

# who else has this variation?

table(substring(hh, 15, 15))
 A   G 
168   8 


## human data, how many SNPs without ancestral sequences?

aa=readLines("Results/out_sup_CLTCL1.fasta")

hh=trimws(aa[seq(2,length(aa),2)])
pops=trimws(aa[seq(1,length(aa),2)])


 tt=tajima(hh[1:5044]); tt


sort(table(hh[1:5044]))


## which one is the M1316V in the human data?
# it is the 9th
head -n 9 Results/anno_CLTCL1.csv | tail -n 1
22,19184095,rs1061325,T,C,missense

# how about raw data for archaic?










