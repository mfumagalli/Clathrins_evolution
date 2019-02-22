
fin=commandArgs(T)

# read haplo
ss=readLines(fin)
pops=substring(trimws(ss[seq(1,length(ss),2)]),2)
haps=trimws(ss[seq(2,length(ss),2)])

upop=setdiff(sort(unique(pops)), c("CHIMP","ALTAI","DENISOVA"))

alles=names(rev(sort(table(haps)))[1:2])

cat("Pop.,Obs.Homo1,Obs.Hetero,Obs.Homo2,Exp.Homo1,Exp.Hetero,Exp.Homo2,RatioHetero,p-value\n")

for (i in 1:length(upop)) {

	#cat(upop[i],"\n")
	ind=ceiling(which(pops==upop[i])/2)
	nind=length(ind)/2
	nind_eff=0; f1=0; f2=0; # samples with both alleles (MAJOR HAPLOTYPES)
	homo1=hetero=homo2=0

	for (j in ind) {

		j2=j*2; j1=j2-1
		if (prod(haps[j1:j2] %in% alles)) {
			nind_eff=nind_eff+1;
			tmp1=sum(haps[j1:j2] %in% alles[1])
			tmp2=sum(haps[j1:j2] %in% alles[2])
			f1=f1+tmp1;
			f2=f2+tmp2;
			if (tmp1==2) homo1=homo1+1
			if (tmp2==2) homo2=homo2+1
			if (tmp1==1) hetero=hetero+1	
		}
	}

	if(!((homo1+homo2+hetero) == nind_eff)) cat("Error")
	if (!((f1+f2)/2 == nind_eff)) cat("Error")

	obs=c(homo1,hetero,homo2)
	exp=c((f1/(f1+f2))^2*nind_eff, 2*(f1/(f1+f2))*(f2/(f1+f2))*nind_eff,(f2/(f1+f2))^2*nind_eff)

	pv <- chisq.test(matrix(c(obs,exp), nrow=2, ncol=3, byrow=T))$p.value

	#cat("OBS:", obs, "\n")
	#cat("EXP:", exp, "\n")
	#cat("ratio:", obs[2]/exp[2], "\n\n")

	cat(upop[i], obs, exp, obs[2]/exp[2], pv, sep=",", "\n")

}



