

args=commandArgs(T)
msdir=args[1]
rdir=args[2]
fout=args[3] # name file output for command
paramout=args[4] # name parameters file output
niter=as.numeric(args[5]) # nr of replicates
rm(args)

options(scipen=999)

# ms on Europeans only

# Gutenkunst model
# -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1

#with Nref=7310

mu=1.5e-8
L=1e5
Nref=7310
theta=4*Nref*mu*L

n_ceu=198
n_chb=206

msmstemplate = paste(msdir, " -N ", Nref," -ms ", n_ceu+n_chb, " 1 -I 3 0 ", n_ceu, " ", n_chb, " -t ", theta, " -r ", theta, " ", L, " -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1 -Sp 0.75 -SI [seltime] 3 [initialfrequencyP1] [initialfrequencyP2] [initialfrequencyP3] -Sc 0 2 [halfselstrengthP2] [selstrengthP2] [halfselstrengthP2] | gzip > msms.txt.gz", sep="",collapse="")


"%+%" = function(x,y){paste(x,y,sep="")}
complete_gsub = function(command,paramlist) {
	        for (p in names(paramlist)) command = gsub('[' %+% p %+% ']',paramlist[[p]],command,fixed=T)
        command
}

# initialise files and template
fp = file(fout,open='w');
fpar = file(paramout, open="w");
msmscommand = msmstemplate;

pams_head=paste(c("sel_time","initial_freq_AFR","initial_freq_CEU","initial_freq_CHB","sel_strength"), sep="", collapse=" ")
write(pams_head,file=fpar)

cat(c("S_5","Ss_5","TD_5","FLDs_5","FLFs_5","H1_5","H2_5","H2H1_5","FST_5",
      "S_50","Ss_50","TD_50","FLDs_50","FLFs_50","H1_50","H2_50","H2H1_50","FST_50",
      "S_100","Ss_100","TD_100","FLDs_100","FLFs_100","H1_100","H2_100","H2H1_100","FST_100"), "\n", file="summary_stats_BAL.txt")


# for each iteration
for (j in 1:niter) {

	msmscommand = msmstemplate

	paramlist = list()
	# time is pastward in units of 4N e generations
	t0=20000/(4*Nref*25)
	paramlist[['seltime']] = round(runif(1,min=t0,max=t0)*1e5)/1e5 #

	f0=0.20
	paramlist[['initialfrequencyP1']] = f0
	paramlist[['initialfrequencyP2']] = f0
	paramlist[['initialfrequencyP3']] = f0
	
	sc=runif(1,min=0,max=0.02)*(2*Nref) # 
	paramlist[['selstrengthP2']] = round(sc)
	paramlist[['halfselstrengthP2']] = round(paramlist[['selstrengthP2']]/2.0)

	# write msms commands
	msmscommand = complete_gsub(msmscommand,paramlist)
	write(msmscommand,file=fp,append=T)

	# write Rscript
	cmd=paste(rdir, "Scripts/msms2stats.R", "msms.txt.gz", L, "fake_5.txt fake_50.txt fake_100.txt", "summary_stats_BAL.txt", n_ceu, n_chb)
	write(cmd, file=fp, append=T)
	write("rm msms.txt.gz", file=fp, append=T)

	# write parameters
	pams=c(t0,f0,f0,f0,sc)
	write(paste(pams,sep="",collapse=" "),file=fpar, append=T)
	rm(pams)

}


close(fp)
close(fpar)






