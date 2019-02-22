
# read snp tables from UCSC and save R objects
Rscript Scripts/readSNPtable.R

# create haplotypes for 1000G data
bash run.sh CLTCL1 22 19167732 19279164
bash run.sh CLTC 17 57697050 57774317

# get unique set of genes
Rscript Scripts/getUniqueGenes.R

# get null distribution of genes
Rscript Scripts/getNull.R 0.05 Results/null.txt

NG=`wc -l Results/null.txt | cut -f 1 -d " "`

# calculate haplotypes and summary stats for control genes
for I in `seq 1 $NG`
do
	GENE=`head -n $I Results/null.txt | tail -n 1 | cut -f 1`
	CHROM=`head -n $I Results/null.txt | tail -n 1 | cut -f 2`
	START=`head -n $I Results/null.txt | tail -n 1 | cut -f 3`
	END=`head -n $I Results/null.txt | tail -n 1 | cut -f 4`

	echo $I $GENE $CHROM $START $END

	bash run.sh $GENE $CHROM $START $END >> log.txt
done

## Chimp data
# create haplotypes for chimp data
bash run_chimp.sh CLTCL1 22 19167732 19279164
bash run_chimp.sh CLTC 17 57697050 57774317

## Apes
# apes are on hg18: CLTCL1 22 17546987 17659239
bash run_ape.sh CLTCL1 22 17536987 17669239 # aa 1641
bash run_ape.sh CLTC 17 55041832 55139099 # aa



