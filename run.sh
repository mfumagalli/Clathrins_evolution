
# INFO:
# This scripts takes a gene name and its coordinates and retrieve haplotypic data for humans and archaic humans for nonsyno mutations only.
# It uses VCF files from 1000 Genomes

GENE=$1
CHROM=$2
START=$3
END=$4

# for instance
# CLTCL1 22:19167732-19279164
# CLTC 17:57697050-57774317

echo $GENE $CHROM:$START-$END

# paths to where the chromosome-specific VCF files are
VCF=~/Data/1000G/VCF/ALL.chr$CHROM.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# path to vcflib
LIB=~/Software/vcflib/bin
# path to where the info on all control genes will be stored
DIR=~/Data/CLTCL1/Genes

## filter human
echo filtering?
if [ ! -f $DIR/$GENE.vcf ];
then	
	echo human
	if [ ! $GENE = "FMO6P" ]; # FMO6P has a weird file
	then
		$LIB/vcffilter -f "VT = SNP & QUAL > 1 & AC > 1 & DP > 5000" -r $CHROM:$START-$END $VCF | grep EX_TARGET | tail -n +2 > tmp
	else
		$LIB/vcffilter -f "VT = SNP & QUAL > 1 & AC > 1 & DP > 5000" -r $CHROM:$START-$END $VCF | tail -n +2 > tmp
	fi
	zcat $VCF | head -n 253 > head.txt
	cat head.txt tmp > $DIR/$GENE.vcf
	rm tmp
	rm head.txt
fi

## filter Altai
if [ ! -f $DIR/$GENE.Altai.vcf ];
then
	echo altai
	ALTAI=~/Data/Altai/AltaiNea.hg19_1000g.$CHROM.mod.indexed.vcf.gz
	$LIB/vcffilter -f "QUAL > 1 & DP > 10" -r $CHROM:$START-$END $ALTAI > $DIR/$GENE.Altai.vcf
fi

## filter Denisova
if [ ! -f $DIR/$GENE.Denisova.vcf ];
then
	echo denisova
	DENI=~/Data/Denisova/DenisovaPinky.hg19_1000g.$CHROM.mod.indexed.vcf.gz
	$LIB/vcffilter -f "QUAL > 1 & DP > 10" -r $CHROM:$START-$END $DENI > $DIR/$GENE.Denisova.vcf
fi

## merge Altai, Denisova, Inuit
## select only missense variants
## create a fasta file

# vcf2fasta.R name_gene # name_gene is for naming all result files
NSVCF=`wc -l $DIR/$GENE.vcf | cut -f 1 -d " "`

if [ $NSVCF -gt 253 ];
then
	echo processing...
	Rscript Scripts/vcf2fasta_noInuit.R $GENE $CHROM $START $END

	## create a fasta with unique haplotypes and record their frequency (in traits file)
	Rscript Scripts/fasta2haplos.R Results/Genes/out_sup_$GENE.fasta 0 0
	# Written these files: Results/out_sup_CLTCL1.fasta_0_0.haps.fasta Results/out_sup_CLTCL1.fasta_0_0.traits.txt
	Rscript Scripts/fasta2haplos.R Results/Genes/out_sup_$GENE.fasta 20 0
	# Written these files: Results/out_sup_CLTCL1.fasta_20_0.haps.fasta Results/out_sup_CLTCL1.fasta_20_0.traits.txt
	Rscript Scripts/fasta2haplos.R Results/Genes/out_sup_$GENE.fasta 0 3 # only at least freq=3
	# Written these files: Results/out_sup_CLTCL1.fasta_0_3.haps.fasta Results/out_sup_CLTCL1.fasta_0_3.traits.txt

	Rscript Scripts/fasta2haplos.R Results/Genes/out_pop_$GENE.fasta 0 0

	## write table and differences among haplotypes
	Rscript Scripts/writeTable.R Results/Genes/out_sup_$GENE.fasta_0_0.haps.fasta Results/Genes/out_sup_$GENE.fasta_0_0.traits.txt
	Rscript Scripts/writeTable.R Results/Genes/out_sup_$GENE.fasta_20_0.haps.fasta Results/Genes/out_sup_$GENE.fasta_20_0.traits.txt

	echo calculating...
	Rscript Scripts/haplos2stats.R $GENE $START $END

fi


## get protein position
# cut -d "," -f 3 Results/anno_$GENE.csv 
# then run it on polyphen http://genetics.bwh.harvard.edu/pph2/bgi.shtml
# then save it as: Results/polyphen_CLTCL1.short.txt Results/polyphen_CLTCL1.full.txt

# plot summary stats

#Rscript Scripts/calcPval.R $GENE

# hetero
#Rscript Scripts/calcHetero.R Results/Genes/out_sup_$GENE.fasta > Results/out_sup_$GENE.hetero.csv
#Rscript Scripts/calcHetero.R Results/Genes/out_pop_$GENE.fasta > Results/out_pop_$GENE.hetero.csv

# hetero
#Rscript Scripts/calcHetero_allelic.R Results/Genes/out_sup_$GENE.fasta 9 > Results/out_sup_$GENE.hetero.allelic.csv
#Rscript Scripts/calcHetero_allelic.R Results/Genes/out_pop_$GENE.fasta 9 > Results/out_pop_$GENE.hetero.allelic.csv




