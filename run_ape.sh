
# on hg18!

GENE=$1
CHROM=$2
START=$3
END=$4

echo $GENE $CHROM $START $END

LIB=~/Software/vcflib/bin
DIR=~/Data/Ape
SHAPEIT=~/Software/shapeit/bin/shapeit
GMAP=~/Data/HapMap/Recomb/genetic_map_chr${CHROM}_b36.txt

echo filtering...
if [ ! -f ~/Data/CLTCL1/Apes/$GENE.vcf ];
then

	$LIB/vcfcombine $DIR/Gorilla.vcf.gz $DIR/Pan_paniscus.vcf.gz $DIR/Pan_troglodytes.vcf.gz $DIR/Pongo_abelii.vcf.gz $DIR/Pongo_pygmaeus.vcf.gz -r chr$CHROM:$START-$END > $DIR/tmp.vcf
	$LIB/vcffilter -f "AC > 1" -f "QUAL > 32" -f "FS < 27" -f "DP > 50" -f "DP < 7000" -f "MQ > 25" $DIR/tmp.vcf > $DIR/tmp2.vcf
	Rscript Scripts/fillHomoRef.R $DIR/tmp2.vcf 36 > ~/Data/CLTCL1/Apes/$GENE.vcf
	\rm $DIR/tmp.vcf $DIR/tmp2.vcf

fi

echo phasing...
$SHAPEIT -phase --input-vcf ~/Data/CLTCL1/Apes/$GENE.vcf -O ~/Data/CLTCL1/Apes/$GENE.phased --thread 2 --burn 50 --prune 20 --main 100 --window 0.5 --effective-size 20000
$SHAPEIT -convert --input-haps ~/Data/CLTCL1/Apes/$GENE.phased --output-vcf ~/Data/CLTCL1/Apes/$GENE.phased.vcf &> /dev/null
\rm shapeit_*

echo protein...
Rscript Scripts/makeProteinApe.R $GENE
#less -S Results/Apes/$GENE.anno.txt

echo fasta...
Rscript Scripts/createFastaApe.R $GENE

echo haplotypes...
Rscript Scripts/fasta2haplos_apes.R Results/Apes/$GENE.pop.fasta 0 0 0
Rscript Scripts/fasta2haplos_apes.R Results/Apes/$GENE.pop.fasta 0 0 1

echo tables...
Rscript Scripts/writeTable.R Results/Apes/$GENE.pop.fasta_0_0.haps.fasta Results/Apes/$GENE.pop.fasta_0_0.traits.txt



