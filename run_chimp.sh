
# on hg19!

GENE=$1
CHROM=$2
START=$3
END=$4

echo $GENE $CHROM $START $END

DIR=~/Data/Chimp/hg19
LIB=~/Software/vcflib/bin
DATA=$DIR/Freebayes_hg19_Concat_1Mb_Clean_Callable_HWE.vcf.gz

SHAPEIT=~/Software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit
GMAP=~/Data/HapMap/Recomb/hg19/genetic_map_GRCh37_chr${CHROM}.fixed.txt

echo filtering...
if [ ! -f ~/Data/CLTCL1/Chimp/$GENE.vcf ];
then

	$LIB/vcffilter -f "AC > 1" -f "QUAL > 32" -f "DP > 50" -f "DP < 7000" -r chr$CHROM:$START-$END $DATA > $DIR/tmp3.vcf
        Rscript Scripts/fillHomoRef.R $DIR/tmp3.vcf 57 > ~/Data/CLTCL1/Chimp/$GENE.vcf
        \rm $DIR/tmp3.vcf

fi

echo phasing...

$SHAPEIT -phase --input-vcf ~/Data/CLTCL1/Chimp/$GENE.vcf -O ~/Data/CLTCL1/Chimp/$GENE.phased --thread 2 --burn 50 --prune 20 --main 100 --window 0.5 --effective-size 20000 &> /dev/null

$SHAPEIT -convert --input-haps ~/Data/CLTCL1/Chimp/$GENE.phased --output-vcf ~/Data/CLTCL1/Chimp/$GENE.phased.vcf &> /dev/null
\rm shapeit_*

echo protein...
Rscript Scripts/makeProteinChimp.R $GENE
echo `wc -l Results/Chimp/$GENE.anno.txt`

echo fasta...
Rscript Scripts/createFastaChimp.R $GENE

echo haplotypes...
# it's called _apes but it is valid for chimp too
Rscript Scripts/fasta2haplos_apes.R Results/Chimp/$GENE.pop.fasta 0 0 0
Rscript Scripts/fasta2haplos_apes.R Results/Chimp/$GENE.pop.fasta 0 0 1

echo tables...
Rscript Scripts/writeTable.R Results/Chimp/$GENE.pop.fasta_0_0.haps.fasta Results/Chimp/$GENE.pop.fasta_0_0.traits.txt






