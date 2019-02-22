        
#NR=`cat mytmp.bim | wc -l`

NR=121

> mytmp.anc

for N in `seq 1 1 $NR`
do
	#echo $N
	head -n $N mytmp.bim | tail -n 1 > tmp
	CHROM=`cut -f 1 tmp`
	POS=`cut -f 4 tmp`
	A1=`cut -f 5 tmp`
	A2=`cut -f 6 tmp`
	rm tmp

	#echo $A1
	#echo $A2

	ANC=`~/Software/samtools/samtools faidx ~/Data/Human/human_ancestor_$CHROM.fa $CHROM:$POS-$POS | head -n 2 | tail -n 1 | tr "[:lower:]" "[:upper:]"`

	#echo $ANC

	if [ ! "$ANC" == "$A1" ]; then
		if [ ! "$ANC" == "$A2" ]; then	 

			if [ "$ANC" == "A" ]; then
				ANC="T"
			fi

                        if [ "$ANC" == "T" ]; then
                                ANC="A"
                        fi

			if [ "$ANC" == "G" ]; then
                                ANC="C"
                        fi

			if [ "$ANC" == "C" ]; then
                        	ANC="G"
                        fi

		fi
	fi
	
	if [ ! "$ANC" == "$A1" ]; then
                if [ ! "$ANC" == "$A2" ]; then
			ANC="N"
		fi
	fi

	echo "$ANC" >> mytmp.anc

done



