#!/bin/bash

finalDir=RNA_trim/
filesDir=FastQ/RNA_SBC_downsample/
adapters=../../Programs/trimmomatic/adapters/TruSeq3-PE-2.fa
trimmomatic=../../Programs/trimmomatic/trimmomatic-0.38.jar

finalDir=../$finalDir

cd $filesDir

for r1 in *_R1*
do
    r2=`echo "$r1" | sed -e 's/_R1/_R2/g'`
	ident=`echo "$r1" | sed -e 's/_R1.*//'`
	echo ${ident}

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -j 6 \
		-o ${finalDir}"$ident"_R1_temp.fastq.gz -p ${finalDir}"$ident"_R2_temp.fastq.gz \
		$r1 $r2

	java -jar $trimmomatic PE -phred33 -threads 6 \
		${finalDir}"$ident"_R1_temp.fastq.gz \
		${finalDir}"$ident"_R2_temp.fastq.gz \
                ${finalDir}"$ident"_R1_trim.fastq.gz ${finalDir}"$ident"_forward_unpaired.fastq.gz \
		${finalDir}"$ident"_R2_trim.fastq.gz ${finalDir}"$ident"_reverse_unpaired.fastq.gz \
		HEADCROP:13 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50 

	rm ${finalDir}"$ident"_R1_temp.fastq.gz
	rm ${finalDir}"$ident"_R2_temp.fastq.gz
done
