#!/bin/bash
#PBS -N Transcript_Quantification_Salmon
#PBS -l nodes=2:ppn=16
#PBS -l mem=29gb
#PBS -q short
#PBS -j oe
#PBS -o /data1/victorhm/assembly/output/

input=/data1/victorhm/assembly/input
output=/data1/victorhm/assembly/quant_useEM
salmon=/data1/victorhm/assembly/salmon/bin/salmon
ref=/data1/victorhm/assembly/default.Trinity.fasta

$salmon index -t $ref -i ${ref}.index -k 31

cd $input

for R1 in $(ls *R1*)

do
	sample=$(echo $R1 | sed 's/_trim.fastq.gz//')

	sampleR1=$(echo ${R1})

	sampleR2=$(echo ${R1} | sed 's/_R1_trim.fastq.gz/_R2_trim.fastq.gz/')

	$salmon quant -p 30 -i ${ref}.index -l ISR \
	    -1 ${sampleR1} -2 ${sampleR2} \
            -g ${ref}.map \
	    --validateMappings \
	    --gcBias \
	    --useEM \
	    -o ${output}/quants_out/${sample}
done
