#!/bin/bash
#PBS -N blastp_job
#PBS -l nodes=1:ppn=15
#PBS -q short
#PBS -j oe
#PBS -o log.out

module load blast/2.6.0

blastp 	-db Swissprot/swissprot_db \
	-query analysis/default.Trinity.fasta.transdecoder.pep \
	-max_target_seqs 20 -outfmt 6 -evalue 1e-5 -num_threads 15 > analysis/blastp_out.outfmt6
