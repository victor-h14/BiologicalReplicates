#!/bin/bash
#PBS -N hmmer_script
#PBS -l nodes=1:ppn=16
#PBS -l mem=10gb
#PBS -q short
#PBS -j oe
#PBS -o log.out

hmmscan --cpu 16 --domtblout analysis/hmmer.out \
	analysis/Pfam-A.hmm \
	analysis/default.Trinity.fasta.transdecoder.pep >analysis/hmmer.log
