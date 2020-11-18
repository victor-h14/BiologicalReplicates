#!/bin/bash
#PBS -N signalp_script
#PBS -l nodes=1:ppn=16
#PBS -l mem=10gb
#PBS -q short
#PBS -j oe
#PBS -o log.out

signalp/signalp-5.0/bin/signalp -fasta analysis/default.Trinity.fasta.transdecoder.pep -gff3 > analysis/signalp.out

