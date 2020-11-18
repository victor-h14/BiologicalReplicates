#!/bin/bash
#PBS -N trinity_job
#PBS -l nodes=1:ppn=15
#PBS -q bigmem
#PBS -j oe
#PBS -o assembly/output/error.out

#Load bowtie2 and samtools
module load bowtie2/2.3.2
module load samtools/1.5

#export environment path
export PATH=${PATH}:assembly/jellyfish-2.2.10/bin
export PATH=${PATH}:assembly/salmon-0.12.0_linux_x86_64/bin
export LD_LIBRARY_PATH=${PATH}:assembly/jellyfish-2.2.10/lib
export MANPATH=${PATH}:assembly/jellyfish-2.2.10/share/man
export PKG_CONFIG_PATH=${PATH}:assembly/jellyfish-2.2.10/lib/pkgconfig


#Trinity - Local directory:
Trinity=assembly/trinity/Trinity    
Left=`echo assembly/input/*_R1_trim.fastq.gz | tr ' ' ,`
Right=`echo assembly/input/*_R2_trim.fastq.gz | tr ' ' ,`

$Trinity --seqType fq --max_memory 200G \
    --left $Left \
    --right $Right \
    --CPU 15 \
    --SS_lib_type RF \
    --full_cleanup \
    --normalize_by_read_set --verbose \
    --output assembly/output-trinity/
