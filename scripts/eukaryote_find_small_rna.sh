#!/bin/sh
#$ -V
## -N eukaryote_find_small_rna.sh
##$ -l h_vmem=10G
##$ -pe smp 1
#$ -m abe





###parse_BAMfile
echo "perl $scriptDir/scripts/eukaryote_find_small_rna.pl $workdir $reffile $gfffile $sample\n"
perl $scriptDir/scripts/eukaryote_find_small_rna.pl $workdir $reffile  $gfffile $sample

