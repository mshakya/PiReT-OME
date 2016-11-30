#!/bin/sh
##$ -V
## -N prokaryote_find_small_rna
##$ -l h_vmem=10G
##$ -pe smp 1
##$ -m abe
##$ -M sfeng@lanl.gov
##$ -j y
##$ -o /users/203270/scratch/prokaryote_find_small_rna.log



###parse_BAMfile
echo "perl $scriptDir/scripts/prokaryote_find_small_rna.pl $workdir $reffile $gfffile $sample\n"
perl $scriptDir/scripts/prokaryote_find_small_rna.pl $workdir $reffile $gfffile $sample


