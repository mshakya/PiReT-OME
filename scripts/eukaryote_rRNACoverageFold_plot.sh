#!/bin/sh
#$ -V
## -N eukaryote_rRNACoverageFold_$sample
##$ -l h_vmem=10G
##$ -pe smp 1
##$ -m abe
##$ -M sfeng@lanl.gov
##$ -j y
##$ -o /users/203270/scratch/eukaryote_rRNACoverageFold.log



###parse_BAMfile
echo "perl $scriptDir/scripts/eukaryote_rRNACoverageFold_plot.pl  $sample   $workdir"
perl $scriptDir/scripts/eukaryote_rRNACoverageFold_plot.pl  $sample   $workdir

