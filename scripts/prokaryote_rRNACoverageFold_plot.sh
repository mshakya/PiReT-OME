#!/bin/sh
##$ -V
## -N prokaryote_rRNACoverageFold_plot_$sample
#$ -l h_vmem=10G
##$ -pe smp 1
##$ -m abe
##$ -M sfeng@lanl.gov
##$ -j y
##$ -o /users/203270/scratch/prokaryote_rRNACoverageFold_plot.log


###parse_BAMfile
echo "perl $scriptDir/prokaryote_rRNACoverageFold_plot.pl  $sample   $workdir"
perl $scriptDir/prokaryote_rRNACoverageFold_plot.pl  $sample   $workdir
