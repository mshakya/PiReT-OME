#!/code/bin/sh
## -N htseq_$sample
##$ -V
##$ -l h_vmem=20G
##$ -pe smp 1
##$ -m abe
##$ -M sfeng@lanl.gov
##$ -j y
##$ -o /users/203270/scratch/htseq.log


echo"$scriptDir/htseq-count.pl $workdir $sample $test\n"

$scriptDir/htseq-count.pl $workdir $sample $test


