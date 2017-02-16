#!/bin/sh
##$ -V
## -N trim_readmapping_$sample
##$ -l h_vmem=10G
##$ -pe smp 10
##$ -m abe
##$ -M sfeng@lanl.gov
##$ -j y
##$ -o /users/203270/scratch/trim_readmapping.log




date
#map to ref
echo "perl $scriptDir/rRNA_reads_mapping.pl -cpu $numCPU  -p1  $workdir/$sample/trimming_results/$sample.1.trimmed.fastq -p2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq -prefix $sample -index $indexref    -o $workdir"
perl $scriptDir/rRNA_reads_mapping.pl -test $test -cpu $numCPU  -p1  $workdir/$sample/trimming_results/$sample.1.trimmed.fastq -p2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq -prefix $sample -index $indexref -o $workdir

###classification

date
#echo "perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU"
#cd $workdir/$sample/mapping_results/
#perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU  

date
