#!/bin/sh
##$ -V
## -N trim_readmapping_$sample
##$ -l h_vmem=10G
##$ -pe smp 10
##$ -m abe
##$ -M sfeng@lanl.gov
##$ -j y
##$ -cwd
#TODO: change it to a writable directory
##$ -o /users/203270/scratch/trim_readmapping.log

echo "perl $scriptDir/illumina_fastq_QC.pl  -min_L 60 -n 5 -q 15  $workdir/refseq/artifact.fna -lc 0.7 -t $numCPU  -prefix $sample  -d  $workdir/$sample/trimming_results/ -p  $rawreads"

perl $scriptDir/illumina_fastq_QC.pl  -min_L 60 -n 5 -q 15  $workdir/refseq/artifact.fna -lc 0.7  -t $numCPU  -prefix $sample  -d  $workdir/$sample/trimming_results/ -p  $rawreads

#echo "/opt/apps/bin/perl  $scriptDir/illumina_fastq_QC.pl  -min_L 60 -n 5 -q 15  $workdir/refseq/artifact.fna -lc 0.7 -t $numCPU  -prefix $sample  -d  $workdir/$sample/trimming_results/ -p  $rawreads"

#/opt/apps/bin/perl  $scriptDir/illumina_fastq_QC.pl  -min_L 60 -n 5 -q 15  $workdir/refseq/artifact.fna -lc 0.7  -t $numCPU  -prefix $sample  -d  $workdir/$sample/trimming_results/ -p  $rawreads

#
#date
##map to ref
echo "perl $scriptDir/rRNA_reads_mapping.pl -cpu $numCPU  -p1  $workdir/$sample/trimming_results/$sample.1.trimmed.fastq -p2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq -prefix $sample -index $indexref    -o $workdir"
perl $scriptDir/rRNA_reads_mapping.pl -test $test -cpu $numCPU  -p1  $workdir/$sample/trimming_results/$sample.1.trimmed.fastq -p2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq -prefix $sample -index $indexref -o $workdir


#
####classification
#
#date
##echo "perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU"
##cd $workdir/$sample/mapping_results/
##perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU  
#
#date
