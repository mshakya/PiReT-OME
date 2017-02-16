#!/bin/sh


#QC step
echo "perl $scriptDir/illumina_fastq_QC.pl  -min_L 60 -n 5 -q 15  $workdir/refseq/artifact.fna -lc 0.7 -t $numCPU  -prefix $sample  -d  $workdir/$sample/trimming_results/ -p  $rawreads"
perl $scriptDir/illumina_fastq_QC.pl  -min_L 60 -n 5 -q 15  $workdir/refseq/artifact.fna -lc 0.7  -t $numCPU  -prefix $sample  -d  $workdir/$sample/trimming_results/ -p  $rawreads

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
