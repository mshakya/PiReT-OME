#!/bin/sh

date
#map to ref
# echo "perl $scriptDir/rRNA_reads_mapping.pl -cpu $numCPU  -p1  $workdir/$sample/trimming_results/$sample.1.trimmed.fastq -p2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq -prefix $sample -index $indexref    -o $workdir"
echo "perl $scriptDir/reads_mapping.pl \
-P1 $workdir/$sample/trimming_results/$sample.1.trimmed.fastq \
-P2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq \
-I $indexref -K $test -E $eukarya_fasta -B $prokaryote_fasta \
-W $workdir"

perl $scriptDir/reads_mapping.pl \
-P1 $workdir/$sample/trimming_results/$sample.1.trimmed.fastq \
-P2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq \
-I $indexref -K $test -E $eukarya_fasta -B $prokaryote_fasta \
-W $workdir
 
#parse BAM files
# echo "perl $scriptDir/parse_BAMfile.pl  -bamfile $rawreads -sample $sample -o $workdir"
# perl $scriptDir/parse_BAMfile.pl  -bamfile $rawreads -sample $sample -o $workdir


###classification

#echo "perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU"
#cd $workdir/$sample/mapping_results/
#perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU  

