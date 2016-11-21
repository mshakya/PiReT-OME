#!/bin/sh
##$ -V
## -N parse_BAMfile_$sample
##$ -l h_vmem=10G
##$ -pe smp 10
##$ -m abe
##$ -M sfeng@lanl.gov
##$ -j y
##$ -o /users/203270/scratch/parse_BAMfile.log



###parse_BAMfile
echo "perl $scriptDir/parse_BAMfile.pl  -bamfile $rawreads -sample $sample    -o $workdir"
perl $scriptDir/parse_BAMfile.pl  -bamfile $rawreads -sample $sample    -o $workdir

###classification

date
#echo "perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU"
#cd $workdir/$sample/mapping_results/
#perl $scriptDir/../../scripts/bwa_sam2classification.pl -p1 unmapped.$sample.1.fastq -p2 unmapped.$sample.2.fastq -u unmapped.$sample.unpaired.fastq -cpu $numCPU  

