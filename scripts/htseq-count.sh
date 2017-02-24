#!/bin/sh


##$PBS_O_WORKDIR=$workdir
##$PBS_O_PATH=$path

#cd $PBS_O_WORKDIR
echo "perl $scriptDir/htseq-count.pl $workdir $sample $test\n"
perl $scriptDir/htseq-count.pl $workdir $sample $test
