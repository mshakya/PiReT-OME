#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use FindBin qw($Bin);


#TODO: test and see if this ENV can be removed since it is set on runPipeline
$ENV{PATH} = "$Bin:$Bin/bin/:$ENV{PATH}";

umask 000;


my ($pairedReadsFile1, $indexFile, $pairedReadsFile2, $unpairedReadsFile, $workdir, $sampleDir, $mapDir,  $ref);
my $Bowtie2Opts="--fast"; 
my $numCPU=4; #number of threads [4]
my $prefix="Unmapped_reads";
my $outputHost=0;
my $outsam=0;
my $ref;#="/users/218819/scratch/data/databases/human_chromosomes/all_chromosome.fasta";
my $outFasta=0;
my $test='prokaryote';
my $htseq='gene';

my $testcoverage=2;


GetOptions( 
			'p1=s'   	=> \$pairedReadsFile1, 
            'p2=s'   	=> \$pairedReadsFile2,
            'u=s'   	=> \$unpairedReadsFile,
            'index=s' 	=> \$indexFile,
            'prefix=s' 	=> \$prefix,
            'test=s' 	=> \$test,
            'host'  	=> \$outputHost,
            'sam' 		=> \$outsam,
            'fasta' 	=> \$outFasta,
            'cpu=i'   	=> \$numCPU, # bwa option
            'Bowtie2Opts=s'   => \$Bowtie2Opts,    # bwa mem options
            'o=s'   	=> \$workdir,
#            'geneopt=s' =>\$htseq,       #count reads based on 'gene' or 'CDS' or 'tRNA' or 'mRNA' in annotation file, default ='gene';
            'help|?'   	=> sub{&Usage(1)}
);


# First, check if the working directory exists
if ( ! -e $workdir) {print "no working dir $workdir found\n";}

# Second, make directories based on sample name within working dir
$sampleDir=join '/', ($workdir, "$prefix");
mkdir $sampleDir if ( ! -e $sampleDir);

# Third, throw, a complain, if cant make directory
if ( ! -e $sampleDir) {print "cannot make dir $sampleDir\n";}

# Fourth, make directory for adding the mapping results
$mapDir=join '/', ($sampleDir, 'mapping_results');
mkdir $mapDir if ( ! -e $mapDir);
if ( ! -e $mapDir) {print "cannot make dir $mapDir\n";}


my $headfile="$workdir/prokaryote.fa.fai";
my %seqln;

#unless(&file_check("$outDir2/coverage.fa.fai")<0 ) {$testcoverage=2;}


# check, if the prokaryotic index file exist,
#TODO: Ask Shihai, how is this affected when eukaryotic option is chosen
open (GENOIN, "$headfile") or die "$headfile does not exist $!";

# Process
 while (<GENOIN>) {

         chomp;

       my @line = split /\s+/,$_;
     for (my $i=1; $i<=$line[1]; $i++) {
     $seqln{$line[0]}=$line[1];
      }
    }
close GENOIN;

# reassining prefix to sample variables
my $sample=$prefix;


# Run the subroutine checkfiles
&checkFiles($indexFile, $pairedReadsFile1,$pairedReadsFile2,$unpairedReadsFile);


# subroutine for Usage of this script
sub Usage
{
     my $Bowtie2op=shift;
     print <<"END";
 Usage: perl $0 [options] -p1 reads1.fastq -p2 reads2.fastq -u reads.fa -ref reference.fa -index reference  -o out_directory
        Input File:
        -ref          reference sequences in fasta  

        -index             Index filename prefix (minus trailing .X.bt2
                      NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.

        -u            Unpaired reads, Comma-separated Single end reads  
       
        -p1            leftSequenceFile1,leftSequenceFile2 Comma-separated Left-end reads  

        -p2            rightSequenceFile1,rightSequenceFile2 Comma-separated right-end reads 

        -test 
        Output:
        -o            Output directory.
 
        -prefix       Output File Name Prefix [ref_clean] 
 
        -fasta        <boolean> Output in fasta format instead of fastq.
 
        -host         <boolean> Output Host id list.   

        -sam          <boolean>Output sam file
  
        Options:
        -Bowtie2Opts   <String in quote> see "Bowtie2 options" 
                         ex: '--sensitive-local' 
   
        -cpu          number of CPUs [4]
 
        -h            print "Bowtie2 " options
END
 
if ($Bowtie2op)
{ 
print <<"END";

Bowtie2 options:

 Input:
  -q                 query input files are FASTQ .fq/.fastq (default)
  --qseq             query input files are in Illumina's qseq format
  -f                 query input files are (multi-)FASTA .fa/.mfa
  -r                 query input files are raw one-sequence-per-line
  -c                 <m1>, <m2>, <r> are sequences themselves, not files
  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
  -u/--upto <int>    stop after first <int> reads/pairs (no limit)
  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
  --phred33          qualities are Phred+33 (default)
  --phred64          qualities are Phred+64
  --int-quals        qualities encoded as space-delimited integers

 Presets:                 Same as:
  For --end-to-end:
   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50
   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50
   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 
   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

  For --local:
   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75 (default)
   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 
   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50


END
 }
exit;
}


# create the outDir,
#TODO: why is working directory being created here?
# mkdir $outDir if ( ! -e $outDir);


#run the mapping subroutine and exit
&runMapping($indexFile,$pairedReadsFile1,$pairedReadsFile2,$unpairedReadsFile,$mapDir);
print "mapping is run and donw";
exit(0);




# subroutine for checking if all files exist
sub checkFiles
{
    my $indexFile_r=shift;
    my $queryPairedFile_r1=shift;
    my $queryPairedFile_r2=shift;  
    my $queryUnpairedFile=shift;
    my @queryPairedFile1 = split /\,/, $queryPairedFile_r1;
    my @queryPairedFile2 = split /\,/, $queryPairedFile_r2;
    my @queryUnpairedFile = split /\,/, $queryUnpairedFile;

              $indexFile_r .= ".1.ht2l";
            if ( -z $indexFile_r || ! -e $indexFile_r ) { die "Index file $indexFile_r file is empty or does not exist";}

     if ($#queryPairedFile1 != $#queryPairedFile2 ) { die "mate1 and mate2 have diffrent number of files ";}
    foreach my $queryPairedFile1 (@queryPairedFile1)
    {
        if (defined $queryPairedFile1)
        {
            if ( -z $queryPairedFile1 || ! -e $queryPairedFile1 ) { die "mate1 $queryPairedFile1 file is empty";}
        }
    }

     foreach my $queryPairedFile2 (@queryPairedFile2)
    {
        if (defined $queryPairedFile2)
        {
            if ( -z $queryPairedFile2 || ! -e $queryPairedFile2 ) { die "mate2 $queryPairedFile2 file is empty";}
        }
    }

    

    foreach my $queryUnpairedFile (@queryUnpairedFile)
    {
        if (defined $queryUnpairedFile)
        {
            if ( -z $queryUnpairedFile || ! -e $queryUnpairedFile) { die "$queryUnpairedFile file is empty";}
        }
    } 
}

##############################################################################
#subroutine for running mapping
sub runMapping 
{
	# inputs
    my $IndexFile=shift;
    my $queryPairedFile_r1=shift;
    my $queryPairedFile_r2=shift;
    my $queryUnpairedFile=shift;	
	#outputs
    my $outputDir=shift; #output directory
    my $mappingLogFile="$outputDir/$prefix.mapping.log";
    my $statsfile="$outputDir/$prefix.stats.text";
    my $unalignedNonPairedFile="$outputDir/unmapped.$prefix.unpaired.fastq";
    my $unalignedMate1File = "$outputDir/unmapped.$prefix.1.fastq";
    my $unalignedMate2File = "$outputDir/unmapped.$prefix.2.fastq";
    my $alignedNonPairedFile="$outputDir/unpaired.ref.fastq";
    my $refIdList= "$outputDir/$prefix.refId.txt";
    my $numUnmappedPairedFile=0; # non ref reads number
    my $numUnmappedUnpairedFile=0; # non ref reads number
    my $numTotalReadsPaired=0;
    my $numTotalUnmappedReadsPaired=0;
    my $numTotalReadsUnpaired=0;
    my $numTotalUnmappedReadsUnpaired=0;
	#sf
	my $numMapped=0;
	my %eukarya_reads;
	my %prokaryote_reads;
	my %pair_prokaryote_reads;
	my %pair_eukarya_reads;
	my $eukarya_reads=0;
	my $prokaryote_reads=0;
	my $pair_prokaryote_reads=0;
	my $pair_eukarya_reads=0;
	#outputs
	open (OUTNON,">$workdir/$sample/mapping_results/paired.eukarya_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/paired.eukarya_ref.sam\n";
	open (OUTALL,">$workdir/$sample/mapping_results/paired.prokaryote_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/paired.prokaryote_ref.sam\n";
	open (BADMAPEU,">$workdir/$sample/mapping_results/Notproperpaired.eukarya_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/Notproperpaired.eukarya_ref.sam\n";
	open (BADMAPPRO,">$workdir/$sample/mapping_results/Notproperpaired.prokaryote_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/Notproperpaired.prokaryote_ref.sam\n";
	open (OUTFW,">$workdir/$sample/mapping_results/forward.prokaryote_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/forward.prokaryote_ref.sam\n";
	open (OUTBW,">$workdir/$sample/mapping_results/backward.prokaryote_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/backward.prokaryote_ref.sam\n";
	open (OUTEUFW,">$workdir/$sample/mapping_results/forward.eukarya_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/forward.eukarya_ref.sam\n";
	open (OUTEUBW,">$workdir/$sample/mapping_results/backward.eukarya_ref.sam") or die "$! cannot open $workdir/$sample/mapping_results/backward.eukarya_ref.sam\n";
	#unlinked unPaired files
    unlink $unalignedNonPairedFile if ( -s $unalignedNonPairedFile);
    #print in screen the status
	print colored ("Running reads mapping to reference sequence ... and log file: $mappingLogFile",'yellow'),"\n";
	# if reads are paired
	if ($queryPairedFile_r1 && $queryPairedFile_r2)
    	{
        open (my $unalignedMate1_fh, ">$unalignedMate1File") or die "$! $unalignedMate1File";
        open (my $unalignedMate2_fh, ">$unalignedMate2File") or die "$! $unalignedMate2File";
        open (my $unalignedNonPaired_fh, ">$unalignedNonPairedFile") or die "$! $unalignedNonPairedFile";
        open (my $refId_fh, ">$refIdList") or die "$refIdList $!\n" if ($outputHost);
     	# spliting paired files
       	my @queryPairedFile1 = split /\,/, $queryPairedFile_r1;
       	my @queryPairedFile2 = split /\,/, $queryPairedFile_r2;
		# loop through query paired file
        for ( my $i=0; $i <=  $#queryPairedFile1; $i++)
        	{
           	my $command; 
           	my $queryPairedFile1=$queryPairedFile1[$i];
           	my $queryPairedFile2=$queryPairedFile2[$i];
			# if splice site file is present
         	if (-e "$workdir/differential_gene/eukarya/splice_sites_gff.txt") 
           	{
				&executeCommand("hisat2 $Bowtie2Opts --known-splicesite-infile $workdir/differential_gene/eukarya/splice_sites_gff.txt -p $numCPU -x $IndexFile -1 $queryPairedFile1 -2 $queryPairedFile2  2>$mappingLogFile > $workdir/$sample/mapping_results/mapped.sam");
			}
			else {
            	$command = "hisat2 $Bowtie2Opts -p $numCPU -x $IndexFile -1 $queryPairedFile1 -2 $queryPairedFile2  2>$mappingLogFile > $workdir/$sample/mapping_results/mapped.sam "; 
           		&executeCommand($command);
			}
			# open more files
			open (my $fh, "$workdir/$sample/mapping_results/mapped.sam" ) or die "$! hisat2 $Bowtie2Opts failed\n";
            while (<$fh>)
            {
                chomp;
                next if (/^\@/);
                my @samFields=split /\t/,$_;
                my $samline = $_;
                if ($samFields[10] eq "*" and !$outFasta) {
                    $samFields[10] = "f" x length($samFields[9]);
                }
                # bit operation [and] on the flag 
                if (($samFields[1] & 4) and ($samFields[1] & 8)) { 
					# both paired reads unmapped
                    if ($samFields[1] & 64) {
						# the read is the first read in a pair
                        $samFields[0] =~ s/\/\d$//;
                            print $unalignedMate1_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                    if ($samFields[1] & 128) {
						# the read is the second read in a pair
                        $samFields[0] =~ s/\/\d$//;
                        print $unalignedMate2_fh  "@".$samFields[0]."/2\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                    $numUnmappedPairedFile++;
                }
                elsif($samFields[1] & 4)  # query is unmapped
                {
                    if ($samFields[1] & 64)
                    {
                        $samFields[0] =~ s/\/\d$//;
                        if ($outFasta){
                           # print $unalignedNonPaired_fh  ">".$samFields[0]."/1\n".$samFields[9]."\n";
                        }
                        else {
                            print $unalignedNonPaired_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                        }
                    }
                    if ($samFields[1] & 128){
                        $samFields[0] =~ s/\/\d$//;
                            print $unalignedNonPaired_fh  "@".$samFields[0]."/2\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                    $numUnmappedPairedFile++;
                }
                else  #mapped reads
                {

         			if ($seqln{$samFields[2]}){
        				if ($samFields[1]==99||$samFields[1]==147||$samFields[1]==83||$samFields[1]==163 ) { print OUTALL "$samline\n"; $pair_prokaryote_reads++; $pair_prokaryote_reads{$samFields[2]}++;}#  print OUT "$samline\n";}
                 		else {print BADMAPPRO "$samline\n";}
                		if ($samFields[1]==99||$samFields[1]==147 ) { if ($testcoverage==2) { print OUTFW "$samline\n";}}
        				if ($samFields[1]==83||$samFields[1]==163 ) { if ($testcoverage==2) { print OUTBW "$samline\n";}}
        			}
					else {
        				if ($samFields[1]==99||$samFields[1]==147||$samFields[1]==83||$samFields[1]==163 ) {print OUTNON "$samline\n"; $pair_eukarya_reads++; $pair_eukarya_reads{$samFields[2]}++;} 
        				else {print BADMAPEU "$samline\n";}
         		if ($samFields[1]==99||$samFields[1]==147 ) { if ($testcoverage==2) { print OUTEUFW "$samline\n";}}
        if ($samFields[1]==83||$samFields[1]==163 ) { if ($testcoverage==2) { print OUTEUBW "$samline\n";}}
       }

      }
     }
            close $fh;
            my @statReadspaired = &parseMappingLog($mappingLogFile);
            $numTotalReadsPaired += $statReadspaired[0];
            $numTotalUnmappedReadsPaired += $statReadspaired[1];
            die "\nERROR: $queryPairedFile1 and $queryPairedFile2 Paired reads have different names.\n" if (`grep "paired reads have different names" $mappingLogFile`);
        } # end foreach 
	close $unalignedMate1_fh;
    close $unalignedMate2_fh;
    close $unalignedNonPaired_fh;
    close BADMAPEU;
    close BADMAPPRO;
    }

   	close OUTFW;
    close OUTBW;
    close OUTEUFW;
    close OUTEUBW;
    close OUTNON;
    close OUTALL; 

	my $trimmedreads = $numTotalReadsUnpaired + $numTotalReadsPaired*2;
	my $totalrawreads = $trimmedreads;

 	if (-e "$sampleDir/trimming_results/$prefix.stats.txt"){
		$totalrawreads = &parsetrimmingfile("$sampleDir/trimming_results/$prefix.stats.txt");
	}	

	my $unmappedreads = $numTotalUnmappedReadsUnpaired + $numTotalUnmappedReadsPaired*2;
	my $mappedreads = $trimmedreads-$unmappedreads; 
	open (OUTSTAT,">$statsfile") or die "$! cannot open $statsfile\n";
	print OUTSTAT "total_reads\t$trimmedreads\n";
	print OUTSTAT "total_Unmapped_reads\t$unmappedreads\n";
	print OUTSTAT "total_Mapped_reads\t$mappedreads\n";
	print OUTSTAT "proper_paired_prokaryote_reads\t$pair_prokaryote_reads\n";
	print OUTSTAT "proper_paired_eukarya_reads\t$pair_eukarya_reads\n";
	foreach (sort keys %pair_prokaryote_reads) {
		print OUTSTAT "proper_paried_reads_in_prokaryote_chromos:\t$_\t$pair_prokaryote_reads{$_}\n";
	}
	foreach (sort keys %pair_eukarya_reads) {
		print OUTSTAT "proper_paried_reads_in_eukarya_chromos:\t$_\t$pair_eukarya_reads{$_}\n"}
	close OUTSTAT;
}
############################################################
sub parseMappingLog 
{
    my $log=shift;
    my $numReads=0;
    my $unmapped=0;
    open (my $fh, $log) or die "$! open $log failed\n";
    while (<$fh>)
    {  
        if ($_=~ /reads\; of these\:/)
        {
         my  @tmpline = split /\s+/, $_  ;
           $numReads=$tmpline[0];

        }
        if ($_=~ /aligned 0 times/)
        {
           my @tmpline=split /\(/, $_;
           $unmapped=$tmpline[0];
           $unmapped=~s/\s+//g;
        last;
        }

    }
    close $fh;
    return ($numReads, $unmapped);
}

############################################################
sub parsetrimmingfile 

{

 my $log=shift;
    my $numReads=0;
    open (my $fh, $log) or "$! open $log failed\n";
    while (<$fh>)
    {
        if ($_=~ /Reads \#\: /)
        {
         my  @tmpline = split /\s+/, $_  ;
           $numReads=$tmpline[-1];
         last;
        }

    }
    close $fh;
    return $numReads;
}

################################################################
sub executeCommand 
{
    my $command = shift;
	print $command;
	`$command`;
}
################################################################
sub file_check 
{
    my $file=shift;
    my $exist=-1;
    if (-e $file) {$exist=1};
    if (-z $file) {$exist=-1};
    return $exist;
}
 
