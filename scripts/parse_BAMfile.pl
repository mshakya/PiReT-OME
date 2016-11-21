#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$Bin/bin/:$ENV{PATH}";

umask 000;


my ($outDir1,$outDir,$bamfile, $outDir2,  $ref);
my $sample="sample1";
my $outputHost=0;
my $outsam=0;
my $ref;#="/users/218819/scratch/data/databases/human_chromosomes/all_chromosome.fasta";
my $outFasta=0;



GetOptions( 
            'bamfile'=> \$bamfile,
            'sample=s' => \$sample,
            'o=s'   => \$outDir2,
);


if ( ! -e $outDir2) {print "no working dir $outDir2\n";}

$outDir1=join '/', ($outDir2, "$sample");
mkdir $outDir1 if ( ! -e $outDir1);
if ( ! -e $outDir1) {print "can not make dir $outDir1\n";}

$outDir=join '/', ($outDir1, 'mapping_results');
mkdir $outDir if ( ! -e $outDir);
if ( ! -e $outDir) {print "can not make dir $outDir\n";}


#sfeng
my $workdir=$outDir2;
my $headfile="$workdir/prokaryote.fa.fai";
my %seqln;

#unless(&file_check("$outDir2/coverage.fa.fai")<0 ) {$testcoverage=2;}

open (GENOIN, "$headfile") or die "$headfile does not exist $!";
 while (<GENOIN>) {

         chomp;

       my @line = split /\s+/,$_;
     for (my $i=1; $i<=$line[1]; $i++) {
     $seqln{$line[0]}=$line[1];
      }
    }
close GENOIN;


#sfeng

mkdir $outDir if ( ! -e $outDir);

&readBAM($bamfile,$outDir,$sample);

exit(0);


sub readBAM
{
    my $BAMfile=shift;
    my $outputDir=shift;
    my $prefix=shift;


    my $statsfile="$outputDir/$prefix.stats.text";
    my $unalignedNonPairedFile="$outputDir/unmapped.$prefix.unpaired.fastq";
    my $unalignedMate1File = "$outputDir/unmapped.$prefix.1.fastq";
    my $unalignedMate2File = "$outputDir/unmapped.$prefix.2.fastq";
    my $alignedNonPairedFile="$outputDir/unpaired.ref.fastq";

    my $mapfile= "$outputDir/$prefix.error.txt";
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



open (OUTNON,">$workdir/$sample/mapping_results/paired.eukarya_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/paired.eukarya_ref.sam\n";
open (OUTALL,">$workdir/$sample/mapping_results/paired.prokaryote_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/paired.prokaryote_ref.sam\n";
open (BADMAPEU,">$workdir/$sample/mapping_results/Notproperpaired.eukarya_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/Notproperpaired.eukarya_ref.sam\n";
open (BADMAPPRO,">$workdir/$sample/mapping_results/Notproperpaired.prokaryote_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/Notproperpaired.prokaryote_ref.sam\n";
open (OUTFW,">$workdir/$sample/mapping_results/forward.prokaryote_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/forward.prokaryote_ref.sam\n";
open (OUTBW,">$workdir/$sample/mapping_results/backward.prokaryote_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/backward.prokaryote_ref.sam\n";
open (OUTEUFW,">$workdir/$sample/mapping_results/forward.eukarya_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/forward.eukarya_ref.sam\n";
open (OUTEUBW,">$workdir/$sample/mapping_results/backward.eukarya_ref.sam") or die "$! canot open $workdir/$sample/mapping_results/backward.eukarya_ref.sam\n";

#sf
        open (my $unalignedMate1_fh, ">$unalignedMate1File") or die "$! $unalignedMate1File";
        open (my $unalignedMate2_fh, ">$unalignedMate2File") or die "$! $unalignedMate2File";
        open (my $unalignedNonPaired_fh, ">$unalignedNonPairedFile") or die "$! $unalignedNonPairedFile";
       
           my @BAMfile=split /,/, $BAMfile;
          my %mappedreads;
          my %unmappedreads;
         for ( my $i=0; $i <=  $#BAMfile; $i++)
        {
           my  $tmpfile=$BAMfile[$i];
            unless (-e $tmpfile) {next;}


          my $command = "samtools view -h $tmpfile  2>$mapfile | ";
            open (my $fh, "$command") or die "$! $command  failed\n";
            while (<$fh>)
            {
                chomp;
                next if (/^\@/);
                my @samFields=split /\t/,$_;
                my $samline=$_;
                if ($samFields[10] eq "*" and !$outFasta)
                {
                    $samFields[10] = "f" x length($samFields[9]);
                }
                  # bit operation [and] on the flag 
                if (($samFields[1] & 4) and ($samFields[1] & 8))  # both paired reads unmapped
                {
                    if ($samFields[1] & 64) # the read is the first read in a pair
                    {
                        $samFields[0] =~ s/\/\d$//;
                            print $unalignedMate1_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                    if ($samFields[1] & 128) # the read is the second read in a pair
                    {
                        $samFields[0] =~ s/\/\d$//;
                            print $unalignedMate2_fh  "@".$samFields[0]."/2\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                    $unmappedreads{$samFields[0]}=1;
                }
                elsif($samFields[1] & 4)  # query is unmapped
                {
                    if ($samFields[1] & 64)
                    {
                        $samFields[0] =~ s/\/\d$//;
                        if ($outFasta)
                        {
                           # print $unalignedNonPaired_fh  ">".$samFields[0]."/1\n".$samFields[9]."\n";
                        }
                        else
                        {
                            print $unalignedNonPaired_fh  "@".$samFields[0]."/1\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                        }
                    }
                    if ($samFields[1] & 128)
                    {
                        $samFields[0] =~ s/\/\d$//;
                            print $unalignedNonPaired_fh  "@".$samFields[0]."/2\n".$samFields[9]."\n+\n".$samFields[10]."\n";
                    }
                    $unmappedreads{$samFields[0]}=1;
                }
                else  #mapped reads
                {
         if ($seqln{$samFields[2]})
        {
        if ($samFields[1]==99||$samFields[1]==147||$samFields[1]==83||$samFields[1]==163 ) { print OUTALL "$samline\n"; $pair_prokaryote_reads++; $pair_prokaryote_reads{$samFields[2]}++;}#  print OUT "$samline\n";}
                 else {print BADMAPPRO "$samline\n";}
                if ($samFields[1]==99||$samFields[1]==147 )  { print OUTFW "$samline\n";}
        if ($samFields[1]==83||$samFields[1]==163 )  { print OUTBW "$samline\n";}
        } else {
        if ($samFields[1]==99||$samFields[1]==147||$samFields[1]==83||$samFields[1]==163 ) {print OUTNON "$samline\n"; $pair_eukarya_reads++; $pair_eukarya_reads{$samFields[2]}++;} 
        else {print BADMAPEU "$samline\n";}
                 if ($samFields[1]==99||$samFields[1]==147 )  { print OUTEUFW "$samline\n";}
        if ($samFields[1]==83||$samFields[1]==163 ) { print OUTEUBW "$samline\n";}
       }
               $mappedreads{$samFields[0]}=1;
                }
            }
            close $fh;
          
            $numTotalReadsPaired += (keys(%unmappedreads)+keys(%mappedreads))*2;
            $numTotalUnmappedReadsPaired += (keys(%unmappedreads))*2;
        } # end foreach 
        close $unalignedMate1_fh;
        close $unalignedMate2_fh;
        close $unalignedNonPaired_fh;

          close BADMAPEU;
          close BADMAPPRO;
           close OUTFW;
        close OUTBW;
           close OUTEUFW;
        close OUTEUBW;

    close OUTNON;
    close OUTALL; 


my $trimmedreads=$numTotalReadsUnpaired + $numTotalReadsPaired;
my $totalrawreads=$trimmedreads;
 if (-e "$outDir1/trimming_results/$prefix.stats.txt") 
{
$totalrawreads =&parsetrimmingfile("$outDir1/trimming_results/$prefix.stats.txt");
}
my $unmappedreads=$numTotalUnmappedReadsUnpaired + $numTotalUnmappedReadsPaired;
my $mappedreads=$trimmedreads-$unmappedreads; 
open (OUTSTAT,">$statsfile") or die "$! canot open $statsfile\n";
print OUTSTAT "total_reads\t$trimmedreads\n";
print OUTSTAT "total_Unmapped_reads\t$unmappedreads\n";
print OUTSTAT "total_Mapped_reads\t$mappedreads\n";
print OUTSTAT "proper_paired_prokaryote_reads\t$pair_prokaryote_reads\n";
print OUTSTAT "proper_paired_eukarya_reads\t$pair_eukarya_reads\n";
foreach (sort keys %pair_prokaryote_reads) {
print OUTSTAT "proper_paried_reads_in_prokaryote_chromos:\t$_\t$pair_prokaryote_reads{$_}\n";
}
foreach (sort keys %pair_eukarya_reads) {
print OUTSTAT "proper_paried_reads_in_eukarya_chromos:\t$_\t$pair_eukarya_reads{$_}\n";
}

close OUTSTAT;

}

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

sub parsetrimmingfile 

{

 my $log=shift;
    my $numReads=0;
    open (my $fh, $log) or die "$! open $log failed\n";
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

sub executeCommand 
{
    my $command = shift;
    system($command) == 0
         || die "the command $command failed\n";
}

sub file_check 
{
    my $file=shift;
    my $exist=-1;
    if (-e $file) {$exist=1};
    if (-z $file) {$exist=-1};
    return $exist;
}


  
