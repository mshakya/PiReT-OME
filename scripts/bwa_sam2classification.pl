#! /usr/bin/env perl -w
use strict;

use FindBin qw($Bin);
use POSIX qw(strftime);


$|=1;
$ENV{PATH} = "$Bin";
use lib "$Bin";

use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use gi2lineage;

my $outfile='unmapped.reads.classification.txt';
my $outfile1='unmapped.ginumber.txt';
my ($IndexFile, $bowti2opt, $queryPairedFile2, $queryPairedFile1, $unpairedFile, $numCPU);
my %gi;

GetOptions(
              'p1=s'   => \$queryPairedFile1,
              'p2=s'   => \$queryPairedFile2,
            'u=s'   => \$unpairedFile,
             'cpu=i' => \$numCPU,
             'index_ref_bt2=s' => \$IndexFile,
             'btopt=s' => \$bowti2opt,
            'help|?'   => sub{&Usage()}
);


$IndexFile='/users/203270/scratch/bowtie2_db/bt2_NCBI-Bacteria-Virus';
$numCPU=1;
$bowti2opt='--very-fast';

unless ($queryPairedFile1 && $queryPairedFile2  )  {&Usage();}

sub Usage
{
     print <<"END";
 Usage: perl $0 [options] bwa_sam2classification.pl  -index_ref_bt2 indexfile -p1 queryPairedFile1 -p2 queryPairedFile2 -U unpairedFile -cpu numCPU -btopt optbowti2option 
                 
 example: perl ./bwa_sam2classification.pl  -p1 unmapped.IS_24h_3.1.fastq -p2 unmapped.IS_24h_3.2.fastq -u unmapped.IS_24h_3.unpaired.fastq -cpu 10 

END
exit;
}


if(&file_check($queryPairedFile1)<0 ) { die "failed: The file $queryPairedFile1 doesn't exist or empty.\n";}
if(&file_check($queryPairedFile2)<0 ) { die "failed: The file $queryPairedFile2 doesn't exist or empty.\n";}
my $tmpcheckfile=join "", ($IndexFile, '.2.ht2l');
my $tmpcheckfile1=join "", ($IndexFile, '.2.bt2');
if(&file_check($tmpcheckfile)<0 && &file_check($tmpcheckfile1)<0) { die "failed: The index file $IndexFile doesn't exist or empty.\n";}
my $command;

if ( $unpairedFile) 
 {
 if(&file_check($unpairedFile)<0 ) { die "failed: The file $unpairedFile doesn't exist or empty.\n";} 
  else {
$command = "bowtie2  $bowti2opt  -p $numCPU -x $IndexFile -1 $queryPairedFile1 -2 $queryPairedFile2 -U $unpairedFile  2> classificationLogFile| ";
     }
 } else {
 $command = "bowtie2  $bowti2opt  -p $numCPU -x $IndexFile -1 $queryPairedFile1 -2 $queryPairedFile2   2> classificationLogFile| ";
 } 

my @rank=qw(phylum class order family genus species);
#loadTaxonomy("preload");

open (FILEIN, "$command") or die "failed: $! $command failed\n";
while (<FILEIN>)
 {	chomp;
	my @fields = split /\t/, $_;
	if( $#fields >5 ){
        my ($gi) = $fields[2] =~ /gi\|(\d+)/;
        if ($gi) {
        $gi{$gi}++;
       }
    }
 }
close FILEIN;

open (FILEOUT, ">$outfile") or die "failed: $! failed\n";
 print FILEOUT "#reads";
foreach (my $i=$#rank; $i>=0; $i--)     
 {
 print  FILEOUT "\t$rank[$i]";
 }
 print FILEOUT"\n";


open (FILEOUT1, ">$outfile1") or die "failed: $! $outfile1 failed\n"; 
my %lineage;
foreach (sort {$gi{$b} <=> $gi{$a}} keys %gi ) 
 {
 my $gi=$_;
print FILEOUT1"$gi\n";
 my @lineage;
foreach (my $i=$#rank; $i>=0; $i--) 
 {
 my $name = gi2rank($gi,$rank[$i]);
  if ($name) 
  {
 push @lineage, $name;
  }
 }
 my $lineage=join "\t", @lineage;
 $lineage{$lineage}+=$gi{$gi};
}


 foreach (sort {$lineage{$b} <=> $lineage{$a}} keys %lineage )
 {
 if ($lineage{$_}>=100) 
   {
 print FILEOUT"$lineage{$_}\t$_\n";
   }
  }

close FILEOUT;
close FILEOUT1;


sub file_check 
{
    #check file exist and non zero size
    my $file=shift;
    my $exist=-1;
    if (-e $file) {$exist=1};
    if (-z $file) {$exist=-1};
    return $exist;
}

