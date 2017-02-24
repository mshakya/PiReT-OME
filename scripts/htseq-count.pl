#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use FindBin qw($Bin);

$ENV{PATH} = "$Bin/../bin/:$ENV{PATH}";


umask 000;

my $workdir=$ARGV[0];
my $sample=$ARGV[1];
my $test=$ARGV[2];

my $process_log_file="$workdir/process.log";
open ( LOG, ">>", $process_log_file) or die "failed: failed to write $process_log_file\n$!";

if (-e "$workdir/$sample/mapping_results/finish_htseq-count.txt") { `rm $workdir/$sample/mapping_results/finish_htseq-count.txt`;}

if ( $test eq 'both' || $test eq 'prokaryote')
{
if (-e "$workdir/$sample/mapping_results/paired.prokaryote_ref.sam" ) 
{
`samtools view -bt $workdir/prokaryote.fa.fai  $workdir/$sample/mapping_results/paired.prokaryote_ref.sam  > $workdir/$sample/mapping_results/paired.prokaryote_ref.bam`;
# unlink "$workdir/$sample/mapping_results/paired.prokaryote_ref.sam" ;
}

if( -e "$workdir/$sample/mapping_results/paired.prokaryote_ref.bam")
{
`samtools sort -n $workdir/$sample/mapping_results/paired.prokaryote_ref.bam -o $workdir/$sample/mapping_results/sortedname.paired.prokaryote_ref.bam`;
# unlink "$workdir/$sample/mapping_results/paired.prokaryote_ref.bam";
}

`samtools view -h $workdir/$sample/mapping_results/sortedname.paired.prokaryote_ref.bam > $workdir/$sample/mapping_results/sortedname.paired.prokaryote_ref.sam`;

 opendir(DIR, "$workdir/differential_gene/prokaryote/") or die $!;
 while (my $tmpdir = readdir(DIR))
   {
       next unless -d "$workdir/differential_gene/prokaryote/$tmpdir";
       next if ($tmpdir =~ /^\./);
       my $tmpoutfile="$workdir/sum_gene_count/tmp_count/prokaryote/$tmpdir/$sample.prokaryote.name.htseq.locus_tag.txt";
 print LOG "htseq-count  -t gene -q -i locus_tag $workdir/$sample/mapping_results/sortedname.paired.prokaryote_ref.sam  $workdir/differential_gene/prokaryote/$tmpdir/prokaryote.gff  > $tmpoutfile\n";
 `htseq-count -t gene -q -i locus_tag $workdir/$sample/mapping_results/sortedname.paired.prokaryote_ref.sam  $workdir/differential_gene/prokaryote/$tmpdir/prokaryote.gff  > $tmpoutfile`;
  }
# unlink "$workdir/$sample/mapping_results/sortedname.paired.prokaryote_ref.sam";
 }

if ( $test eq 'both' || $test eq 'eukarya')
{

if (-e "$workdir/$sample/mapping_results/paired.eukarya_ref.sam" )    
{
`samtools view -bt $workdir/eukarya.fa.fai  $workdir/$sample/mapping_results/paired.eukarya_ref.sam  > $workdir/$sample/mapping_results/paired.eukarya_ref.bam`;
# unlink "$workdir/$sample/mapping_results/paired.eukarya_ref.sam";
}

if (-e "$workdir/$sample/mapping_results/paired.eukarya_ref.bam" )    
{
`samtools sort -n $workdir/$sample/mapping_results/paired.eukarya_ref.bam -o $workdir/$sample/mapping_results/sortedname.paired.eukarya_ref.bam`;
# unlink "$workdir/$sample/mapping_results/paired.eukarya_ref.bam";
}
`samtools view -h $workdir/$sample/mapping_results/sortedname.paired.eukarya_ref.bam > $workdir/$sample/mapping_results/sortedname.paired.eukarya_ref.sam`;

 opendir(DIR, "$workdir/differential_gene/eukarya/") or die $!;
 while (my $tmpdir = readdir(DIR))
  {
    next unless -d "$workdir/differential_gene/eukarya/$tmpdir";
    next if ($tmpdir =~ /^\./);
    my $tmpoutfile="$workdir/sum_gene_count/tmp_count/eukarya/$tmpdir/$sample.eukarya.name.htseq.locus_tag.txt";
  print LOG "htseq-count  -t gene -q -i gene_id  $workdir/$sample/mapping_results/sortedname.paired.eukarya_ref.sam $workdir/differential_gene/eukarya/$tmpdir/eukarya.gff  > $tmpoutfile\n";
  `htseq-count  -t gene  -q -i gene_id  $workdir/$sample/mapping_results/sortedname.paired.eukarya_ref.sam $workdir/differential_gene/eukarya/$tmpdir/eukarya.gff  > $tmpoutfile`;
  }
 #unlink "$workdir/$sample/mapping_results/sortedname.paired.eukarya_ref.sam";
 }

open (LOGFILE, ">$workdir/$sample/mapping_results/finish_htseq-count.txt") or die "  can not open $workdir/$sample/mapping_results/finish_htseq-count.txt $!";
print LOGFILE "$sample is done with htseq-count\n";
