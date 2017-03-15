#!/usr/bin/perl

package Mapping;
use Exporter;
@ISA = ('Exporter');
@EXPORT = ('hello');


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

