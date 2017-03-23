#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Map::Mapping;
use Getopt::Long;
use File::Basename;

# perl $scriptDir/rRNA_reads_mapping.pl-test $test -cpu $numCPU
#-p1  $workdir/$sample/trimming_results/$sample.1.trimmed.fastq
#-p2 $workdir/$sample/trimming_results/$sample.2.trimmed.fastq
#-sample $sample
#-index $indexref
#-o $workdir

my ($pairedReadsFile1,  $indexFile, $pairedReadsFile2,
    $unpairedReadsFile, $workdir,   $euk_ref,
    $hisat2opts,        $sample,    $prok_ref,
    $sampleDir,         $mapDir,    $splicesite
);
my $version = "develop";
my $kingdom = 'prokaryote';
my $numCPU  = 1;

GetOptions(
    'P1|paired_reads1=s' => \$pairedReadsFile1,
    'P2|paired_reads2=s' => \$pairedReadsFile2,
    'U|unpaired_read=s'  => \$unpairedReadsFile,
    'I|index_ref_ht2=s'  => \$indexFile,
    'S|sample=s'         => \$sample,
    'K|kingdom=s'        => \$kingdom,
    'P|numCPU=i'         => \$numCPU,                 # bwa option
    'hisat2opt=s'        => \$hisat2opts,             # bwa mem options
    'W|workdir=s'        => \$workdir,
    'E|euk_ref=s'        => \$euk_ref,
    'B|prok_ref=s'       => \$prok_ref,
    'T|splicesite=s'     => \$splicesite,
    'V|version'          => sub { printVersion() },
    'help|h|?'           => sub { &Usage() }
);

sub Usage {
    print <<"END";

DESCRIPTION\n
    Maps (using hisat2) reads to given reference index (creates one if its not
    present)
 
Usage: 
                    
OPTIONS
   Generic Program Information
        -h, --help 
            Print a usage message briefly summarizing these command-line
            options, then exits.

        -V, --version
            Print the version number to the standard output stream. 
            This version number should be included in all bug reports.

        -P1, --paired_reads1
            Fastq file containing first read of the pair.

        -P2, --paired_reads2
            Fastq file contatining second read of the pair.

        -U, --unpaired_read
            Fastq file of read that are not paired.

        -W, --workdir
            Path to a directory where the whole analysis will be stored. 
            Must have write permission.

        -I, --index_ref_ht2
            Hisat2 mapping index file. 
            If the mapping index is already generated, one can pass that as well.

        -K, --kingdom
            `both` (for eukarya and prokaryote), `prokaryote`, or `eukarya`. 
            Default: prokaryote

        -E, --euk_ref
            Path to eukarya genome reference fasta

        -B, --prok_ref
            Path to bacterial genome reference fasta

        -P, --numCPU
            number of cpu to be used.
            Default: 1

        -S, --sample
            name of the sample

        -T, --splicesite
            tab delimited file with information on splice sites

END
    exit;
}

unless ( $indexFile
    && $workdir
    && $kingdom
    && ( $euk_ref || $prok_ref )
    && ( ( $pairedReadsFile1 & $pairedReadsFile2 ) || $unpairedReadsFile ) )
{
    &Usage;
}

if ( $kingdom eq 'both' ) {
    if ( -e $indexFile . ".ht2l" ) {
        print "Index already exists";
    }
    else {
        Mapping::createHisatIndex(
            f1        => $euk_ref,
            f2        => $prok_ref,
            numCPU    => $numCPU,
            out_index => $indexFile
        );
    }
}
elsif ( $kingdom eq 'eukarya' ) {
    if ( -e $indexFile . ".ht2l" ) {
        print "Index already exists";
    }
    Mapping::createHisatIndex(
        f1        => $euk_ref,
        numCPU    => $numCPU,
        out_index => $indexFile
    );
}
elsif ( $kingdom eq 'prokaryote' ) {
    if ( -e $indexFile . ".ht2l" ) {
        print "Index already exists";
    }
    Mapping::createHisatIndex(
        f1        => $prok_ref,
        numCPU    => $numCPU,
        out_index => $indexFile
    );
}
else {
    print "Error!, not the correct option for kingdom"
        . "The correct options are: both, eukarya, and prokaryote";
}

# First, check if the working directory exists
if ( !-e $workdir ) { print "no working dir $workdir found\n"; }

# Second, make directories based on sample name within working dir
$sampleDir = join '/', ( $workdir, $sample );
mkdir $sampleDir if ( !-e $sampleDir );

# Third, throw, a complain, if cant make directory
if ( !-e $sampleDir ) { print "cannot make dir $sampleDir\n"; }

# Fourth, make directory for adding the mapping results
$mapDir = join '/', ( $sampleDir, 'mapping_results' );
mkdir $mapDir if ( !-e $mapDir );
if ( !-e $mapDir ) { print "cannot make dir $mapDir\n"; }

# Mapping::open_samfiles( $kingdom, $mapDir );

Mapping::runMapping(
    r1             => $pairedReadsFile1,
    r2             => $pairedReadsFile2,
    hisat2options  => "--fast",
    IndexFile      => $indexFile,
    numCPU         => $numCPU,
    splicesite     => $splicesite,
    outsam         => "$mapDir/mapped.sam",
    mappingLogFile => "$mapDir/mapped.log"
);

open( UMAP1, ">", "$mapDir/unmap_1.fastq" ) or die "Damn. $!";
open( UMAP2, ">", "$mapDir/unmap_2.fastq" ) or die "Damn. $!";

if ( $kingdom eq 'both' ) {
    open( PR_EUK, ">", "$mapDir/paired_euk.sam" ) or die "Damn. $!";
    open( PR_PROK, ">", "$mapDir/paired_prok.sam" )
        or die "Damn. $!";
    open( NP_EUK, ">", "$mapDir/notproper_euk.sam" )
        or die "Damn. $!";
    open( NP_PROK, ">", "$mapDir/notproper_prok.sam" )
        or die "Damn. $!";
    open( FW_PROK, ">", "$mapDir/fwd_prok.sam" )
        or die "Damn. $!";
    open( RV_PROK, ">", "$mapDir/rev_prok.sam" )
        or die "Damn. $!";
    open( FW_EUK, ">", "$mapDir/fwd_euk.sam" ) or die "Damn. $!";
    open( RV_EUK, ">", "$mapDir/rev_euk.sam" ) or die "Damn. $!";
    open( FW_EUK_FQ, ">", "$mapDir/fwd_euk_unmapped.fastq" ) or die "Damn. $!";
    open( RV_EUK_FQ, ">", "$mapDir/rev_euk_unmapped.fastq" ) or die "Damn. $!";
    open( FH, "$mapDir/mapped.sam" ) or die "Damn. $!";
    
    while (<FH>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    my @samFields = split /\t/, $samline;
    Mapping::parsePairedUnmapped( $samline, UMAP1, UMAP2 );
    if ( $samFields[2] eq )
    
    Mapping::parsePairedMapped(
        samline    => $samline,
        paired_out => $Pmap,
        nonproper  => $Nmap
    );
    
    Mapping::parseSingles()
    }
}
elif( $kingdom eq 'prokaryote' ) {
    open( PR_PROK, ">", "$mapDir/paired_prok.sam" )
        or die "Damn. $!";
    open( NP_PROK, ">", "$mapDir/notproper_prok.sam" )
        or die "Damn. $!";
    open( FW_PROK, ">", "$mapDir/fwd_prok.sam" )
        or die "Damn. $!";
    open( RV_PROK, ">", "$mapDir/rev_prok.sam" )
        or die "Damn. $!";
    open( FW_PROK_FQ, ">", "$mapDir/fwd_prok_unmapped.fastq" ) or die "Damn. $!";
    open( RV_PROK_FQ, ">", "$mapDir/rev_prok_unmapped.fastq" ) or die "Damn. $!";
    open( FH, "$mapDir/mapped.sam" ) or die "Damn. $!";
    
    while (<FH>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Mapping::parsePairedUnmapped( $samline, UMAP1, UMAP2 );
    Mapping::parsePairedMapped(
        samline    => $samline,
        paired_out => PR_PROK,
        nonproper  => NP_PROK
    );
    Mapping::parseSingles(
        samline  => $samline,
        unMapFwd => $FR,
        unMapRev => $RR,
        MapFwd   => $FM,
        MapRev   => $RM
    );
}
elif( $kingdom eq 'eukarya' ) {
    open( PR_EUK, ">", "$mapDir/paired_euk.sam" ) or die "Damn. $!";
    open( NP_EUK, ">", "$mapDir/notproper_euk.sam" )
        or die "Damn. $!";
    open( FW_EUK, ">", "$mapDir/fwd_prok.sam" ) or die "Damn. $!";
    open( RV_EUK, ">", "$mapDir/rev_euk.sam" ) or die "Damn. $!";
    open( FW_EUK_FQ, ">", "$mapDir/fwd_euk_unmapped.fastq" ) or die "Damn. $!";
    open( RV_EUK_FQ, ">", "$mapDir/rev_euk_unmapped.fastq" ) or die "Damn. $!";
    Mapping::parsePairedUnmapped( $samline, UMAP1, UMAP2 );
    Mapping::parsePairedMapped(
        samline    => $samline,
        paired_out => PR_EUK,
        nonproper  => NP_EUK
}

open( FH, "$mapDir/mapped.sam" ) or die "Damn. $!";
while (<FH>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Mapping::parsePairedUnmapped( $samline, $unmap1, $unmap2 );
    Mapping::parsePairedMapped(
        samline    => $samline,
        paired_out => $Pmap,
        nonproper  => $Nmap
    );
}
close $Nmap;
close $Pmap;
close $unmap1;
close $unmap2;
close FH;

sub printVersion {
    print basename($0), " version: $version\n";
    exit;
}
