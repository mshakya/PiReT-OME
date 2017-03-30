#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use PiReT::Map;
use Getopt::Long;
use File::Basename;

my ($pairedReadsFile1,  $indexFile,     $pairedReadsFile2,
    $unpairedReadsFile, $workdir,       $euk_ref,
    $hisat2opts,        $sample,        $prok_ref,
    $sampleDir,         $mapDir,        $splicesite,
    $samindex,          $samindex_prok, $samindex_euk
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
        Map::createHisatIndex(
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
    Map::createHisatIndex(
        f1        => $euk_ref,
        numCPU    => $numCPU,
        out_index => $indexFile
    );
}
elsif ( $kingdom eq 'prokaryote' ) {
    if ( -e $indexFile . ".ht2l" ) {
        print "Index already exists";
    }
    Map::createHisatIndex(
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

Map::runMapping(
    r1             => $pairedReadsFile1,
    r2             => $pairedReadsFile2,
    hisat2options  => "--fast",
    IndexFile      => $indexFile,
    numCPU         => $numCPU,
    splicesite     => $splicesite,
    outsam         => "$mapDir/mapped.sam",
    mappingLogFile => "$mapDir/mapped.log"
);

open( my $UMAP1, ">", "$mapDir/unmap_1.fastq" ) or die "Damn. $!";
open( my $UMAP2, ">", "$mapDir/unmap_2.fastq" ) or die "Damn. $!";

if ( $kingdom eq 'both' ) {
    open( my $PR_EUK, ">", "$mapDir/paired_euk.sam" ) or die "Damn. $!";
    open( my $PR_PROK, ">", "$mapDir/paired_prok.sam" )
        or die "Damn. $!";
    open( my $NP_EUK, ">", "$mapDir/notproper_euk.sam" )
        or die "Damn. $!";
    open( my $NP_PROK, ">", "$mapDir/notproper_prok.sam" )
        or die "Damn. $!";
    open( my $FW_PROK, ">", "$mapDir/fwd_prok.sam" )
        or die "Damn. $!";
    open( my $RV_PROK, ">", "$mapDir/rev_prok.sam" )
        or die "Damn. $!";
    open( my $FW_EUK, ">", "$mapDir/fwd_euk.sam" ) or die "Damn. $!";
    open( my $RV_EUK, ">", "$mapDir/rev_euk.sam" ) or die "Damn. $!";
    open( my $FW_EUK_FQ, ">", "$mapDir/fwd_euk_unmapped.fastq" )
        or die "Damn. $!";
    open( my $RV_EUK_FQ, ">", "$mapDir/rev_euk_unmapped.fastq" )
        or die "Damn. $!";
    open( my $FW_PROK_FQ, ">", "$mapDir/fwd_prok_unmapped.fastq" )
        or die "Damn. $!";
    open( my $RV_PROK_FQ, ">", "$mapDir/rev_prok_unmapped.fastq" )
        or die "Damn. $!";
    open( FH, "$mapDir/mapped.sam" ) or die "Damn. $!";

    while (<FH>) {
        chomp;
        next if (/^\@/);
        my $samline = $_;
        my @samFields = split /\t/, $samline;

        Map::parsePairedUnmapped( $samline, $UMAP1, $UMAP2 );

        $samindex_euk = join('.', $euk_ref, 'fai');
        if ( !-e $samindex_euk ){
            Map::createFAI($samindex_euk);
            if ( !-e $samindex_euk) {die "Cannot create $samindex_euk $!"};
        }

        $samindex_prok = join('.', $prok_ref, 'fai');
        if ( !-e $samindex_prok){
            Map::createFAI($samindex_prok);
            if ( !-e $samindex_prok) {die "Cannot create $samindex_prok $!"};
        }

        my %prok_id = Map::parseFAI($samindex_prok);
        my %euk_id  = Map::parseFAI($samindex_euk);
        if ( $euk_id{ $samFields[2] } ) {
            Map::parsePairedMapped(
                samline    => $samline,
                paired_out => $PR_PROK,
                nonproper  => $NP_PROK
            );
            Map::parseSingles(
                samline  => $samline,
                unMapFwd => $FW_EUK_FQ,
                unMapRev => $RV_EUK_FQ,
                MapFwd   => $FW_EUK,
                MapRev   => $RV_EUK
            );
        }
        elsif ( $prok_id{ $samFields[2] } ) {
            Map::parsePairedMapped(
                samline    => $samline,
                paired_out => $PR_EUK,
                nonproper  => $NP_EUK
            );
            Map::parseSingles(
                samline  => $samline,
                unMapFwd => $FW_PROK_FQ,
                unMapRev => $RV_PROK_FQ,
                MapFwd   => $FW_PROK,
                MapRev   => $RV_PROK
            );
        }
    }

}

elsif ( $kingdom eq 'prokaryote' ) {
    open( my $PR_PROK, ">", "$mapDir/paired_prok.sam" )
        or die "Damn. $!";
    open( my $NP_PROK, ">", "$mapDir/notproper_prok.sam" )
        or die "Damn. $!";
    open( my $FW_PROK, ">", "$mapDir/fwd_prok.sam" )
        or die "Damn. $!";
    open( my $RV_PROK, ">", "$mapDir/rev_prok.sam" )
        or die "Damn. $!";
    open( my $FW_PROK_FQ, ">", "$mapDir/fwd_prok_unmapped.fastq" )
        or die "Damn. $!";
    open( my $RV_PROK_FQ, ">", "$mapDir/rev_prok_unmapped.fastq" )
        or die "Damn. $!";
    open( FH, "$mapDir/mapped.sam" ) or die "Damn. $!";

    while (<FH>) {
        chomp;
        next if (/^\@/);
        my $samline = $_;
        Map::parsePairedUnmapped( $samline, $UMAP1, $UMAP2 );
        Map::parsePairedMapped(
            samline    => $samline,
            paired_out => $PR_PROK,
            nonproper  => $NP_PROK
        );
        Map::parseSingles(
            samline  => $samline,
            unMapFwd => $FW_PROK_FQ,
            unMapRev => $RV_PROK_FQ,
            MapFwd   => $FW_PROK,
            MapRev   => $RV_PROK
        );
    }
}
elsif ( $kingdom eq 'eukarya' ) {
    open( my $PR_EUK, ">", "$mapDir/paired_euk.sam" ) or die "Damn. $!";
    open( my $NP_EUK, ">", "$mapDir/notproper_euk.sam" )
        or die "Damn. $!";
    open( my $FW_EUK, ">", "$mapDir/fwd_prok.sam" ) or die "Damn. $!";
    open( my $RV_EUK, ">", "$mapDir/rev_euk.sam" )  or die "Damn. $!";
    open( my $FW_EUK_FQ, ">", "$mapDir/fwd_euk_unmapped.fastq" )
        or die "Damn. $!";
    open( my $RV_EUK_FQ, ">", "$mapDir/rev_euk_unmapped.fastq" )
        or die "Damn. $!";
    open( FH, "$mapDir/mapped.sam" ) or die "Damn. $!";
    while (<FH>) {
        chomp;
        next if (/^\@/);
        my $samline = $_;
        Map::parsePairedUnmapped( $samline, $UMAP1, $UMAP2 );
        Map::parsePairedMapped(
            samline    => $samline,
            paired_out => $PR_EUK,
            nonproper  => $NP_EUK
        );
        Map::parseSingles(
            samline  => $samline,
            unMapFwd => $FW_EUK_FQ,
            unMapRev => $RV_EUK_FQ,
            MapFwd   => $FW_EUK,
            MapRev   => $RV_EUK
        );
    }
}

sub printVersion {
    print basename($0), " version: $version\n";
    exit;
}
