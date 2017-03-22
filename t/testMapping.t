#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";

# use Map::Mapping;
# use Test::Simple tests =>1;
use Test::More;

# use Test::File;

# Verify module can be included via "use" pragma
BEGIN { use_ok('Map::Mapping') }

# Verify module can be included via "require" pragma
require_ok('Map::Mapping');

# Verify if mapping produces an expected sam file
Mapping::runMapping(
    r1             => "data/1.trimmed.fastq",
    r2             => "data/2.trimmed.fastq",
    hisat2options  => "--fast",
    IndexFile      => "data/euk_prok_index",
    numCPU         => 1,
    splicesite     => "data/splice_sites_gff.txt",
    outsam         => "results/mapped.sam",
    mappingLogFile => "results/mapped.log"
);

my $map_cnt = count_lines("results/mapped.sam");
is( $map_cnt, 62, "runMapping() IS test" );

# Verify if reads that did not map are being parsed correctly
open( my $unmap1, ">", "results/unmap_1.fastq" ) or die "Damn. $!";
open( my $unmap2, ">", "results/unmap_2.fastq" ) or die "Damn. $!";
open( FH, "results/mapped.sam" ) or die "Damn. $!";
while (<FH>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Mapping::parsePairedUnmapped( $samline, $unmap1, $unmap2 );
}
close $unmap1;
close $unmap2;
close FH;

my $un_cnt1 = count_lines("results/unmap_1.fastq");
my $un_cnt2 = count_lines("results/unmap_2.fastq");

is( $un_cnt1, 76, "parsePairedUnmapped() read 1 test" );
is( $un_cnt2, 76, "parsePairedUnmapped() read 2 test" );

# Verify if reads are properly mapped ---> <----
open( my $Pmap, ">", "results/paired.mapped.sam" ) or die "Damn. $!";
open( FH1, "results/mapped.sam" ) or die "Damn. $!";
while (<FH1>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Mapping::parsePairedMapped(
        samline    => $samline,
        paired_out => $Pmap
    );
}
close $Pmap;
close FH1;

my $pd_cnt = count_lines("results/paired.mapped.sam");
is( $pd_cnt, 12, "parsePairedmapped() test" );

# Verify if mapped single reads (Forward/Reverse) are correctly parsed
open( my $FR, ">", "results/forward.unmapped.fastq" ) or die "Damn. $!";
open( my $RR, ">", "results/reverse.unmapped.fastq" ) or die "Damn. $!";
open( FH1, "results/mapped.sam" ) or die "Damn. $!";
while (<FH1>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Mapping::parseSingleUnMapped(samline  => $samline,
                                 unMapFwd => $FR,
                                 unMapRev => $RR
    );
}
close $FR;
close $RR;
close FH1;

my $fwd_cnt = count_lines("results/forward.unmapped.fastq");
my $rev_cnt = count_lines("results/reverse.unmapped.fastq");
is( $fwd_cnt, 4, "parseSingleUnmapped() test" );
is( $rev_cnt, 4, "parseSingleUnmapped() test" );

# Function to count lines
sub count_lines {
    my $fn = shift;

    my $cnt;
    open( FH, $fn ) or die "Damn. $!";
    $cnt++ while <FH>;
    close FH;
    return $cnt;
}

done_testing();
