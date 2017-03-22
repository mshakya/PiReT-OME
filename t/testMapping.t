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

# Verify if indexes are created correctly
Mapping::createHisatIndex(
    f1 => "data/euk_test.fna",
    f2 => "data/prok_test.fna",
    numCPU => 1,
    out_index => "results/euk_prok_index");

my $in_cnt = count_lines("results/euk_prok_index.6.ht2l");
is( $in_cnt, 1852, "createHisatIndex() IS test");

# Verify if mapping produces an expected sam file
Mapping::runMapping(
    r1             => "data/1.trimmed.fastq",
    r2             => "data/2.trimmed.fastq",
    hisat2options  => "--fast",
    IndexFile      => "results/euk_prok_index",
    numCPU         => 1,
    splicesite     => "data/splice_sites_gff.txt",
    outsam         => "results/mapped.sam",
    mappingLogFile => "results/mapped.log"
);

my $map_cnt = count_lines("results/mapped.sam");
is( $map_cnt, 69, "runMapping() IS test" );

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

is( $un_cnt1, 72, "parsePairedUnmapped() read 1 test" );
is( $un_cnt2, 72, "parsePairedUnmapped() read 2 test" );

# Verify if reads are properly mapped ---> <----
open( my $Pmap, ">", "results/paired.mapped.sam" ) or die "Damn. $!";
open( my $Nmap, ">", "results/Nonproper.mapped.sam" ) or die "Damn. $!";
open( FH1, "results/mapped.sam" ) or die "Damn. $!";
while (<FH1>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Mapping::parsePairedMapped(
        samline    => $samline,
        paired_out => $Pmap,
        nonproper => $Nmap
    );
}
close $Pmap;
close $Nmap;
close FH1;

my $pd_cnt = count_lines("results/paired.mapped.sam");
is( $pd_cnt, 8, "parsePairedmapped() test" );

# Verify if single reads (mapped/unmaped Forward/Reverse) are correctly parsed
open( my $FR, ">", "results/forward.unmapped.fastq" ) or die "Damn. $!";
open( my $RR, ">", "results/reverse.unmapped.fastq" ) or die "Damn. $!";
open( my $FM, ">", "results/forward.mapped.sam" )     or die "Damn. $!";
open( my $RM, ">", "results/reverse.mapped.sam" )     or die "Damn. $!";
open( FH1, "results/mapped.sam" ) or die "Damn. $!";
while (<FH1>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Mapping::parseSingles(
        samline  => $samline,
        unMapFwd => $FR,
        unMapRev => $RR,
        MapFwd   => $FM,
        MapRev   => $RM
    );
}
close $FR;
close $RR;
close $FM;
close $RM;
close FH1;

my $fwdU_cnt = count_lines("results/forward.unmapped.fastq");
my $revU_cnt = count_lines("results/reverse.unmapped.fastq");
my $fwdM_cnt = count_lines("results/forward.mapped.sam");
my $revM_cnt = count_lines("results/reverse.mapped.sam");
is( $fwdU_cnt, 4, "parseSingles() test" );
is( $revU_cnt, 16, "parseSingles() test" );
is( $fwdM_cnt, 4, "parseSingles() test" );
is( $revM_cnt, 1, "parseSingles() test" );


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
