#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../ext/lib/perl5";
use Test::More;

# use Test::File;

# Verify module can be included via "use" pragma
BEGIN { use_ok('PiReT::Map') }

# Verify module can be included via "require" pragma
require_ok('PiReT::Map');

# Verify if indexes are created correctly
Map::createHisatIndex(
    f1        => "t/data/euk_test.fna",
    f2        => "t/data/prok_test.fna",
    numCPU    => 1,
    out_index => "t/results/euk_prok_index"
);

my $in_cnt = count_lines("t/results/euk_prok_index.6.ht2l");
is( $in_cnt, 2195, "createHisatIndex() IS test" );

# Verify if mapping produces an expected sam file
Map::runMapping(
    r1             => "t/data/1.trimmed.fastq",
    r2             => "t/data/2.trimmed.fastq",
    hisat2options  => "--fast",
    IndexFile      => "t/results/euk_prok_index",
    numCPU         => 1,
    splicesite     => "t/data/splice_sites_gff.txt",
    outsam         => "t/results/mapped.sam",
    mappingLogFile => "t/results/mapped.log"
);

my $map_cnt = count_lines("t/results/mapped.sam");
is( $map_cnt, 71, "runMapping() IS test" );

# Verify if reads that did not map are being parsed correctly
open( my $unmap1, ">", "t/results/unmap_1.fastq" ) or die "Damn. $!";
open( my $unmap2, ">", "t/results/unmap_2.fastq" ) or die "Damn. $!";
open( FH, "t/results/mapped.sam" ) or die "Damn. $!";
while (<FH>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Map::parsePairedUnmapped( $samline, $unmap1, $unmap2 );
}
close $unmap1;
close $unmap2;
close FH;

my $un_cnt1 = count_lines("t/results/unmap_1.fastq");
my $un_cnt2 = count_lines("t/results/unmap_2.fastq");

is( $un_cnt1, 56, "parsePairedUnmapped() read 1 test" );
is( $un_cnt2, 56, "parsePairedUnmapped() read 2 test" );

# Verify if reads are properly mapped ---> <----
open( my $Pmap, ">", "t/results/paired.mapped.sam" )    or die "Damn. $!";
open( my $Nmap, ">", "t/results/Nonproper.mapped.sam" ) or die "Damn. $!";
open( FH1, "t/results/mapped.sam" ) or die "Damn. $!";
while (<FH1>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Map::parsePairedMapped(
        samline    => $samline,
        paired_out => $Pmap,
        nonproper  => $Nmap
    );
}
close $Pmap;
close $Nmap;
close FH1;

my $pd_cnt = count_lines("t/results/paired.mapped.sam");
is( $pd_cnt, 14, "parsePairedmapped() test" );

# Verify if single reads (mapped/unmaped Forward/Reverse) are correctly parsed
open( my $FR, ">", "t/results/forward.unmapped.fastq" ) or die "Damn. $!";
open( my $RR, ">", "t/results/reverse.unmapped.fastq" ) or die "Damn. $!";
open( my $FM, ">", "t/results/forward.mapped.sam" )     or die "Damn. $!";
open( my $RM, ">", "t/results/reverse.mapped.sam" )     or die "Damn. $!";
open( FH1, "t/results/mapped.sam" ) or die "Damn. $!";
while (<FH1>) {
    chomp;
    next if (/^\@/);
    my $samline = $_;
    Map::parseSingles(
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

my $fwdU_cnt = &count_lines("t/results/forward.unmapped.fastq");
my $revU_cnt = &count_lines("t/results/reverse.unmapped.fastq");
my $fwdM_cnt = &count_lines("t/results/forward.mapped.sam");
my $revM_cnt = &count_lines("t/results/reverse.mapped.sam");
is( $fwdU_cnt, 4,  "parseSingles() test" );
is( $revU_cnt, 16, "parseSingles() test" );
is( $fwdM_cnt, 4,  "parseSingles() test" );
is( $revM_cnt, 1,  "parseSingles() test" );

# Verify if counting function works
Map::sumMaps(
    fastq => "t/data/1.trimmed.fastq",
    ppm   => "t/results/paired.mapped.sam",
    npm   => "t/results/Nonproper.mapped.sam",
    fm    => "t/results/forward.mapped.sam",
    rm    => "t/results/reverse.mapped.sam",
    um    => "t/results/unmap_1.fastq",
    out   => "t/results/stats_table.tab"
);
my $stat_cnt = &count_lines("t/results/stats_table.tab");
is( $stat_cnt, 8, "sumMaps() IS test" );

# Verify if the parseFAI function is working as it should
my %fai_dic = Map::parseFAI("t/data/prokaryote.fa.fai");
my $ind     = $fai_dic{"gi|50196905|ref|NC_007530.2|"};
is( $ind, '5227419', "parseFAI() IS test" );

#Verify if thr createFAI is working
Map::createFAI("t/data/euk_test.fna\n");
Map::createFAI("t/data/prok_test.fna\n");
my $fai_euk  = &count_lines("t/data/euk_test.fna.fai");
my $fai_prok = &count_lines("t/data/prok_test.fna.fai");
is( $fai_euk,  8, "createFAI() IS test with euk\n" );
is( $fai_prok, 2, "createFAI() IS test with prok" );

# Verify if parseMapProk is working
Map::executeCommand("mkdir -p t/results/prok_map_test\n");
Map::parseMapProk(
    mapDir  => "t/results/prok_map_test",
    samFile => "t/results/mapped.sam"
);
my $pp_prok = &count_lines("t/results/prok_map_test/paired_prok.sam");
is( $pp_prok, 14, "parseMapProk() IS test" );

# Verify if parseMapEuk is working
Map::executeCommand("mkdir -p t/results/euk_map_test\n");
Map::parseMapEuk(
    mapDir  => "t/results/euk_map_test",
    samFile => "t/results/mapped.sam"
);
my $pp_euk = &count_lines("t/results/euk_map_test/paired_euk.sam");
is( $pp_euk, 14, "parseMapEuk() IS test" );

# Verify if parseMapBoth is working
Map::executeCommand("mkdir -p t/results/both_map_test\n");
Map::parseMapEuk(
    mapDir  => "t/results/both_map_test",
    samFile => "t/results/mapped.sam",
    ref_prok => "t/data/prok_test.fna",
    ref_euk=> "t/data/euk_test.fna"
);

my $pp_euk_both = &count_lines("t/results/both_map_test/paired_euk.sam");
is( $pp_euk_both, 14, "parseMapBoth() IS test" );


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
